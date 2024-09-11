#!/usr/bin/env python3
"""Compare images in two folders."""

import argparse
import enum
import json
import re
import subprocess
import tempfile
import textwrap

from pathlib import Path

import colors


class MissingVariable(Exception):
    """Simple class to declare that variables not defined by a datafile."""

    def __init__(self, variable: str, datafile: str):
        self.variable = variable
        self.datafile = datafile

    def __str__(self):
        return f"Missing from {self.datafile}: {self.variable}"


class ExecutionError(Exception):
    """Simple class to indicate an error during execution."""

    def __init__(self, error: str):
        self.error = error

    def __str__(self):
        return "Error output:\n" + self.error


class MismatchedPlot(Exception):
    """Indicate that a plot changes when the data source is changed."""

    def __init__(self, image_name):
        self.image_name = image_name

    def __str__(self):
        return f"Full and restricted plots of {self.image_name} don't match"


class ValidationError(Exception):
    """Simple class to indicate a validation error."""

    def __str__(self):
        original = super().__str__()
        lines = original.splitlines()
        for index, line in enumerate(lines):
            if index == 0:
                lines[index] = colors.yellow(line)
            else:
                lines[index] = colors.faint(line)

        return "\n".join(lines)


class Usage(enum.IntEnum):
    """Enumeration to detail how a variable from a datafile is used."""

    ACCESSED = 1
    REDUNDANT = 2
    DIFFERENT = 3


def simplify_symbol(symbol: str) -> str:
    """Simplify an R symbol."""
    if re.match(r"^[\w|\.]+$", symbol):
        return symbol

    if match := re.match(r"^([\w|\.]+)\$", symbol):
        return simplify_symbol(match.group(1))

    if match := re.match(r"^([\w|\.]+)\[", symbol):
        return simplify_symbol(match.group(1))

    if match := re.match(r"^(?:names|colnames|rownames)\(([\w|\.]+)\)$", symbol):
        return simplify_symbol(match.group(1))

    return symbol


def make_mosaic(restricted: Path, output: Path):
    """Make a mosaic of the restricted, full, and reference images."""
    reference = restricted.parent.parent / "reference_figures" / restricted.name

    assert reference.is_file()

    ref_dims = subprocess.check_output([
        "magick",
        "identify",
        "-format",
        "%h:%w",
        reference,
    ]).decode("utf-8")
    ref_height, ref_width = (int(item) for item in ref_dims.split(":"))

    if restricted.is_file():
        restricted_args = [restricted, "-resize", f"x{ref_height}"]
    else:
        restricted_args = [
            "(",
            "-size",
            f"{ref_width}x{ref_height}",
            "xc:gray",
            ")",
        ]

    subprocess.check_call([
        "magick",
        *restricted_args,
        reference,
        "+append",
        output,
    ])


class Figure:
    """Class to compare Figure usage."""

    @classmethod
    def from_logs(cls, logfile: Path) -> list["Figure"]:
        """Parse Figures from the logfile."""
        log_text = logfile.read_text(encoding="utf-8")
        # Just get the logs from the last session
        last_text = log_text.split("Starting up\n")[-1]

        symbol_re = re.compile(
            r"""
            (?P<action>(?:Redundant|Accessed|Different))
            \|
            (?P<figure>[^|]+)
            \|
            file:(?P<source>[^|]+)
            \|
            (?P<variable>[^|]+)
            $""",
            flags=re.VERBOSE,
        )

        figure_start_re = re.compile(r"^.*33;1m(Figure[^.]+\.R)")
        source_start_re = re.compile(r"^.*32;1mSourcing ([^.]+\.\w+)")
        problem_re = re.compile(r"^.*31;1mProblem with (Figure[^.]+\.R) !")

        current_figure_lines = []

        figures = {}

        active_figure = ""
        active_source = ""

        for line in last_text.splitlines():
            if source_match := source_start_re.search(line.strip()):
                active_source = source_match.group(1)
                continue

            if figure_match := figure_start_re.search(line.strip()):
                active_figure = figure_match.group(1)
                current_figure_lines = []

                if active_figure not in figures:
                    figures[active_figure] = cls(Path(active_figure))

                # Determine if the source is a full or restricted dataset
                if "data" in Path(active_source).parts:
                    if not figures[active_figure].restricted_dataset:
                        figures[active_figure].restricted_dataset = active_source
                    assert figures[active_figure].restricted_dataset == active_source

                if "outlier" in Path(active_source).parts:
                    if not figures[active_figure].full_dataset:
                        figures[active_figure].full_dataset = active_source
                    assert figures[active_figure].full_dataset == active_source
                continue

            if problem_match := problem_re.search(line.strip()):
                assert active_figure
                assert active_source

                assert problem_match.group(1) == active_figure
                if active_figure not in figures:
                    figures[active_figure] = cls(Path(active_figure))

                assert active_source not in figures[active_figure].errors
                figures[active_figure].errors[active_source] = current_figure_lines[:]
                continue

            if symbol := symbol_re.search(line.strip()):
                action, figurefile, source, variable = symbol.groups()
                if figurefile not in figures:
                    figures[figurefile] = cls(Path(figurefile))

                source_dict = figures[figurefile].variables.setdefault(source, {})
                source_dict[variable] = max(
                    Usage[action.upper()], source_dict.get(variable, Usage.ACCESSED)
                )

                # Determine if the source is a full or restricted dataset
                if "data" in Path(source).parts:
                    if not figures[figurefile].restricted_dataset:
                        figures[figurefile].restricted_dataset = source
                    assert figures[figurefile].restricted_dataset == source

                if "outlier" in Path(source).parts:
                    if not figures[figurefile].full_dataset:
                        figures[figurefile].full_dataset = source
                    assert figures[figurefile].full_dataset == source
            else:
                current_figure_lines.append(line)

        return list(figures.values())

    def __init__(self, sourcefile: Path):
        self.sourcefile = sourcefile

        self.variables: dict[str, dict[str, Usage]] = {}
        self.restricted_dataset = ""

        self.errors: dict[str, list[str]] = {}

        if not (name_match := re.match(r"Figure(\d+)(\w+)", sourcefile.stem)):
            raise ValueError("Invalid filename!")

        number, letter = name_match.groups()

        self.restricted_base = sourcefile.parent.parent.parent / "restricted_figures"
        self.comparison_base = sourcefile.parent.parent.parent / "comparison_figures"
        template = f"Figure_{number}_{letter}*.png"

        self.expected_images = set(
            path.name for path in self.restricted_base.glob(template)
        )

        if not self.expected_images:
            self.expected_images.add(f"Figure_{number}_{letter}.png")

        # Scrape the defined symbols from the script
        self.defined_symbols = set()

        with self.sourcefile.open(encoding="utf-8") as infile:
            for line in infile:
                if match := re.match(r"^\s*([^#\n]+?)\s*<-", line):
                    self.defined_symbols.add(simplify_symbol(match.group(1).strip()))

    def validate_restricted(self):
        """Validate the restricted dataset."""
        # Second level - were the data sources parsed?
        if not self.restricted_dataset:
            raise ValidationError("Restricted dataset not parsed")

        for imagename in sorted(self.expected_images):
            restricted = self.restricted_base / imagename

            if not restricted.is_file():
                raise ValidationError(f"Restricted {imagename} doesn't exist")

            make_mosaic(
                restricted,
                self.comparison_base / ("comp-" + restricted.name)
            )

        if self.restricted_dataset in self.errors:
            error_lines = "\n".join(self.errors[self.restricted_dataset])
            if undefined_match := re.search(r"object '([^']+)' not found", error_lines):
                raise MissingVariable(undefined_match.group(1), self.restricted_dataset)

            raise ExecutionError(error_lines)

        # Third level - are there modified variables from the restricted
        # dataset?
        bad_symbols = [
            (variable, usage)
            for (variable, usage) in self.variables[self.restricted_dataset].items()
            if usage >= Usage.REDUNDANT
        ]

        bad_symbols.sort(key=lambda x: x[1])

        if bad_symbols:
            raise ValidationError(
                "Variables misused in restricted dataset: %s" % bad_symbols
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile", type=Path)

    args = parser.parse_args()

    missing_variables: dict[str, dict[str, str]] = {}
    execution_errors: dict[str, dict[str, str]] = {}
    mismatched_images = []

    for figure in Figure.from_logs(args.logfile):
        print(colors.bold(f"{figure.sourcefile}:"))
        print("  Restricted dataset:\t", end="")
        try:
            figure.validate_restricted()
            print(colors.green("GOOD"))

        except MissingVariable as err:
            print(colors.yellow(err))
            datadict = missing_variables.setdefault(err.datafile, {})
            datadict[figure.sourcefile.name] = err.variable

        except ExecutionError as err:
            figuredict = execution_errors.setdefault(figure.sourcefile.name, {})
            figuredict[figure.restricted_dataset] = err.error
            print(textwrap.indent(colors.faint(err), "\t\t\t").lstrip())

        except ValidationError as err:
            print(err)

        print()

    with open("missing_variables.json", mode="w", encoding="utf-8") as outfile:
        json.dump(missing_variables, outfile, indent=2)

    with open("execution_errors.json", mode="w", encoding="utf-8") as outfile:
        json.dump(execution_errors, outfile, indent=2)
