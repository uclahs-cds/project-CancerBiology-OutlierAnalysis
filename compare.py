#!/usr/bin/env python3
"""Compare images in two folders."""

import argparse
import subprocess
import re
import enum

from pathlib import Path


class ValidationError(Exception):
    """Simple class to indicate a validation error."""


class Usage(enum.IntEnum):
    """Enumeration to detail how a variable from a datafile is used."""

    ACCESSED = 1
    REDUNDANT = 2
    DIFFERENT = 3


class Figure:
    """Class to compare Figure usage."""

    @classmethod
    def from_logs(cls, logfile: Path) -> list["Figure"]:
        """Parse Figures from the logfile."""
        log_text = logfile.read_text(encoding="utf-8")
        # Just get the logs from the last session
        last_text = log_text.split(" Starting up\n")[-1]

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

        figures = {}

        for line in last_text.splitlines():
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

        return list(figures.values())

    def __init__(self, sourcefile: Path):
        self.sourcefile = sourcefile

        self.variables: dict[str, dict[str, Usage]] = {}
        self.full_dataset = ""
        self.restricted_dataset = ""

        if not (name_match := re.match(r"Figure(\d+)(\w+)", sourcefile.stem)):
            raise ValueError("Invalid filename!")

        number, letter = name_match.groups()

        self.full_base = sourcefile.parent.parent.parent / "full_figures"
        self.restricted_base = sourcefile.parent.parent.parent / "restricted_figures"
        template = f"Figure_{number}_{letter}*.png"

        self.expected_images = set(path.name for path in self.full_base.glob(template))
        self.expected_images |= set(
            path.name for path in self.restricted_base.glob(template)
        )

        if not self.expected_images:
            self.expected_images.add(f"Figure_{number}_{letter}.png")

    def validate(self):
        """Compare the full and restricted images."""
        # First level - do the images match?
        for imagename in sorted(self.expected_images):
            restricted = self.restricted_base / imagename
            full = self.full_base / imagename

            if not restricted.is_file():
                if not full.is_file():
                    raise ValidationError(
                        f"Neither restricted nor full {imagename} exist"
                    )
                raise ValidationError(f"Restricted {imagename} doesn't exist")

            if not full.is_file():
                raise ValidationError(f"Full {imagename} doesn't exist")

            pixel_difference = int(
                subprocess.run(
                    ["compare", "-metric", "AE", restricted, full, "NULL:"],
                    capture_output=True,
                    check=False,
                )
                .stderr.decode("utf-8")
                .strip()
            )

            if pixel_difference:
                raise ValidationError(f"Full and restricted {imagename} don't match")

        # Second level - were the data sources parsed?
        if not self.restricted_dataset:
            raise ValidationError("Restricted dataset not parsed")

        if not self.full_dataset:
            raise ValidationError("Full dataset not parsed")

        # Third level - are there redundant or modified variables from the
        # restricted dataset?
        for variable, usage in self.variables[self.restricted_dataset].items():
            if usage > Usage.ACCESSED:
                raise ValidationError(
                    f"Variable `{variable}` {usage.name} in restricted dataset"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("logfile", type=Path)

    args = parser.parse_args()
    # compare_folders(args.first, args.second)
    for figure in Figure.from_logs(args.logfile):
        print(f"{figure.sourcefile}:\t", end="")
        try:
            figure.validate()
            print("GOOD")
        except ValidationError as err:
            print(err)
