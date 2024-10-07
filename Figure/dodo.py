#!/usr/bin/env python3
"""DoIt dodo file."""

import datetime
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import colors
import doit
import requests

from doit.reporter import ConsoleReporter
from doit.exceptions import TaskFailed, UnmetDependency


class CustomJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, set):
            return list(o)
        if isinstance(o, Path):
            return str(o)
        return super().default(o)


class MissingVariableException(TaskFailed):
    """Indicate that a variable was missing."""


class FailedFigureException(TaskFailed):
    """Indicate that at least one figure failed."""

    def __init__(self, goodfigs, badfigs):
        super().__init__("Figures did not compare identical")
        self.goodfigs = goodfigs
        self.badfigs = badfigs


class MyReporter(ConsoleReporter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.explicit_tasks = []
        self.unmet_map = {}
        self.all_tasks = {}
        self.task_statuses = {}
        self.task_starts = {}

        # Find the last report
        self.prior_duration_data = {}
        self.current_duration_data = {}

        prior_datafiles = sorted(Path(__file__).parent.glob("2024*.json"))
        if prior_datafiles:
            with prior_datafiles[-1].open(encoding="utf-8") as infile:
                prior_data = json.load(infile)

            for name, task in prior_data["tasks"].items():
                duration = task["values"].get("duration", None)
                if duration is not None:
                    self.prior_duration_data[name] = float(duration)

    def initialize(self, tasks, selected_tasks):
        self.all_tasks = tasks
        self.explicit_tasks = selected_tasks[:]

    def _save_duration(self, task):
        try:
            requests.post(
                "http://localhost:65035",
                json={"name": task.title()},
                timeout=0.2
            )

        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            pass
        if task.title() in self.task_starts:
            duration = (datetime.datetime.now() - self.task_starts.pop(task.title())).total_seconds()
            task.values["duration"] = duration
            self.current_duration_data[task.title()] = duration

    def complete_run(self):
        super().complete_run()

        result = {
            "status": self.task_statuses,
            "tasks": {
                name: task.pickle_safe_dict() for (name, task) in self.all_tasks.items()
            },
        }
        for task in result["tasks"].values():
            task.pop("io")

        if self.task_starts:
            print("LEFTOVER START TIMES:", self.task_starts)

        with Path(datetime.datetime.now().isoformat() + ".json").open(
            mode="w", encoding="utf-8"
        ) as outfile:
            json.dump(result, outfile, cls=CustomJSONEncoder, indent=2)

        for task_name in sorted(self.prior_duration_data.keys() & self.current_duration_data.keys()):
            prior_duration = self.prior_duration_data[task_name]
            new_duration = self.current_duration_data[task_name]
            if max(prior_duration, new_duration) > 5:
                print(task_name, prior_duration, new_duration)


    def execute_task(self, task):
        self.task_starts[task.title()] = datetime.datetime.now()

        threshold = 3
        prior_duration = self.prior_duration_data.get(task.title(), 0)

        try:
            requests.post(
                "http://localhost:65035",
                json={"name": task.title(), "seconds": prior_duration},
                timeout=0.2
            )
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            pass

        # if task.actions and not task.title().startswith("summary"):
        if task.actions:
            if prior_duration == 0:
                self.outstream.write(f" -> {task.title()}\n")
            elif threshold <= prior_duration:
                self.outstream.write(f" -> {task.title()} (expect {prior_duration:0.1f} seconds)\n")

    def skip_uptodate(self, task):
        self.task_statuses[task.title()] = "up-to-date"
        if any(task.title().startswith(item) for item in self.explicit_tasks):
            self.outstream.write(colors.faint(f" -- {task.title()}") + "\n")

    def add_success(self, task):
        self.task_statuses[task.title()] = "success"
        self._save_duration(task)
        if task.actions:
            self.outstream.write(colors.green(f"    {task.title()}") + "\n")

    def add_failure(self, task, fail):
        """called when execution finishes with a failure"""
        self.task_statuses[task.title()] = "failure"
        self._save_duration(task)
        super().add_failure(task, fail)

    def skip_ignore(self, task):
        """skipped ignored task"""
        self.task_statuses[task.title()] = "ignore"

    def _write_failure(self, result, write_exception=True):
        task = result["task"]
        exception = result["exception"]

        if isinstance(exception, MissingVariableException):
            self.outstream.write(colors.red(f" xx {task.title()}"))
            self.outstream.write(f" - {exception.message}\n")

        elif isinstance(exception, FailedFigureException):
            self.outstream.write(
                colors.red(f" // {task.title()}") + " - Mismatched images\n"
            )
        elif isinstance(exception, UnmetDependency):
            self.unmet_map[task.title()] = exception.message.split()
            if any(task.title().startswith(item) for item in self.explicit_tasks):
                self.outstream.write(
                    colors.yellow(f" // {task.title()}") + " - Unmet dependencies\n"
                )
        else:
            self.outstream.write(colors.red(f" xx {task.title()}") + " - FAILED\n")
            self.outstream.write(f"{exception.message}\n")
            print(exception)


DOIT_CONFIG = {
    "reporter": MyReporter,
    "verbosity": 2,
    "continue": True,
    "process": 4,
    "default_tasks": ["summary", "digest"],
}


def run_r(task):
    """Run R and set the appropriate values."""
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as instrumented_file:
        with open(task.meta["scriptfile"], "r", encoding="utf-8") as infile:
            for index, line in enumerate(infile):
                if line.startswith("#"):
                    instrumented_file.write(
                        f"message(sprintf('Time: %s | Line: %d', Sys.time(), {index})); "
                    )
                instrumented_file.write(line)

        instrumented_file.flush()

        proc = subprocess.run(
            # ["Rscript", task.meta["scriptfile"]],
            ["Rscript", instrumented_file.name],
            shell=False,
            capture_output=True,
            check=False,
        )

    task.values["returncode"] = proc.returncode
    task.values["stdout"] = proc.stdout.decode("utf-8")
    task.values["stderr"] = proc.stderr.decode("utf-8")

    if proc.returncode != 0:
        if undefined_match := re.search(
            r"object '([^']+)' not found", proc.stderr.decode("utf-8")
        ):
            return MissingVariableException(undefined_match.group(1))
        return TaskFailed(proc.stderr.decode("utf-8"))

    return True


def compare_action(task):  # input_image, output_image, diff_image):
    """Compare the two images."""
    input_image = task.meta["input"]
    output_image = task.meta["output"]
    diff_image = task.meta["diff"]

    with tempfile.TemporaryDirectory() as tempdir:
        width, height = subprocess.check_output(
            ["magick", "identify", "-format", "%wx%h", input_image]
        ).decode("utf-8").strip().split("x")
        aspect_ratio = int(width) / int(height)

        temp_diff = Path(tempdir, diff_image.name)
        compare_proc = subprocess.run(
            [
                "compare",
                "-metric",
                "AE",
                input_image,
                output_image,
                temp_diff,
            ],
            capture_output=True,
            check=False,
        )

        subprocess.check_call([
            "magick",
            input_image,
            temp_diff,
            output_image,
            "+append" if aspect_ratio < 1.5 else "-append",
            diff_image,
        ])

    pixel_difference = int(
        float(
            re.search(r"^([0-9.]+)", compare_proc.stderr.decode("utf-8").strip()).group(
                1
            )
        )
    )

    task.values["pixel_difference"] = pixel_difference

    return {"cached_answer": pixel_difference == 0}


def parse_script(scriptfile: Path) -> tuple[str, list[Path], list[Path]]:
    """Parse a script to determine (datafile, saved_vars, loaded_vars)."""
    datafile = "/Users/nwiltsie/src/project-CancerBiology-OutlierAnalysis/Figure/untracked_data/data/2024-10-08_Figure1_2_3_4_min_input.rda"

    text = scriptfile.read_text(encoding="utf-8")

    load_regex = re.compile(r"load\.multiple\.computed\.variables\(c\(([^)]+)\)\)")
    if load_match := load_regex.search(text):
        load_vars = json.loads(
            ("[" + load_match.group(1) + "]").replace("\n", " ").replace("'", '"')
        )
    else:
        load_vars = []

    cache_regex = re.compile(r"cache\.multiple\.computed\.variables\(c\(([^)]+)\)\)")
    if cache_match := cache_regex.search(text):
        cache_vars = json.loads(
            ("[" + cache_match.group(1) + "]").replace("\n", " ").replace("'", '"')
        )
    else:
        cache_vars = []

    cachedir = Path("variable-cache")

    load_files = [cachedir / (varname + ".rda") for varname in load_vars]
    cache_files = [cachedir / (varname + ".rda") for varname in cache_vars]

    return (datafile, cache_files, load_files)


def make_output_directory(targets):
    for target in targets:
        Path(target).mkdir(exist_ok=True)


def task_json_directory():
    """Create the output directory."""
    output_directory = Path(__file__).parent / "output" / "json"

    return {
        "actions": [make_output_directory],
        "clean": True,
        "targets": [output_directory],
        "task_dep": ["output_directory"],
    }

def task_diff_directory():
    """Create the output directory."""
    output_directory = Path(__file__).parent / "output" / "current_diffs"

    return {
        "actions": [make_output_directory],
        "clean": True,
        "targets": [output_directory],
        "task_dep": ["output_directory"],
    }


def task_output_directory():
    """Create the output directory."""
    output_directory = Path(__file__).parent / "output"

    return {
        "actions": [make_output_directory],
        "clean": True,
        "targets": [output_directory],
    }


def task_run_scripts():
    """Run each of the R scripts."""
    output_directory = Path(__file__).parent / "output"

    all_datafiles = set()
    all_cached_vars = set()

    all_scriptfiles = list(Path(__file__).parent.glob("Figure*/Figure*.R"))
    all_scriptfiles.extend(Path(__file__).parent.glob("01.Analysis/*.R"))

    for scriptfile in sorted(all_scriptfiles):
        plot_files = sorted(
            Path(__file__).parent.glob(f"outfigures-reference/{scriptfile.stem}*.png")
        )

        pairs = [(item, output_directory / item.name) for item in plot_files]

        try:
            datafile, targets, deps = parse_script(scriptfile)
        except json.decoder.JSONDecodeError:
            print(scriptfile)
            raise

        if datafile:
            datafile = Path(__file__).parent / "untracked_data" / "data" / datafile
            all_datafiles.add(datafile)

        all_datafiles.update(targets)
        all_datafiles.update(deps)

        all_cached_vars.update(targets)

        # Yield the task to run the script
        all_targets = [pair[1] for pair in pairs]
        all_targets.extend(pair[1].with_suffix(".pdf") for pair in pairs)
        all_targets.append(output_directory / (scriptfile.stem + ".txt"))
        all_targets.extend(targets)

        yield {
            "name": scriptfile.stem,
            "actions": [run_r],
            "meta": {"scriptfile": scriptfile, "pairs": pairs},
            "file_dep": [
                scriptfile,
                "common_functions.R",
                datafile,
                *deps,
            ],
            "task_dep": ["output_directory"],
            "targets": all_targets,
            "clean": True,
        }

    # Remove all of the undeclared variable files
    existing_cached_vars = set(Path("variable-cache").glob("*.rda"))
    for old_var in existing_cached_vars - all_cached_vars:
        print(colors.red(f"Removing outdated variable {old_var.stem}"))
        old_var.unlink()

    json_directory = Path(__file__).parent / "output" / "json"

    all_jsons = []

    for datafile in sorted(Path(datafile) for datafile in all_datafiles):
        data_json = json_directory / (datafile.stem + ".json")
        all_jsons.append(data_json)

        yield {
            "basename": "digest",
            "name": "compute-" + datafile.stem,
            "meta": {"rda": datafile, "json": data_json},
            "file_dep": [datafile],
            "actions": [get_digests],
            "targets": [data_json],
            "task_dep": ["json_directory"],
            "clean": True,
        }

    yield {
        "basename": "digest",
        "name": "consolidate",
        "file_dep": all_jsons,
        "actions": [compare_digests],
        "task_dep": ["output_directory"],
        "uptodate": [False],
        "clean": True,
    }


def compare_digests(task):
    """Compare all of the digests for multi-defined variables."""
    variable_values = {}
    for dep in task.file_dep:
        with Path(dep).open(encoding="utf-8") as infile:
            dep_variables = json.load(infile)

        for varname, value in dep_variables.items():
            value = value[0]
            variable_values.setdefault(varname, {})
            variable_values[varname].setdefault(value, [])
            variable_values[varname][value].append(Path(dep).stem)

    broken = False

    for varname, values in variable_values.items():
        if len(values) > 1:
            task.values[varname] = True
            broken = True
            print(varname)
            for index, (_, stems) in enumerate(values.items()):
                for stem in stems:
                    print(f"  {index}: {stem}")

    return not broken


def get_digests(task):
    """Get all of the digests."""
    datafile = task.meta["rda"]
    jsonfile = task.meta["json"]

    rscript = f"""\
        load(file.path("{datafile}"), new_environment <- new.env());
        vars <- ls(envir = new_environment);
        hash.list <- setNames(
            lapply(vars, function(var) {{
                digest::digest(get(var, envir = new_environment))
                }}),
            vars
            );
        jsonlite::write_json(hash.list, "{jsonfile}");
    """

    proc = subprocess.run(["Rscript", "-e", rscript], capture_output=True, check=True)

    task.values["returncode"] = proc.returncode
    task.values["stdout"] = proc.stdout.decode("utf-8")
    task.values["stderr"] = proc.stderr.decode("utf-8")

    return proc.returncode == 0


def task_compare_outputs():
    """Compare the outputs of each of the R scripts"""
    output_directory = Path(__file__).parent / "output"
    diff_directory = Path(__file__).parent / "output" / "current_diffs"

    for script_task in task_run_scripts():
        all_diffs = []

        pairs = script_task.get("meta", {}).get("pairs", [])
        if not pairs:
            continue

        for input_image, output_image in pairs:
            diff_image = (
                output_directory / f"{output_image.stem}_diff{output_image.suffix}"
            )
            all_diffs.append(diff_image)

            yield {
                "basename": script_task["name"],
                "name": output_image.stem,
                "meta": {
                    "input": input_image,
                    "output": output_image,
                    "diff": diff_image,
                },
                "actions": [compare_action],
                "file_dep": [input_image, output_image],
                "targets": [diff_image],
                "task_dep": ["output_directory", "diff_directory"],
                "clean": True,
            }

        yield {
            "basename": "summary",
            "name": script_task["name"],
            "meta": {
                "diffdir": diff_directory,
                "cleanglob": f"{script_task['name']}*.png",
            },
            "actions": [final_status],
            "getargs": {"values": (script_task["name"], "cached_answer")},
            "file_dep": all_diffs,
            "task_dep": ["output_directory", "diff_directory"],
        }


def task_style_files():
    """Style all of the R files."""
    for scriptfile in sorted(Path(__file__).parent.glob("Figure*/Figure*.R")):
        yield {
            "name": scriptfile.stem,
            "actions": [["style-file", scriptfile]],
            "file_dep": [scriptfile],
        }
    for scriptfile in sorted(Path(__file__).parent.glob("01.Analysis/*.R")):
        yield {
            "name": scriptfile.stem,
            "actions": [["style-file", scriptfile]],
            "file_dep": [scriptfile],
        }


def final_status(task, dependencies, values):
    goodfigs = []
    badfigs = []

    diff_directory = Path(task.meta["diffdir"])

    dep_map = {
        re.sub(r"_diff$", "", Path(item).stem): Path(item) for item in dependencies
    }

    for figure, result in values.items():
        copied_file = diff_directory / dep_map[figure].name

        if result:
            goodfigs.append(figure)
            copied_file.unlink(missing_ok=True)
        else:
            badfigs.append(figure)
            shutil.copyfile(dep_map[figure], copied_file)

    if badfigs:
        return FailedFigureException(goodfigs, badfigs)

    return True


if __name__ == "__main__":
    doit.run(globals())
