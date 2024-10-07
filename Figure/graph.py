#!/usr/bin/env python3

import argparse
import datetime
import itertools
import json
import re
import subprocess

from collections import defaultdict, Counter
from pathlib import Path

import graphviz
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def get_variables(names: list[str]) -> list[str]:
    """Extract the variable names from this list."""
    varnames = []
    for name in names:
        path = Path(name)
        if path.suffix == ".rda" and "2024-" not in path.stem:
            varnames.append(path.stem)

    return varnames


def keep_name(name: str) -> bool:
    """Return True if the name should be included."""
    # return not ("run_script" not in name and "patches" not in name and "summary:" not in name)
    return not ("run_scripts:" not in name and "patches" not in name)


def make_plots(datafile: Path):
    """Make all of the plots."""
    with datafile.open(encoding="utf-8") as infile:
        data = json.load(infile)

    plot_dependency_graph(data)
    filenames = plot_duration_graph(data)
    filenames.extend(plot_line_speed(data))
    filenames.extend(plot_file_lengths(data))
    filenames.extend(plot_ref_counts(data))

    if filenames:
        subprocess.check_call(["open", *filenames])


def plot_ref_counts(data) -> list[str]:
    """Plot the number of times a variable is referenced."""
    use_count = Counter()
    common_vars = set()

    for name, task in data["tasks"].items():
        if not name.startswith("run_scripts:"):
            continue
        use_count.update(get_variables(task["file_dep"]))

        if (targets := get_variables(task["targets"])) and "Figure" not in name:
            common_vars.update(targets)

    labels = []
    counts = []
    colors = []

    for varname, count in use_count.most_common():
        labels.append(varname)
        counts.append(count)
        colors.append("blue" if varname in common_vars else "red")

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.set_dpi(200)
    ax.barh(
        labels,
        counts,
        color=colors
    )

    ax.invert_yaxis()
    ax.set_xlabel("Use count")
    ax.grid()
    ax.set_axisbelow(True)

    plt.savefig("variables.png", bbox_inches="tight")
    plt.close(fig)
    return ["variables.png"]



def plot_file_lengths(data) -> list[str]:
    """Plot the file lengths."""
    datapoints = []

    for name, task in data["tasks"].items():
        name_match = re.match(r"^run_scripts:Figure(\d\w+)$", name)

        if not name_match:
            continue

        fignum = name_match.group(1)
        scriptfile = Path(task["meta"]["scriptfile"])
        with scriptfile.open(mode="r", encoding="utf-8") as infile:
            datapoints.append((fignum, len(infile.readlines())))

    datapoints.sort()

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.set_dpi(200)
    ax.barh(
        [point[0] for point in datapoints],
        [point[1] for point in datapoints],
    )

    ax.invert_yaxis()
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel("Lines of code")

    plt.savefig("lengths.png", bbox_inches="tight")
    plt.close(fig)
    return ["lengths.png"]
    # subprocess.check_call(["open", "durations.png"])


def plot_duration_graph(data) -> list[str]:
    """Make the duration barplot."""
    datapoints = []

    fig_colors = {
        "1": "#a6cee3",
        "2": "#1f78b4",
        "3": "#b2df8a",
        "4": "#33a02c",
    }

    for name, task in data["tasks"].items():
        if not name.startswith("run_scripts:"):
            continue

        if name_match := re.match(r"^run_scripts:Figure(\d)(\w+)$", name):
            fignum = name_match.group(1)
            subfig = name_match.group(2)
            color = fig_colors[fignum]
        else:
            color = "blue"

        try:
            duration = task["values"]["duration"]
        except KeyError:
            duration = 60
            color = "red"

        datapoints.append((
            re.sub("^run_scripts:(?:Figure)?", "", name),
            duration,
            color,
        ))

    datapoints.sort()

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.set_dpi(200)
    ax.barh(
        [point[0] for point in datapoints],
        [point[1] for point in datapoints],
        # label=[point[1] for point in datapoints],
        color=[point[2] for point in datapoints]
    )

    ax.invert_yaxis()
    ax.set_xlabel("Runtime (s)")
    ax.grid()
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(60))

    plt.savefig("durations.png", bbox_inches="tight")
    plt.close(fig)
    return ["durations.png"]
    # subprocess.check_call(["open", "durations.png"])


def plot_line_speed(data: dict) -> list[str]:
    """Make a plot of line speed."""
    regex = re.compile(r"Time: (.{19}\.\d+) \| Line: (\d+)")

    all_points = []

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for name, task in data["tasks"].items():
        if not name.startswith("run_scripts"):
            continue

        duration = task["values"].get("duration", 0)
        if not duration:
            continue

        stderr = task["values"].get("stderr", "")

        points = []

        for match in regex.findall(stderr):
            points.append((
                datetime.datetime.fromisoformat(match[0]),
                int(match[1]) + 1
            ))

        if points:
            all_points.append((duration, re.sub("^run_scripts:(?:Figure)?", "", name), points))
            continue

    all_points.sort(key=lambda x: x[0])
    if len(all_points) > 3:
        all_points = all_points[-3:]

    if not all_points:
        return []

    fig, axs = plt.subplots(nrows=len(all_points), ncols=1, figsize=(20, 10), squeeze=False)
    fig.set_dpi(200)
    for ax_array, (_, label, points) in zip(axs, all_points):
        print(f"Plotting {label}")
        ax = ax_array[0]
        patches = []

        for index, (alpha, beta) in enumerate(itertools.pairwise(points)):
            delta_time = (beta[0] - alpha[0]).total_seconds()
            delta_lines = beta[1] - alpha[1]
            time_per_line = delta_time / delta_lines

            patches.append(Rectangle(
                (alpha[1], 0),
                delta_lines,
                time_per_line,
                color=color_cycle[index % len(color_cycle)],
                linewidth=0,
            ))

        # fig, ax = plt.subplots(figsize=(20, 10))
        pc = PatchCollection(patches, match_original=True)
        ax.add_collection(pc)
        ax.autoscale_view()
        ax.grid()
        ax.set_axisbelow(True)
        ax.set_ylim(0, )
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.set(title = label)

    filename = "linespeed.png"
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)
    return [filename]


def plot_dependency_graph(data: dict):
    """Plot the thing."""
    dot = graphviz.Digraph(engine="dot", format="png")
    dot.attr(rankdir="LR")

    target_vars = {
        "depmap.drug.info.match.sanger.dup",
        "cas.effect.breast.05.na",
        "cas.effect.breast.05.na",
        "gene.dependency.breast.t.num.match.05.na",
        "gene.dependency.diff.matrix.05.overlap",
        "gene.effect.diff.matrix.05.overlap",
        "gene.rnai.breast.t.num.match.05.na",
        "gene.rnai.diff.matrix.05.overlap",
        "rnai.05.box",
        "rnai.score.05.overlap.minus.05",
    }

    cached_vars = {}
    seen_vars = set()

    image_status = defaultdict(list)

    bad_vars = data["tasks"]["digest:consolidate"]["values"].keys()

    # Find the image statuses
    for name, task in data["tasks"].items():
        if not re.match(r"^Figure.*:Figure.*$", name):
            continue

        figure = name.split(":")[0]
        try:
            image_status[figure].append(task["values"]["pixel_difference"])
        except KeyError:
            image_status[figure].append(-1)

    for name, task in data["tasks"].items():
        if not keep_name(name):
            continue

        name = name.replace(":", "/")

        for varname in get_variables(task["targets"]):
            if varname in bad_vars:
                cached_vars[varname] = (name, "red")
            else:
                cached_vars[varname] = (name, "white")

            seen_vars.add(varname)
        
        try:
            if task["values"]["returncode"] == 0:
                if task["values"]["duration"] > 180:
                    color = "orange"
                else:
                    color = "#99d594"
            else:
                color = "#fc8d59"
                print(name)
                print(task["values"]["stderr"])

        except KeyError:
            color = "gray"

        kwargs = {
            "label": name.rsplit("/", maxsplit=1)[-1],
            "shape": "rect",
            "style": "rounded,filled",
            "fillcolor": color,
        }

        if name.startswith("run_scripts/") and color != "gray":
            figure = name.split("/")[-1]
            statuses = image_status.get(figure)
            if statuses:
                if max(statuses) > 0:
                    kwargs["color"] = "red"
                    kwargs["penwidth"] = "10"
                    # color += ":#fc8d59"
                elif min(statuses) < 0:
                    kwargs["color"] = "black"
                    kwargs["penwidth"] = "10"

        if "Figure" not in name:
            kwargs["shape"] = "polygon"
            kwargs["fillcolor"] = "lightgreen"

        dot.node(name, **kwargs)

        def varnamesortkey(item):
            return list(reversed(item[0].split(".")))

        for varname, (source, color) in sorted(list(cached_vars.items()), key=varnamesortkey):
            if source == name:
                cached_vars.pop(varname)
                dot.node(varname, fillcolor=color, style="filled", shape="cds")
                dot.edge(name, varname)

        for varname in get_variables(task["file_dep"]):
            if varname in cached_vars:
                source, color = cached_vars.pop(varname)
                dot.node(varname, fillcolor=color, style="filled", shape="cds")
                dot.edge(source, varname)
            elif varname not in seen_vars:
                dot.node(varname, shape="cds", style="filled", fillcolor="#ffffbf")
                seen_vars.add(varname)

            dot.edge(varname, name)

    for varname, (source, color) in cached_vars.items():
        dot.node(varname, fillcolor="red", style="filled", shape="cds")
        dot.edge(source, varname)

    dot.render(view=True)
    # dot.render()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--datafile", type=Path)
    args = parser.parse_args()

    if not args.datafile:
        args.datafile = sorted(Path(__file__).parent.glob("2024*.json"))[-1]

    make_plots(args.datafile)
