#!/usr/bin/env python3
"""Compare output plots in two folders."""

import argparse
import re
import subprocess
import tempfile

from pathlib import Path


def compare(original: Path, updated: Path, diff_folder: Path):
    """Compare the plots in two folders."""
    original_pngs = {
        item.name.split("_", maxsplit=1)[-1]: item
        for item in list(original.glob("*.png"))
    }
    updated_pngs = {
        item.name.split("_", maxsplit=1)[-1]: item
        for item in list(updated.glob("*.png"))
    }

    with tempfile.TemporaryDirectory() as tempdir:
        for key in sorted(original_pngs.keys() | updated_pngs.keys()):
            if key not in original_pngs:
                print("ONLY IN UPDATED:")
                print("\t", updated_pngs[key])
                continue

            if key not in updated_pngs:
                print("ONLY IN ORIGINAL:")
                print("\t", original_pngs[key])
                continue

            diff_image = Path(tempdir) / key

            compare_proc = subprocess.run(
                [
                    "compare",
                    "-metric",
                    "AE",
                    original_pngs[key],
                    updated_pngs[key],
                    diff_image,
                ],
                capture_output=True,
                check=False,
            )

            pixel_difference = int(
                float(
                    re.search(
                        r"^([0-9.]+)", compare_proc.stderr.decode("utf-8").strip()
                    ).group(1)
                )
            )

            if pixel_difference:
                print("DIFFERENT:")
                print("\t", original_pngs[key].name)
                print("\t", updated_pngs[key].name)
                diff_folder.mkdir(exist_ok=True)
                diff_image.rename(diff_folder / key)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("original", type=Path)
    parser.add_argument("updated", type=Path)
    parser.add_argument("--diff-dir", type=Path, default=Path("diffs"))

    args = parser.parse_args()

    compare(args.original, args.updated, args.diff_dir)