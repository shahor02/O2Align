#!/usr/bin/env python3
"""Add comments from an mp2param.txt file to Millepede results."""

import argparse
from pathlib import Path


def read_comments(path: Path) -> dict[int, str]:
    """Read label-to-comment mappings from an annotated parameter file."""
    comments: dict[int, str] = {}
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("!"):
                continue

            fields = stripped.split(maxsplit=3)
            try:
                label = int(fields[0])
            except (ValueError, IndexError):
                continue

            # The comment begins with ! and may contain arbitrary whitespace.
            if len(fields) == 4 and fields[3].lstrip().startswith("!"):
                comments[label] = fields[3].lstrip()

    return comments


def result_label(line: str) -> int | None:
    fields = line.split(maxsplit=1)
    if not fields:
        return None
    try:
        return int(fields[0])
    except ValueError:
        return None


def output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_labeled{input_path.suffix}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Append parameter comments to alignment results by label."
    )
    parser.add_argument("results", nargs="?", default="millepede.res")
    parser.add_argument("parameters", nargs="?", default="mp2param.txt")
    args = parser.parse_args()

    results_path = Path(args.results)
    parameters_path = Path(args.parameters)
    comments = read_comments(parameters_path)
    destination = output_path(results_path)

    labels_in_results: set[int] = set()
    missing_comments: set[int] = set()
    with results_path.open(encoding="utf-8") as source, destination.open(
        "w", encoding="utf-8"
    ) as output:
        for line in source:
            label = result_label(line)
            if label is None:
                output.write(line)
                continue

            labels_in_results.add(label)
            comment = comments.get(label, "! N/A")
            if label not in comments:
                missing_comments.add(label)

            content = line.rstrip("\r\n")
            content = content.rstrip()
            # Column 60 is character offset 59 (columns are one-based).
            if len(content) < 59:
                content = content.ljust(59)
            else:
                content += " "
            output.write(f"{content}{comment}\n")

    missing_results = sorted(comments.keys() - labels_in_results)
    print(f"labels seen in {results_path} but missing in {parameters_path}")
    ordered_missing_comments = sorted(missing_comments)
    for start in range(0, len(ordered_missing_comments), 5):
        print("  ".join(map(str, ordered_missing_comments[start : start + 5])))

    print(f"labels seen in {parameters_path} but missing in {results_path}")
    for label in missing_results:
        print(f"{label}  {comments[label]}")


if __name__ == "__main__":
    main()
