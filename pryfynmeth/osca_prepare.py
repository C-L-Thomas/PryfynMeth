from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path


REQUIRED_COLUMNS = {"chr", "pos", "C", "total"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine all .txt methylation-count files in a folder into a "
            "tab-separated sample-by-site OSCA input matrix."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Folder containing per-sample .txt files."
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path,
        help="Output table path."
    )
    parser.add_argument(
        "--missing", default="NA",
        help="Value used when a site is absent or total is zero (default: NA)."
    )
    return parser.parse_args()


def natural_key(value: str) -> tuple:
    """Sort chromosome labels naturally (e.g. chr2 before chr10)."""
    return tuple(int(part) if part.isdigit() else part.lower()
                 for part in re.split(r"(\d+)", value))


def read_sample(path: Path) -> dict[tuple[str, int], float | None]:
    values: dict[tuple[str, int], float | None] = {}

    with path.open("r", encoding="utf-8-sig") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError("file is empty")

        header = header_line.split()
        missing_columns = REQUIRED_COLUMNS.difference(header)
        if missing_columns:
            raise ValueError(
                "missing required column(s): " + ", ".join(sorted(missing_columns))
            )
        column = {name: header.index(name) for name in REQUIRED_COLUMNS}

        for line_number, line in enumerate(handle, start=2):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            try:
                chrom = fields[column["chr"]]
                pos = int(fields[column["pos"]])
                c_count = float(fields[column["C"]])
                total = float(fields[column["total"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"invalid data on line {line_number}: {line.strip()}") from exc

            site = (chrom, pos)
            if site in values:
                raise ValueError(
                    f"duplicate site {chrom}:{pos} on line {line_number}"
                )
            if c_count < 0 or total < 0 or c_count > total:
                raise ValueError(
                    f"invalid counts at {chrom}:{pos}: C={c_count:g}, total={total:g}"
                )
            values[site] = None if total == 0 else c_count / total

    return values


def format_value(value: float | None, missing: str) -> str:
    if value is None or not math.isfinite(value):
        return missing
    return f"{value:.2f}"


def main() -> int:
    args = parse_args()
    if not args.input.is_dir():
        print(f"Error: input folder does not exist: {args.input}", file=sys.stderr)
        return 2

    output_resolved = args.output.resolve()
    input_files = sorted(
        (path for path in args.input.glob("*.txt")
         if path.is_file() and path.resolve() != output_resolved),
        key=lambda path: natural_key(path.name),
    )
    if not input_files:
        print(f"Error: no .txt files found in {args.input}", file=sys.stderr)
        return 2

    samples: list[tuple[str, dict[tuple[str, int], float | None]]] = []
    all_sites: set[tuple[str, int]] = set()
    sample_ids: set[str] = set()

    for path in input_files:
        sample_id = path.stem
        if sample_id in sample_ids:
            print(f"Error: duplicate IID derived from filename: {sample_id}", file=sys.stderr)
            return 2
        try:
            values = read_sample(path)
        except (OSError, ValueError) as exc:
            print(f"Error reading {path}: {exc}", file=sys.stderr)
            return 2
        sample_ids.add(sample_id)
        samples.append((sample_id, values))
        all_sites.update(values)

    sites = sorted(all_sites, key=lambda site: (natural_key(site[0]), site[1]))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    try:
        with args.output.open("w", encoding="utf-8", newline="") as output:
            output.write("IID\t" + "\t".join(f"{chrom}:{pos}" for chrom, pos in sites) + "\n")
            for sample_id, values in samples:
                row = [sample_id]
                row.extend(
                    format_value(values.get(site), args.missing)
                    for site in sites
                )
                output.write("\t".join(row) + "\n")
    except OSError as exc:
        print(f"Error writing {args.output}: {exc}", file=sys.stderr)
        return 2

    print(
        f"Wrote {len(samples)} samples and {len(sites)} sites to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
