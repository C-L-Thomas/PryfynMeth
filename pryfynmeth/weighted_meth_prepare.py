from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class Gene:
    chromosome: str
    start: int
    end: int
    strand: str
    gene_id: str


@dataclass
class GeneCounts:
    raw_c: int = 0
    corrected_c: int = 0
    total: int = 0
    n_observed_sites: int = 0
    n_fdr_significant_sites: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Assign strand-specific methylation sites to genes and calculate "
            "paper-style FDR-corrected and raw weighted methylation per sample."
        )
    )
    parser.add_argument("-i", "--input-folder", required=True, type=Path)
    parser.add_argument("-g", "--genes", required=True, type=Path,
                        help="Gene-position TSV: chromosome,start,end,strand,gene_id")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Output directory; one TSV is written per sample")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="FDR threshold; comparison is strictly below (default: 0.05)")
    parser.add_argument("--min-coverage", type=int, default=1,
                        help="Minimum total coverage for a site to contribute (default: 1)")
    parser.add_argument("--max-coverage", type=int,
                        help="Optional maximum total coverage per site")
    parser.add_argument("--min-sites", type=int, default=1,
                        help="Minimum observed strand-specific sites per sample-gene (default: 1)")
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def load_genes(path: Path) -> dict[str, list[Gene]]:
    genes: dict[str, list[Gene]] = defaultdict(list)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "start", "end", "strand", "gene_id"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Gene file is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            gene = Gene(
                chromosome=row["chromosome"],
                start=int(row["start"]),
                end=int(row["end"]),
                strand=row["strand"],
                gene_id=row["gene_id"],
            )
            genes[gene.chromosome].append(gene)
    for chromosome in genes:
        genes[chromosome].sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))
    return dict(genes)


def sample_name(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in (".txt", ".tsv"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    if name.endswith("_binomial_results"):
        name = name[:-len("_binomial_results")]
    return name


def has_methylation_header(path: Path) -> bool:
    """Return True when a file has the required methylation-table columns."""
    try:
        with open_text(path) as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile):
        return False
    required = {"chr", "pos", "strand", "C", "T", "total", "fdr"}
    return required.issubset(header)


def discover_input_files(folder: Path, genes: Path) -> list[Path]:
    """Find methylation tables while ignoring unrelated files in the folder."""
    excluded = {genes.resolve()}
    candidates: list[Path] = []
    for path in sorted(folder.iterdir()):
        if not path.is_file() or path.resolve() in excluded:
            continue
        lower_name = path.name.lower()
        if not lower_name.endswith((".txt", ".tsv", ".txt.gz", ".tsv.gz")):
            continue
        if has_methylation_header(path):
            candidates.append(path)
    return candidates


def process_sample(
    path: Path,
    genes_by_chr: dict[str, list[Gene]],
    fdr_threshold: float,
    min_coverage: int,
    max_coverage: int | None,
) -> dict[Gene, GeneCounts]:
    results: dict[Gene, GeneCounts] = defaultdict(GeneCounts)
    state_chr: str | None = None
    gene_index = 0
    active: list[Gene] = []
    previous_position = -1
    completed_chromosomes: set[str] = set()

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chr", "pos", "strand", "C", "T", "total", "fdr"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name} is missing columns: {', '.join(sorted(missing))}")

        for line_number, row in enumerate(reader, start=2):
            chromosome = row["chr"]
            position = int(row["pos"])

            if chromosome != state_chr:
                if state_chr is not None:
                    completed_chromosomes.add(state_chr)
                if chromosome in completed_chromosomes:
                    raise ValueError(
                        f"{path.name}:{line_number}: chromosome {chromosome} is not contiguous; "
                        "sort input by chromosome and position"
                    )
                state_chr = chromosome
                gene_index = 0
                active = []
                previous_position = -1
            elif position < previous_position:
                raise ValueError(
                    f"{path.name}:{line_number}: positions are not sorted within {chromosome}"
                )
            previous_position = position

            total = int(row["total"])
            if total < min_coverage or (max_coverage is not None and total > max_coverage):
                continue

            raw_c = int(row["C"])
            if raw_c < 0 or total < raw_c:
                raise ValueError(f"{path.name}:{line_number}: invalid C/total counts")

            try:
                fdr = float(row["fdr"])
            except (TypeError, ValueError):
                fdr = 1.0
            significant = fdr < fdr_threshold
            corrected_c = raw_c if significant else 0

            chromosome_genes = genes_by_chr.get(chromosome, [])
            while gene_index < len(chromosome_genes) and chromosome_genes[gene_index].start <= position:
                active.append(chromosome_genes[gene_index])
                gene_index += 1
            if active:
                active = [gene for gene in active if gene.end >= position]

            for gene in active:
                if row["strand"] != gene.strand:
                    continue
                counts = results[gene]
                counts.raw_c += raw_c
                counts.corrected_c += corrected_c
                counts.total += total
                counts.n_observed_sites += 1
                counts.n_fdr_significant_sites += int(significant)

    return results


def main() -> int:
    args = parse_args()
    if not 0 <= args.fdr <= 1:
        raise ValueError("--fdr must be between 0 and 1")
    if args.min_coverage < 1:
        raise ValueError("--min-coverage must be at least 1")
    if args.min_sites < 1:
        raise ValueError("--min-sites must be at least 1")

    genes_by_chr = load_genes(args.genes)
    files = discover_input_files(args.input_folder, args.genes)
    if not files:
        raise FileNotFoundError(
            f"No methylation tables with the required columns were found in "
            f"{args.input_folder}"
        )

    args.output.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample", "chromosome", "gene_start", "gene_end", "gene_strand",
        "gene_id", "raw_C", "corrected_C", "total",
        "n_observed_sites", "n_fdr_significant_sites", "raw_weighted_methylation",
        "corrected_weighted_methylation",
    ]

    seen_samples: set[str] = set()
    for index, path in enumerate(files, start=1):
        sample = sample_name(path)
        if sample in seen_samples:
            raise ValueError(f"More than one input file produced sample name {sample!r}")
        seen_samples.add(sample)
        print(f"[{index}/{len(files)}] Processing {sample}", file=sys.stderr)
        results = process_sample(
            path, genes_by_chr, args.fdr, args.min_coverage, args.max_coverage
        )
        output_path = args.output / f"{sample}_gene_weighted_methylation.tsv"
        with output_path.open("w", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for gene, counts in sorted(
                results.items(), key=lambda item: (
                    item[0].chromosome, item[0].start, item[0].end, item[0].gene_id
                )
            ):
                if counts.n_observed_sites < args.min_sites:
                    continue
                writer.writerow({
                    "sample": sample,
                    "chromosome": gene.chromosome,
                    "gene_start": gene.start,
                    "gene_end": gene.end,
                    "gene_strand": gene.strand,
                    "gene_id": gene.gene_id,
                    "raw_C": counts.raw_c,
                    "corrected_C": counts.corrected_c,
                    "total": counts.total,
                    "n_observed_sites": counts.n_observed_sites,
                    "n_fdr_significant_sites": counts.n_fdr_significant_sites,
                    "raw_weighted_methylation": counts.raw_c / counts.total,
                    "corrected_weighted_methylation": counts.corrected_c / counts.total,
                })
        print(f"Wrote {output_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(2)
