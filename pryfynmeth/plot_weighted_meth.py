from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_VALUE_COLUMN = "raw_weighted_methylation"


def safe_filename(text: str) -> str:
    """Convert a metadata column name into a filesystem-safe name."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    return cleaned.strip("._") or "factor"


def read_metadata(metadata_file: Path) -> pd.DataFrame:
    """Read tab-, comma-, or whitespace-delimited metadata."""
    try:
        metadata = pd.read_csv(metadata_file, sep=None, engine="python", dtype=str)
    except Exception as exc:
        raise ValueError(f"Could not read metadata file: {metadata_file}") from exc

    metadata.columns = [str(column).strip() for column in metadata.columns]
    if "sample" not in metadata.columns:
        raise ValueError("The metadata must contain a column named 'sample'.")
    if len(metadata.columns) < 2:
        raise ValueError("The metadata must contain at least one factor column.")

    metadata = metadata.apply(lambda column: column.str.strip())
    if metadata["sample"].isna().any() or (metadata["sample"] == "").any():
        raise ValueError("The metadata 'sample' column contains a missing value.")
    if metadata["sample"].duplicated().any():
        duplicates = metadata.loc[metadata["sample"].duplicated(), "sample"].tolist()
        raise ValueError(f"Duplicate samples in metadata: {', '.join(duplicates)}")
    return metadata


def resolve_sample_file(input_folder: Path, sample: str) -> Path:
    """Match a metadata sample to a file, allowing a filename or sample ID."""
    direct = input_folder / sample
    if direct.is_file():
        return direct

    candidates = [
        input_folder / f"{sample}.tsv",
        input_folder / f"{sample}_gene_weighted_methylation.tsv",
    ]
    matches = [path for path in candidates if path.is_file()]
    if len(matches) == 1:
        return matches[0]

    glob_matches = sorted(input_folder.glob(f"{sample}*"))
    glob_matches = [path for path in glob_matches if path.is_file()]
    if len(glob_matches) == 1:
        return glob_matches[0]
    if len(glob_matches) > 1:
        names = ", ".join(path.name for path in glob_matches)
        raise ValueError(f"Sample '{sample}' matches multiple files: {names}")
    raise FileNotFoundError(f"No methylation file found for sample '{sample}'.")


def calculate_sample_means(
    input_folder: Path,
    metadata: pd.DataFrame,
    value_column: str = DEFAULT_VALUE_COLUMN,
) -> pd.DataFrame:
    """Calculate one mean methylation value per biological sample."""
    rows = []
    for _, metadata_row in metadata.iterrows():
        sample = metadata_row["sample"]
        input_file = resolve_sample_file(input_folder, sample)
        table = pd.read_csv(input_file, sep="\t")

        if value_column not in table.columns:
            raise ValueError(
                f"Column '{value_column}' is missing from {input_file.name}."
            )
        values = pd.to_numeric(table[value_column], errors="coerce")
        finite_values = values[np.isfinite(values)]
        if finite_values.empty:
            raise ValueError(
                f"No finite numeric values in '{value_column}' for {input_file.name}."
            )

        row = metadata_row.to_dict()
        row["input_file"] = input_file.name
        row["n_genes_used"] = int(finite_values.size)
        row[f"mean_{value_column}"] = float(finite_values.mean())
        rows.append(row)

    return pd.DataFrame(rows)


def plot_factor(
    sample_means: pd.DataFrame,
    factor: str,
    value_column: str,
    output_file: Path,
    dpi: int = 300,
) -> None:
    """Plot group means +/- SE and overlay biological-sample values."""
    y_column = f"mean_{value_column}"
    plot_data = sample_means.dropna(subset=[factor, y_column]).copy()
    if plot_data.empty:
        raise ValueError(f"No usable values are available for factor '{factor}'.")

    # Preserve the order in which factor levels first appear in the metadata.
    levels = list(pd.unique(plot_data[factor]))
    grouped = plot_data.groupby(factor, sort=False, observed=True)[y_column]
    means = grouped.mean().reindex(levels)
    counts = grouped.count().reindex(levels)
    ses = grouped.std(ddof=1).reindex(levels) / np.sqrt(counts)
    ses = ses.fillna(0.0)

    width = max(6.5, 1.15 * len(levels))
    fig, ax = plt.subplots(figsize=(width, 5.5))
    x = np.arange(len(levels), dtype=float)
    ax.bar(
        x,
        means.to_numpy(),
        yerr=ses.to_numpy(),
        capsize=5,
        width=0.68,
        color="#4C78A8",
        edgecolor="black",
        linewidth=0.8,
        zorder=2,
    )

    rng = np.random.default_rng(12345)
    for position, level in enumerate(levels):
        points = plot_data.loc[plot_data[factor] == level, y_column].to_numpy()
        jitter = rng.uniform(-0.11, 0.11, size=len(points))
        ax.scatter(
            position + jitter,
            points,
            s=30,
            facecolor="white",
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
        )

    ax.set_xticks(x, [str(level) for level in levels])
    ax.set_xlabel(factor.replace("_", " ").title())
    ax.set_ylabel(f"Mean {value_column.replace('_', ' ')}")
    ax.set_title(f"{value_column.replace('_', ' ').title()} by {factor}")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7, zorder=0)
    fig.tight_layout()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_gene_weighted_methylation(
    input_folder: str | Path,
    metadata_file: str | Path,
    output_folder: str | Path,
    value_column: str = DEFAULT_VALUE_COLUMN,
    dpi: int = 300,
) -> pd.DataFrame:
    """Create one PNG per metadata factor and return the sample-level table."""
    input_folder = Path(input_folder)
    metadata_file = Path(metadata_file)
    output_folder = Path(output_folder)

    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input folder does not exist: {input_folder}")
    if not metadata_file.is_file():
        raise FileNotFoundError(f"Metadata file does not exist: {metadata_file}")
    output_folder.mkdir(parents=True, exist_ok=True)

    metadata = read_metadata(metadata_file)
    sample_means = calculate_sample_means(input_folder, metadata, value_column)
    audit_file = output_folder / "sample_mean_methylation.tsv"
    sample_means.to_csv(audit_file, sep="\t", index=False)

    for factor in (column for column in metadata.columns if column != "sample"):
        output_file = output_folder / (
            f"{safe_filename(factor)}_{safe_filename(value_column)}.png"
        )
        plot_factor(sample_means, factor, value_column, output_file, dpi=dpi)

    return sample_means


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate each sample's mean gene-weighted methylation and create "
            "one mean +/- SE bar graph for every metadata factor."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_folder",
        required=True,
        type=Path,
        help="Folder containing gene-weighted methylation TSV files.",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        dest="metadata_file",
        required=True,
        type=Path,
        help="Metadata text/TSV/CSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_folder",
        required=True,
        type=Path,
        help="Folder for PNGs and the sample-level audit TSV.",
    )
    parser.add_argument(
        "--value-column",
        default=DEFAULT_VALUE_COLUMN,
        help=f"Methylation column to summarize (default: {DEFAULT_VALUE_COLUMN}).",
    )
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution (default: 300).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = plot_gene_weighted_methylation(
        args.input_folder,
        args.metadata_file,
        args.output_folder,
        value_column=args.value_column,
        dpi=args.dpi,
    )
    factors = len(result.columns) - 4
    print(f"Processed {len(result)} samples and plotted {factors} metadata factors.")
    print(f"Outputs written to: {args.output_folder.resolve()}")


if __name__ == "__main__":
    main()
