#!/usr/bin/env python3
"""
PCA from per-sample site tables using only NumPy/Pandas (no scikit-learn)

Input files (one per sample), columns (case-insensitive):
chr, pos, strand, C, T, total, p-value, fdr
We use chr/pos/strand/C/total.

Pipeline:
1) Build unique site key "modifier" = chr:pos:strand
2) Value = (C / total) * 100 by default (use --as-fraction for 0–1)
3) Assemble samples × sites matrix (rows = samples, cols = sites)
4) Optional filtering (min_total, min_samples)
5) Handle missing (drop or impute mean/median/most_frequent)
6) Optional scaling (z-score per site)
7) PCA via SVD (NumPy)
8) (NEW) Optional metadata: color points by condition in plots and export joined scores

Example metadata (TSV):
sample	condition
Nanopore_Female_1.txt	1
Nanopore_Female_2.txt	1
...
WGBS_Male_3.txt	4
"""

import argparse
import glob
import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Optional plotting: import matplotlib only if available
HAVE_MPL = True
def _lazy_import_matplotlib():
    global HAVE_MPL
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except Exception:
        HAVE_MPL = False


def parse_args():
    p = argparse.ArgumentParser(description="Run PCA on C/total percentages across samples (no scikit-learn).")
    # data
    p.add_argument("--input-dir", "-i", required=True, help="Directory with per-sample files.")
    p.add_argument("--glob", default="*.tsv", help='Glob for files in input-dir (e.g., "*.tsv").')
    p.add_argument("--delimiter", default="\t", help="Input delimiter for sample files (default: tab).")
    p.add_argument("--as-fraction", action="store_true", help="Use raw fraction (0–1) instead of percentage (0–100).")
    p.add_argument("--min-total", type=int, default=0, help="Drop rows with total <= this (default: 0).")
    p.add_argument("--min-samples", type=int, default=1,
                   help="Keep sites present (non-NaN) in at least this many samples (default: 1).")
    # missing/scaling
    p.add_argument("--dropna", action="store_true",
                   help="Drop samples with any missing values (overrides --impute).")
    p.add_argument("--impute", choices=["mean", "median", "most_frequent"],
                   help="Impute missing values per site (column).")
    p.add_argument("--scale", action="store_true", help="Z-score each site (feature) across samples before PCA.")
    p.add_argument("--n-components", type=int, default=2, help="Number of PCs to compute (default: 2).")
    # plotting
    p.add_argument("--plot", choices=["2d", "3d", "none"], default="2d", help="Plot PCA scores.")
    p.add_argument("--annotate-samples", action="store_true", help="Overlay sample names on plots.")
    # metadata
    p.add_argument("--metadata", "-m", help="Path to metadata file (TSV/CSV) with columns: sample, condition.")
    p.add_argument("--metadata-delimiter", default="\t", help="Delimiter for metadata file (default: tab).")
    p.add_argument("--metadata-sample-col", default="sample", help="Column name in metadata for sample (default: sample).")
    p.add_argument("--metadata-condition-col", default="condition", help="Column name in metadata for condition (default: condition).")
    # output
    p.add_argument("--output-dir", "-o", default="pca_numpy_output", help="Where to save outputs.")
    return p.parse_args()


def standardize_headers(df: pd.DataFrame) -> pd.DataFrame:
    mapping = {c: c.strip().lower().replace("-", "_") for c in df.columns}
    return df.rename(columns=mapping)


def sample_name_from_path(path: str) -> str:
    # use filename stem (no extension) as sample id
    return Path(path).stem


def read_and_summarize_sample(path: str, sep: str, as_fraction: bool, min_total: int) -> pd.Series:
    df = pd.read_csv(path, sep=sep, dtype={"chr": str, "strand": str})
    df = standardize_headers(df)
    for col in ["chr", "pos", "strand", "c", "total"]:
        if col not in df.columns:
            raise ValueError(f"{path} missing required column: {col}")

    df = df[["chr", "pos", "strand", "c", "total"]].dropna()
    df["pos"] = df["pos"].astype(int)
    df["c"] = pd.to_numeric(df["c"], errors="coerce")
    df["total"] = pd.to_numeric(df["total"], errors="coerce")
    df = df.dropna(subset=["c", "total"])
    df = df[df["total"] > min_total]

    if df.empty:
        return pd.Series(dtype=float, name=sample_name_from_path(path))

    # Unique site key
    df["modifier"] = df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":" + df["strand"].astype(str)

    # If duplicate sites in a file, sum counts then compute fraction
    grp = df.groupby("modifier", as_index=False)[["c", "total"]].sum()
    frac = grp["c"] / grp["total"]
    values = frac if as_fraction else frac * 100.0
    return pd.Series(values.values, index=grp["modifier"].values, name=sample_name_from_path(path))


def build_matrix(input_dir: str, pattern: str, sep: str, as_fraction: bool, min_total: int) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        raise FileNotFoundError(f"No files matched {os.path.join(input_dir, pattern)}")
    series = [read_and_summarize_sample(f, sep, as_fraction, min_total) for f in files]
    mat = pd.concat(series, axis=1).T  # rows=samples, cols=sites
    return mat


def filter_sites_by_presence(X: pd.DataFrame, min_samples: int) -> pd.DataFrame:
    if min_samples <= 1:
        return X
    keep = (X.notna().sum(axis=0) >= min_samples)
    Xf = X.loc[:, keep]
    print(f"[filter] kept {Xf.shape[1]} / {X.shape[1]} sites with >= {min_samples} samples present")
    return Xf


def impute_matrix(X: pd.DataFrame, strategy: str) -> pd.DataFrame:
    Xf = X.copy()
    for col in Xf.columns:
        s = Xf[col]
        if strategy == "mean":
            fill = s.mean(skipna=True)
        elif strategy == "median":
            fill = s.median(skipna=True)
        else:  # most_frequent
            mode = s.mode(dropna=True)
            fill = mode.iloc[0] if len(mode) else 0.0
        Xf[col] = s.fillna(fill)
    return Xf


def prepare_for_pca(X: pd.DataFrame, dropna: bool, impute: Optional[str], scale: bool) -> pd.DataFrame:
    # Handle missing
    if dropna:
        Xp = X.dropna(axis=0).copy()  # drop samples with any NaN
        if Xp.shape[0] < 2:
            raise ValueError("Fewer than 2 samples remain after dropping rows with NaNs.")
    else:
        if X.isna().any().any():
            if impute:
                Xp = impute_matrix(X, impute)
            else:
                top = X.isna().sum().sort_values(ascending=False).head(10)
                raise ValueError("Missing values detected. Use --dropna or --impute.\n"
                                 f"Top per-site missing counts:\n{top.to_string()}")
        else:
            Xp = X.copy()

    # Centering is essential; scaling optional
    if scale:
        means = Xp.mean(axis=0)
        stds = Xp.std(axis=0, ddof=0).replace(0, 1.0)
        Xs = (Xp - means) / stds
    else:
        Xs = Xp - Xp.mean(axis=0)

    return Xs


def pca_via_svd(X: pd.DataFrame, n_components: int):
    """
    PCA using thin SVD on centered (or z-scored) matrix X (rows=samples, cols=features).

    X = U Σ V^T
    Scores (per-sample) = U Σ
    Loadings (per-feature) = V
    Explained variance ratio = Σ^2 / (n-1) / sum(Σ^2 / (n-1))
    """
    X_np = X.values.astype(float, copy=False)
    n_samples, n_features = X_np.shape
    kmax = min(n_samples, n_features)
    if kmax < 2:
        raise ValueError("Need at least 2 samples and 2 sites to run PCA.")

    # Thin SVD
    U, S, Vt = np.linalg.svd(X_np, full_matrices=False)
    k = min(n_components, kmax)

    # Scores and loadings
    scores = U[:, :k] * S[:k]  # (n_samples x k)
    loadings = Vt[:k, :].T     # (n_features x k)

    # Explained variance (match sklearn): S^2 / (n_samples - 1)
    ev = (S ** 2) / (n_samples - 1)
    ev_ratio = ev / ev.sum()
    explained = ev_ratio[:k]
    cum_explained = np.cumsum(explained)

    comp_names = [f"PC{i+1}" for i in range(k)]
    scores_df = pd.DataFrame(scores, index=X.index, columns=comp_names)
    loadings_df = pd.DataFrame(loadings, index=X.columns, columns=comp_names)
    return scores_df, loadings_df, explained, cum_explained


def ensure_out_dir(path: str) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def load_metadata_series(path: str, sep: str, sample_col: str, condition_col: str) -> pd.Series:
    """
    Returns a Series indexed by sample *stem* (filename without extension),
    with values = condition (any dtype). Extra rows/cols are ignored.
    """
    md = pd.read_csv(path, sep=sep, dtype=str)
    md = standardize_headers(md)
    if sample_col.lower() not in md.columns or condition_col.lower() not in md.columns:
        raise ValueError(f"Metadata needs columns '{sample_col}' and '{condition_col}' (case-insensitive).")
    s = md[sample_col.lower()].astype(str).apply(lambda x: Path(x).stem)
    cond = md[condition_col.lower()]
    meta = pd.Series(cond.values, index=s, name="condition")
    return meta


def save_tables(out_dir: Path, matrix_samples_x_sites: pd.DataFrame,
                scores_df: pd.DataFrame, loadings_df: pd.DataFrame,
                explained: np.ndarray, cum_explained: np.ndarray,
                meta: Optional[pd.Series] = None) -> None:
    mat_path = out_dir / "matrix_samples_x_sites.csv"
    scores_path = out_dir / "pca_scores_per_sample.csv"
    loadings_path = out_dir / "pca_loadings_per_site.csv"
    var_path = out_dir / "pca_explained_variance.csv"

    matrix_samples_x_sites.to_csv(mat_path, index=True)
    scores_df.to_csv(scores_path, index=True)
    loadings_df.to_csv(loadings_path, index=True)
    pd.DataFrame({
        "component": [f"PC{i+1}" for i in range(len(explained))],
        "explained_variance_ratio": explained,
        "cumulative_explained_variance_ratio": cum_explained,
    }).to_csv(var_path, index=False)

    print(f"[saved] {mat_path}")
    print(f"[saved] {scores_path}")
    print(f"[saved] {loadings_path}")
    print(f"[saved] {var_path}")

    if meta is not None:
        joined = scores_df.copy()
        joined.insert(0, "condition", meta.reindex(scores_df.index).values)
        joined_path = out_dir / "pca_scores_with_condition.csv"
        joined.to_csv(joined_path, index=True)
        print(f"[saved] {joined_path}")


def plot_scores_colored(out_dir: Path, scores_df: pd.DataFrame, cond: pd.Series, mode: str,
                        annotate: bool) -> None:
    """
    Color points by discrete condition values. Unlabeled samples (NaN) shown in gray.
    """
    _lazy_import_matplotlib()
    if not HAVE_MPL:
        print("[warn] matplotlib not available; skipping plots")
        return

    import matplotlib.pyplot as plt  # now safe
    pcs = list(scores_df.columns)
    if len(pcs) < 2 or mode == "none":
        return

    # Align and build masks
    cond = cond.reindex(scores_df.index)
    labeled_mask = cond.notna()
    unlabeled_mask = ~labeled_mask
    cats = pd.Categorical(cond[labeled_mask])
    categories = list(cats.categories)

    # 2D
    if mode in ("2d", "3d"):
        plt.figure(figsize=(7.5, 6.5))
        # plot each category
        for cat in categories:
            m = labeled_mask & (cond == cat)
            plt.scatter(scores_df.loc[m, pcs[0]], scores_df.loc[m, pcs[1]], s=32, alpha=0.9, label=str(cat))
            if annotate:
                for idx in scores_df.loc[m].index:
                    x, y = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]]
                    plt.text(x, y, str(idx), fontsize=8, alpha=0.7)
        # unlabeled
        if unlabeled_mask.any():
            plt.scatter(scores_df.loc[unlabeled_mask, pcs[0]], scores_df.loc[unlabeled_mask, pcs[1]],
                        s=32, alpha=0.5, label="(unlabeled)", marker="x")
            if annotate:
                for idx in scores_df.loc[unlabeled_mask].index:
                    x, y = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]]
                    plt.text(x, y, str(idx), fontsize=8, alpha=0.7)

        plt.xlabel(pcs[0]); plt.ylabel(pcs[1]); plt.title(f"PCA Scores by Condition: {pcs[0]} vs {pcs[1]}")
        plt.grid(True, linestyle="--", alpha=0.35); plt.legend(title="condition", frameon=True)
        plt.tight_layout()
        plt.savefig(out_dir / "pca_scores_2d_by_condition.png", dpi=160)
        plt.close()
        print(f"[saved] {out_dir/'pca_scores_2d_by_condition.png'}")

    # 3D
    if mode == "3d" and len(pcs) >= 3:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        fig = plt.figure(figsize=(8.5, 7.5))
        ax = fig.add_subplot(111, projection="3d")

        for cat in categories:
            m = labeled_mask & (cond == cat)
            ax.scatter(scores_df.loc[m, pcs[0]], scores_df.loc[m, pcs[1]], scores_df.loc[m, pcs[2]],
                       s=32, alpha=0.95, label=str(cat))
            if annotate:
                for idx in scores_df.loc[m].index:
                    x, y, z = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]], scores_df.loc[idx, pcs[2]]
                    ax.text(x, y, z, str(idx), fontsize=8)

        if unlabeled_mask.any():
            ax.scatter(scores_df.loc[unlabeled_mask, pcs[0]], scores_df.loc[unlabeled_mask, pcs[1]],
                       scores_df.loc[unlabeled_mask, pcs[2]], s=32, alpha=0.6, label="(unlabeled)", marker="x")
            if annotate:
                for idx in scores_df.loc[unlabeled_mask].index:
                    x, y, z = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]], scores_df.loc[idx, pcs[2]]
                    ax.text(x, y, z, str(idx), fontsize=8)

        ax.set_xlabel(pcs[0]); ax.set_ylabel(pcs[1]); ax.set_zlabel(pcs[2])
        ax.set_title(f"PCA Scores by Condition: {pcs[0]} vs {pcs[1]} vs {pcs[2]}")
        ax.legend(title="condition")
        fig.tight_layout()
        fig.savefig(out_dir / "pca_scores_3d_by_condition.png", dpi=160)
        print(f"[saved] {out_dir/'pca_scores_3d_by_condition.png'}")


def plot_scree_only(out_dir: Path, explained: np.ndarray) -> None:
    _lazy_import_matplotlib()
    if not HAVE_MPL:
        print("[warn] matplotlib not available; skipping plots")
        return

    import matplotlib.pyplot as plt  # now safe
    plt.figure(figsize=(7, 5))
    xs = np.arange(1, len(explained) + 1)
    plt.bar(xs, explained)
    plt.plot(xs, np.cumsum(explained), marker="o")
    plt.xlabel("Principal Component"); plt.ylabel("Variance Explained Ratio"); plt.title("Scree Plot")
    plt.grid(True, linestyle="--", alpha=0.4); plt.tight_layout()
    plt.savefig(out_dir / "scree_plot.png", dpi=160); plt.close()
    print(f"[saved] {out_dir/'scree_plot.png'}")


def main():
    args = parse_args()
    out_dir = ensure_out_dir(args.output_dir)

    # 1) Build matrix: rows=samples, cols=sites (values = %C or fraction)
    X = build_matrix(args.input_dir, args.glob, args.delimiter, args.as_fraction, args.min_total)

    # 2) Keep sites present in >= min_samples
    X = filter_sites_by_presence(X, args.min_samples)

    # 3) Prepare for PCA (missing + scaling)
    Xp = prepare_for_pca(X, dropna=args.dropna, impute=args.impute, scale=args.scale)

    # 4) PCA via SVD
    scores_df, loadings_df, explained, cum_explained = pca_via_svd(Xp, args.n_components)

    # 5) Save core tables
    meta_series = None
    if args.metadata:
        meta_series = load_metadata_series(
            args.metadata, args.metadata_delimiter, args.metadata_sample_col, args.metadata_condition_col
        )
        # warn about unmatched samples
        unmatched = [s for s in scores_df.index if s not in meta_series.index]
        if unmatched:
            print(f"[warn] {len(unmatched)} samples have no metadata and will be unlabeled in plots:")
            for u in unmatched[:10]:
                print("   -", u)
            if len(unmatched) > 10:
                print("   ...")

    save_tables(out_dir, X, scores_df, loadings_df, explained, cum_explained, meta_series)

    # 6) Plots
    plot_scree_only(out_dir, explained)

    if args.plot != "none":
        if meta_series is not None:
            # Color by condition
            plot_scores_colored(out_dir, scores_df, meta_series, args.plot, args.annotate_samples)
        else:
            # No metadata: simple uncolored plot
            _lazy_import_matplotlib()
            if HAVE_MPL:
                import matplotlib.pyplot as plt
                pcs = list(scores_df.columns)
                if len(pcs) >= 2 and args.plot in ("2d", "3d"):
                    if args.plot == "2d":
                        plt.figure(figsize=(7, 6))
                        plt.scatter(scores_df[pcs[0]], scores_df[pcs[1]], s=28, alpha=0.85)
                        if args.annotate_samples:
                            for idx in scores_df.index:
                                x, y = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]]
                                plt.text(x, y, str(idx), fontsize=8, alpha=0.7)
                        plt.xlabel(pcs[0]); plt.ylabel(pcs[1]); plt.title(f"PCA Scores: {pcs[0]} vs {pcs[1]}")
                        plt.grid(True, linestyle="--", alpha=0.4); plt.tight_layout()
                        plt.savefig(out_dir / "pca_scores_2d.png", dpi=160); plt.close()
                        print(f"[saved] {out_dir/'pca_scores_2d.png'}")
                    else:
                        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
                        fig = plt.figure(figsize=(8, 7)); ax = fig.add_subplot(111, projection="3d")
                        ax.scatter(scores_df[pcs[0]], scores_df[pcs[1]], scores_df[pcs[2]], s=28, alpha=0.9)
                        if args.annotate_samples:
                            for idx in scores_df.index:
                                x, y, z = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]], scores_df.loc[idx, pcs[2]]
                                ax.text(x, y, z, str(idx), fontsize=8)
                        ax.set_xlabel(pcs[0]); ax.set_ylabel(pcs[1]); ax.set_zlabel(pcs[2])
                        ax.set_title(f"PCA Scores: {pcs[0]} vs {pcs[1]} vs {pcs[2]}")
                        fig.tight_layout(); fig.savefig(out_dir / "pca_scores_3d.png", dpi=160)
                        print(f"[saved] {out_dir/'pca_scores_3d.png'}")
            else:
                print("[warn] matplotlib not available; skipping score plot")

    # 7) Console summary
    print("\n[summary]")
    print(f"Samples: {Xp.shape[0]} | Sites: {Xp.shape[1]} | Components: {len(explained)}")
    print("Explained variance ratios:", np.round(explained, 4))
    print("Cumulative explained:", np.round(cum_explained, 4))
    if meta_series is not None:
        counts = meta_series.reindex(scores_df.index).value_counts(dropna=False)
        print("Samples per condition:")
        for k, v in counts.items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
