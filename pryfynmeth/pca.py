#!/usr/bin/env python3
"""
PCA from per-sample site tables using only NumPy/Pandas (no scikit-learn)

- Reads all *.txt (tab-delimited) in --input-dir
- Builds site key "chr:pos:strand" (you may drop strand if desired)
- Value = (C / total) * 100 unless --as-fraction
- Optional metadata: 2 columns [sample, condition]
- Plots (2D/3D) with:
  * PC axis labels showing % variance explained
  * ALWAYS annotated sample labels colored like their points
  * Column-wise label stacking by x (in pixel space)
  * Small x-jitter per stacked label to avoid perfect columns
  * Vertical leader lines point -> label
  * White halo around text for legibility
"""

import argparse
import glob
import os
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd

DEFAULT_PATTERN = "*.txt"

# Optional plotting: force a headless backend if matplotlib is installed
HAVE_MPL = True
def _lazy_import_matplotlib():
    global HAVE_MPL
    try:
        import matplotlib
        matplotlib.use("Agg")             # headless backend for clusters
        import matplotlib.pyplot as plt    # noqa: F401
    except Exception as e:
        HAVE_MPL = False
        print(f"[warn] matplotlib not available; skipping plots ({e})")


def parse_args():
    p = argparse.ArgumentParser(description="Run PCA on C/total percentages across samples (no scikit-learn).")
    # data
    p.add_argument("--input-dir", "-i", required=True, help="Directory with per-sample .txt files (tab-delimited).")
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
    p.add_argument("--label-fontsize", type=int, default=9, help="Font size for sample labels (default: 9).")
    p.add_argument("--axis-fontsize", type=int, default=12, help="Font size for axis titles & ticks (default: 12).")
    # metadata (fixed columns: sample, condition)
    p.add_argument("--metadata", "-m", help="Path to metadata file (TSV/CSV) with columns: sample, condition.")
    p.add_argument("--metadata-delimiter", default="\t", help="Delimiter for metadata file (default: tab).")
    # output
    p.add_argument("--output-dir", "-o", default="pca_numpy_output", help="Where to save outputs.")
    return p.parse_args()


def standardize_headers(df: pd.DataFrame) -> pd.DataFrame:
    mapping = {c: c.strip().lower().replace("-", "_") for c in df.columns}
    return df.rename(columns=mapping)


def sample_name_from_path(path: str) -> str:
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

    # NOTE: if you want strand-collapsed CpGs, change to key = chr:pos
    df["modifier"] = df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":" + df["strand"].astype(str)
    grp = df.groupby("modifier", as_index=False)[["c", "total"]].sum()
    frac = grp["c"] / grp["total"]
    values = frac if as_fraction else frac * 100.0
    return pd.Series(values.values, index=grp["modifier"].values, name=sample_name_from_path(path))


def build_matrix(input_dir: str, sep: str, as_fraction: bool, min_total: int) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(input_dir, DEFAULT_PATTERN)))
    if not files:
        raise FileNotFoundError(f"No files matched {os.path.join(input_dir, DEFAULT_PATTERN)}")
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
    if dropna:
        Xp = X.dropna(axis=0).copy()
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

    if scale:
        means = Xp.mean(axis=0)
        stds = Xp.std(axis=0, ddof=0).replace(0, 1.0)
        Xs = (Xp - means) / stds
    else:
        Xs = Xp - Xp.mean(axis=0)

    return Xs


def pca_via_svd(X: pd.DataFrame, n_components: int):
    X_np = X.values.astype(float, copy=False)
    n_samples, n_features = X_np.shape
    kmax = min(n_samples, n_features)
    if kmax < 2:
        raise ValueError("Need at least 2 samples and 2 sites to run PCA.")

    U, S, Vt = np.linalg.svd(X_np, full_matrices=False)
    k = min(n_components, kmax)

    scores = U[:, :k] * S[:k]
    loadings = Vt[:k, :].T

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


def load_metadata_series(path: str, sep: str) -> pd.Series:
    md = pd.read_csv(path, sep=sep, dtype=str)
    md = standardize_headers(md)
    if "sample" not in md.columns or "condition" not in md.columns:
        raise ValueError("Metadata needs columns 'sample' and 'condition' (case-insensitive).")
    stems = md["sample"].astype(str).apply(lambda x: Path(x).stem)
    cond = md["condition"].astype(str).str.strip()
    meta = pd.Series(cond.values, index=stems, name="condition")
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


# ---------- Label helpers ----------
def _get_scatter_color(sc) -> Tuple[float, float, float, float]:
    fc = sc.get_facecolors()
    if fc is not None and len(fc):
        return tuple(fc[0])
    ec = sc.get_edgecolors()
    if ec is not None and len(ec):
        return tuple(ec[0])
    return (0, 0, 0, 1)


def _shorten_label(name: str) -> str:
    for p in ("Nanopore_", "WGBS_"):
        if name.startswith(p):
            return name[len(p):]
    return name


def _stack_labels_by_x(ax, texts: List["matplotlib.text.Text"], anchors: List[Tuple[float, float]],
                       min_gap_px: float = 4.0, x_bucket_px: float = 18.0,
                       x_jitter_px: float = 6.0, top_margin_px: float = 12.0):
    """
    Deterministic column stacking:
    - Bucket points whose anchor x-values are within ~x_bucket_px in display space.
    - For each bucket, sort by anchor y and place labels in a vertical stack with min_gap_px spacing.
    - Apply small left/right jitter so labels don't form a perfect column.
    - Cap label positions below the top of the axes by top_margin_px.
    """
    import matplotlib.pyplot as plt  # noqa: F401
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    trans = ax.transData
    inv = trans.inverted()

    # anchors in display coords
    pts_disp = trans.transform(np.array(anchors))
    xs = pts_disp[:, 0]
    ys = pts_disp[:, 1]

    # bucket by rounded pixel x
    buckets = {}
    for i, xpix in enumerate(xs):
        key = int(round(xpix / x_bucket_px))
        buckets.setdefault(key, []).append(i)

    ax_bbox = ax.get_window_extent(renderer)
    y_cap = ax_bbox.y1 - top_margin_px  # do not cross plot top

    for key, idxs in buckets.items():
        # sort by anchor y (ascending)
        idxs.sort(key=lambda i: ys[i])

        # get label heights in pixels
        heights = []
        for i in idxs:
            # Use a tiny fudge to ensure bbox exists
            bb = texts[i].get_window_extent(renderer).expanded(1.0, 1.0)
            heights.append(bb.height)

        # stack
        prev_top = None
        for rank, i in enumerate(idxs):
            xpix, ypix = xs[i], ys[i]
            # baseline: just above the point by min_gap_px
            target_bottom = ypix + min_gap_px
            if prev_top is None:
                bottom = target_bottom
            else:
                bottom = max(target_bottom, prev_top + min_gap_px)
            top = bottom + heights[rank]
            # cap at top of axes
            if top > y_cap:
                top = y_cap
                bottom = top - heights[rank]

            # add small left/right jitter by rank
            jitter = ((-1) ** rank) * (x_jitter_px * (rank // 2))

            x_data, y_data = inv.transform((xpix + jitter, bottom))
            texts[i].set_position((x_data, y_data))
            prev_top = top


# ---------- Plotting ----------
def plot_scores_colored(out_dir: Path, scores_df: pd.DataFrame, cond: pd.Series, mode: str,
                        explained: np.ndarray, label_fontsize: int, axis_fontsize: int) -> None:
    _lazy_import_matplotlib()
    if not HAVE_MPL:
        print("[warn] matplotlib not available; skipping plots")
        return

    import matplotlib.pyplot as plt  # now safe
    import matplotlib.patheffects as pe

    pcs = list(scores_df.columns)
    if len(pcs) < 2 or mode == "none":
        return

    cond = cond.reindex(scores_df.index)
    labeled_mask = cond.notna()
    unlabeled_mask = ~labeled_mask
    categories = list(pd.Categorical(cond[labeled_mask]).categories)

    pct1 = f"{(explained[0] * 100):.1f}%" if len(explained) >= 1 else ""
    pct2 = f"{(explained[1] * 100):.1f}%" if len(explained) >= 2 else ""

    # 2D
    if mode in ("2d", "3d"):
        plt.figure(figsize=(9.5, 7.5))
        ax = plt.gca()

        texts: List["matplotlib.text.Text"] = []
        anchors: List[Tuple[float, float]] = []
        colors: List[Tuple[float, float, float, float]] = []

        def _add_group(mask, legend_label):
            sc = ax.scatter(scores_df.loc[mask, pcs[0]], scores_df.loc[mask, pcs[1]],
                            s=36, alpha=0.9, label=str(legend_label))
            col = _get_scatter_color(sc)
            for idx in scores_df.loc[mask].index:
                x = scores_df.loc[idx, pcs[0]]
                y = scores_df.loc[idx, pcs[1]]
                label = _shorten_label(str(idx))
                # start label slightly above point; exact stack happens later
                t = ax.text(x, y, label,
                            fontsize=label_fontsize, color=col, ha="center", va="bottom",
                            path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
                texts.append(t); anchors.append((x, y)); colors.append(col)

        for cat in categories:
            _add_group(labeled_mask & (cond == cat), cat)
        if unlabeled_mask.any():
            _add_group(unlabeled_mask, "(unlabeled)")

        # column stacking by x (in pixel space)
        if texts:
            _stack_labels_by_x(ax, texts, anchors,
                               min_gap_px=4.0, x_bucket_px=18.0,
                               x_jitter_px=6.0, top_margin_px=12.0)

        # leader lines from point -> final label position
        for (x, y), t, col in zip(anchors, texts, colors):
            x2, y2 = t.get_position()
            ax.plot([x, x2], [y, y2], lw=0.8, alpha=0.7, color=col, zorder=0)

        ax.set_xlabel(f"{pcs[0]} ({pct1})", fontsize=axis_fontsize)
        ax.set_ylabel(f"{pcs[1]} ({pct2})", fontsize=axis_fontsize)
        ax.set_title(f"PCA Scores by Condition: {pcs[0]} vs {pcs[1]}", fontsize=axis_fontsize + 1)
        ax.grid(True, linestyle="--", alpha=0.35)
        ax.tick_params(labelsize=max(axis_fontsize - 1, 1))
        ax.legend(title="condition", frameon=True)
        plt.tight_layout()
        plt.savefig(out_dir / "pca_scores_2d_by_condition.png", dpi=160)
        plt.close()
        print(f"[saved] {out_dir/'pca_scores_2d_by_condition.png'}")

    # 3D (unchanged apart from centered labels)
    if mode == "3d" and len(pcs) >= 3:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        fig = plt.figure(figsize=(9, 7.6))
        ax = fig.add_subplot(111, projection="3d")

        for cat in categories:
            m = labeled_mask & (cond == cat)
            sc = ax.scatter(scores_df.loc[m, pcs[0]], scores_df.loc[m, pcs[1]], scores_df.loc[m, pcs[2]],
                            s=32, alpha=0.95, label=str(cat))
            col = _get_scatter_color(sc)
            for idx in scores_df.loc[m].index:
                x, y, z = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]], scores_df.loc[idx, pcs[2]]
                ax.text(x, y, z, str(idx), fontsize=label_fontsize, color=col, ha="center")

        if unlabeled_mask.any():
            sc = ax.scatter(scores_df.loc[unlabeled_mask, pcs[0]], scores_df.loc[unlabeled_mask, pcs[1]],
                            scores_df.loc[unlabeled_mask, pcs[2]], s=32, alpha=0.6, label="(unlabeled)", marker="x")
            col = _get_scatter_color(sc)
            for idx in scores_df.loc[unlabeled_mask].index:
                x, y, z = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]], scores_df.loc[idx, pcs[2]]
                ax.text(x, y, z, str(idx), fontsize=label_fontsize, color=col, ha="center")

        ax.set_xlabel(f"{pcs[0]} ({pct1})", fontsize=axis_fontsize)
        ax.set_ylabel(f"{pcs[1]} ({pct2})", fontsize=axis_fontsize)
        if len(explained) >= 3:
            pct3 = f"{(explained[2] * 100):.1f}%"
            ax.set_zlabel(f"{pcs[2]} ({pct3})", fontsize=axis_fontsize)
        ax.set_title(f"PCA Scores by Condition: {pcs[0]} vs {pcs[1]} vs {pcs[2]}", fontsize=axis_fontsize + 1)
        ax.legend(title="condition")
        plt.tight_layout()
        fig.savefig(out_dir / "pca_scores_3d_by_condition.png", dpi=160)
        print(f"[saved] {out_dir/'pca_scores_3d_by_condition.png'}")


def plot_scree_only(out_dir: Path, explained: np.ndarray, axis_fontsize: int) -> None:
    _lazy_import_matplotlib()
    if not HAVE_MPL:
        print("[warn] matplotlib not available; skipping plots")
        return

    import matplotlib.pyplot as plt  # now safe
    plt.figure(figsize=(7.6, 5.4))
    xs = np.arange(1, len(explained) + 1)
    plt.bar(xs, explained)
    plt.plot(xs, np.cumsum(explained), marker="o")
    plt.xlabel("Principal Component", fontsize=axis_fontsize)
    plt.ylabel("Variance Explained Ratio", fontsize=axis_fontsize)
    plt.title("Scree Plot", fontsize=axis_fontsize + 1)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tick_params(labelsize=max(axis_fontsize - 1, 1))
    plt.tight_layout()
    plt.savefig(out_dir / "scree_plot.png", dpi=160)
    plt.close()
    print(f"[saved] {out_dir/'scree_plot.png'}")


def main():
    args = parse_args()
    out_dir = ensure_out_dir(args.output_dir)

    # 1) Build matrix: rows=samples, cols=sites (values = %C or fraction)
    X = build_matrix(args.input_dir, args.delimiter, args.as_fraction, args.min_total)

    # 2) Keep sites present in >= min_samples
    X = filter_sites_by_presence(X, args.min_samples)

    # 3) Prepare for PCA (missing + scaling)
    Xp = prepare_for_pca(X, dropna=args.dropna, impute=args.impute, scale=args.scale)

    # 4) PCA via SVD
    scores_df, loadings_df, explained, cum_explained = pca_via_svd(Xp, args.n_components)

    # 5) Save core tables (+ metadata if provided)
    meta_series = None
    if args.metadata:
        meta_series = load_metadata_series(args.metadata, args.metadata_delimiter)
        unmatched = [s for s in scores_df.index if s not in meta_series.index]
        if unmatched:
            print(f"[warn] {len(unmatched)} samples have no metadata and will be unlabeled in plots:")
            for u in unmatched[:10]:
                print("   -", u)
            if len(unmatched) > 10:
                print("   ...")

    save_tables(out_dir, X, scores_df, loadings_df, explained, cum_explained, meta_series)

    # 6) Plots
    plot_scree_only(out_dir, explained, args.axis_fontsize)
    if args.plot != "none":
        if meta_series is not None:
            plot_scores_colored(out_dir, scores_df, meta_series, args.plot,
                                explained, args.label_fontsize, args.axis_fontsize)
        else:
            # simple (no metadata) plot with same label treatment
            _lazy_import_matplotlib()
            if HAVE_MPL:
                import matplotlib.pyplot as plt
                import matplotlib.patheeffects as pe
                pcs = list(scores_df.columns)
                if len(pcs) >= 2 and args.plot in ("2d", "3d"):
                    pct1 = f"{(explained[0] * 100):.1f}%" if len(explained) >= 1 else ""
                    pct2 = f"{(explained[1] * 100):.1f}%" if len(explained) >= 2 else ""
                    if args.plot == "2d":
                        plt.figure(figsize=(9.5, 7.5))
                        ax = plt.gca()
                        sc = ax.scatter(scores_df[pcs[0]], scores_df[pcs[1]], s=36, alpha=0.9)
                        col = _get_scatter_color(sc)
                        texts = []; anchors = []
                        for idx in scores_df.index:
                            x, y = scores_df.loc[idx, pcs[0]], scores_df.loc[idx, pcs[1]]
                            t = ax.text(x, y, str(idx), fontsize=args.label_fontsize, color=col,
                                        ha="center", va="bottom",
                                        path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
                            texts.append(t); anchors.append((x, y))
                        if texts:
                            _stack_labels_by_x(ax, texts, anchors)
                        for (x, y), t in zip(anchors, texts):
                            x2, y2 = t.get_position()
                            ax.plot([x, x2], [y, y2], lw=0.8, alpha=0.7, color=col, zorder=0)
                        ax.set_xlabel(f"{pcs[0]} ({pct1})", fontsize=args.axis_fontsize)
                        ax.set_ylabel(f"{pcs[1]} ({pct2})", fontsize=args.axis_fontsize)
                        ax.set_title(f"PCA Scores: {pcs[0]} vs {pcs[1]}", fontsize=args.axis_fontsize + 1)
                        ax.grid(True, linestyle="--", alpha=0.4)
                        ax.tick_params(labelsize=max(args.axis_fontsize - 1, 1))
                        plt.tight_layout()
                        plt.savefig(out_dir / "pca_scores_2d.png", dpi=160)
                        plt.close()
                        print(f"[saved] {out_dir/'pca_scores_2d.png'}")
            else:
                print("[warn] matplotlib not available; skipping score plot")

    # 7) Console summary
    print("\n[summary]")
    print(f"Samples: {Xp.shape[0]} | Sites: {Xp.shape[1]} | Components: {len(explained)}")
    print("Explained variance ratios:", np.round(explained, 4))
    print("Cumulative explained:", np.round(cum_explained, 4))


if __name__ == "__main__":
    main()
