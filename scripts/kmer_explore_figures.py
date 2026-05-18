#!/usr/bin/env python3
"""Create figures from BEST k-mer exploration outputs.

The figures are a presentation layer on top of ``scripts/kmer_explore.py``. They
do not read BAMs or references; they consume the derived CSVs in an exploration
directory and write portable PNG/SVG figures plus a Markdown manifest.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import defaultdict
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, Normalize

try:
    from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
except ImportError:  # pragma: no cover - exercised only when SciPy is absent.
    dendrogram = None
    leaves_list = None
    linkage = None


BASES = ("A", "C", "G", "T")
BASE_COLORS = {
    "A": "#3B78E7",
    "C": "#6A994E",
    "G": "#8E5CF3",
    "T": "#E07A5F",
}
BASE_TO_INDEX = {base: idx for idx, base in enumerate(BASES)}
ERROR_CLASS_COLUMNS = [
    ("mismatches_per_interval", "Mismatch", "#3B78E7"),
    ("non_hp_ins_per_interval", "Non-HP insertion", "#E07A5F"),
    ("non_hp_del_per_interval", "Non-HP deletion", "#6A994E"),
    ("hp_ins_per_interval", "HP insertion", "#F2CC8F"),
    ("hp_del_per_interval", "HP deletion", "#8E5CF3"),
]
SUBSTITUTIONS = [
    f"{ref}>{read}"
    for ref in BASES
    for read in BASES
    if ref != read
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--explore-dir",
        required=True,
        help="Directory created by scripts/kmer_explore.py",
    )
    parser.add_argument(
        "--out-dir",
        help="Output directory for figures. Defaults to <explore-dir>/figures",
    )
    parser.add_argument(
        "--label",
        default="Nanopore k-mer error profile",
        help="Dataset label used in figure titles and the Markdown manifest",
    )
    parser.add_argument(
        "--platform",
        default="Nanopore",
        help="Sequencing platform label used in the Markdown manifest",
    )
    parser.add_argument(
        "--figure-prefix",
        help="Prefix for output figure filenames. Defaults to a slug of --platform.",
    )
    parser.add_argument(
        "--min-intervals",
        type=int,
        default=1000,
        help="Minimum support for point/bar inclusion when a CSV has interval counts",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of contexts or motifs shown in ranked figures",
    )
    parser.add_argument(
        "--annotate-top",
        type=int,
        default=8,
        help="Number of top contexts annotated in scatter plots",
    )
    parser.add_argument(
        "--formats",
        default="png,svg",
        help="Comma-separated output formats supported by Matplotlib",
    )
    parser.add_argument(
        "--max-run-length",
        type=int,
        default=12,
        help="Maximum homopolymer run length to show in the run-profile figure",
    )
    args = parser.parse_args()
    if args.min_intervals < 0:
        parser.error("--min-intervals must be nonnegative")
    if args.top_n <= 0:
        parser.error("--top-n must be positive")
    if args.annotate_top < 0:
        parser.error("--annotate-top must be nonnegative")
    args.formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    if not args.formats:
        parser.error("--formats must include at least one format")
    args.figure_prefix = args.figure_prefix or slugify(args.platform)
    return args


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "kmer"


def figure_stem(file_prefix: str, name: str) -> str:
    return f"{file_prefix}_{name}"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    return default if value == "" else float(value)


def i(row: dict[str, str], key: str, default: int = 0) -> int:
    value = row.get(key, "")
    return default if value == "" else int(float(value))


def passes_support(row: dict[str, str], min_intervals: int, key: str = "intervals") -> bool:
    return i(row, key) >= min_intervals


def percent(value: float) -> float:
    return 100.0 * value


def configure_axes(ax, *, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="x", color="#D8DEE9", linewidth=0.6, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig, out_dir: Path, stem: str, formats: list[str]) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt}"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def choose_supported(rows: list[dict[str, str]], min_intervals: int) -> list[dict[str, str]]:
    supported = [row for row in rows if passes_support(row, min_intervals)]
    return supported or rows


def top_supported_risk_rows(
    rows: list[dict[str, str]],
    min_intervals: int,
    top_n: int,
) -> list[dict[str, str]]:
    rows = choose_supported(rows, min_intervals)
    rows = sorted(
        rows,
        key=lambda row: (
            -f(row, "risk_score", f(row, "error_rate")),
            -f(row, "error_rate"),
            -i(row, "intervals"),
        ),
    )
    return rows[:top_n]


def dominant_error_label(row: dict[str, str]) -> str:
    mismatch_fraction = f(row, "mismatch_fraction")
    non_hp_indel_fraction = f(row, "non_hp_indel_fraction")
    indel_fraction = f(row, "indel_fraction")
    if mismatch_fraction >= max(non_hp_indel_fraction, indel_fraction - non_hp_indel_fraction):
        return "Mismatch-heavy"
    if non_hp_indel_fraction >= 0.5:
        return "Non-HP indel-heavy"
    return "HP/other indel-heavy"


def zscore_columns(matrix: np.ndarray) -> np.ndarray:
    if matrix.size == 0:
        return matrix
    means = matrix.mean(axis=0)
    stds = matrix.std(axis=0)
    stds[stds == 0.0] = 1.0
    return (matrix - means) / stds


def cluster_order(matrix: np.ndarray) -> tuple[list[int], object | None]:
    if len(matrix) <= 2 or linkage is None or leaves_list is None:
        return list(range(len(matrix))), None
    model = linkage(matrix, method="average", metric="euclidean")
    return [int(idx) for idx in leaves_list(model)], model


def draw_base_sequence_heatmap(
    ax,
    kmers: list[str],
    *,
    x_labels: list[str] | None = None,
    title: str | None = None,
) -> None:
    matrix = np.array(
        [[BASE_TO_INDEX.get(base, np.nan) for base in kmer] for kmer in kmers],
        dtype=float,
    )
    cmap = ListedColormap([BASE_COLORS[base] for base in BASES])
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    ax.set_xticks(list(range(matrix.shape[1])), x_labels or list(range(matrix.shape[1])))
    ax.set_yticks(list(range(len(kmers))), kmers)
    if title:
        ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.set_xlabel("Position in read-oriented k-mer")
    ax.set_ylabel("K-mer")
    ax.tick_params(axis="both", length=0)
    for y_idx, kmer in enumerate(kmers):
        for x_idx, base in enumerate(kmer):
            if base in BASE_TO_INDEX:
                ax.text(x_idx, y_idx, base, ha="center", va="center", fontsize=8, color="white")


def top_motif(rows: list[dict[str, str]]) -> tuple[str, int | None] | None:
    for row in rows:
        motif = row.get("motif", "")
        if not motif:
            continue
        motif_type = row.get("motif_type", "")
        if motif_type.startswith("offset_"):
            try:
                return motif, int(motif_type.split("_", 1)[1])
            except ValueError:
                return motif, None
        return motif, None
    return None


def plot_context_error_classes(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    rows = choose_supported(rows, min_intervals)[:top_n]
    if not rows:
        return None

    rows = list(reversed(rows))
    labels = [row["canonical_kmer"] for row in rows]
    y_positions = list(range(len(rows)))
    fig, ax = plt.subplots(figsize=(9.0, max(4.5, 0.35 * len(rows) + 1.5)))
    left = [0.0] * len(rows)
    for column, legend, color in ERROR_CLASS_COLUMNS:
        values = [f(row, column) for row in rows]
        ax.barh(y_positions, values, left=left, label=legend, color=color, height=0.72)
        left = [prev + value for prev, value in zip(left, values)]

    ax.set_yticks(y_positions, labels)
    configure_axes(
        ax,
        title=f"{label}: enriched reverse-complement-collapsed contexts",
        xlabel="Error events per k-mer interval",
        ylabel="Canonical k-mer",
    )
    max_total = max(left) if left else 0.0
    if max_total:
        ax.set_xlim(right=max_total * 1.08)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=9)
    for idx, row in enumerate(rows):
        total = sum(f(row, column) for column, _, _ in ERROR_CLASS_COLUMNS)
        ax.text(
            total + max_total * 0.01,
            idx,
            f"{f(row, 'enrichment'):.2f}x",
            va="center",
            fontsize=8,
            color="#344054",
        )

    return (
        "context_error_classes",
        save_figure(fig, out_dir, figure_stem(file_prefix, "context_error_classes"), formats),
    )


def plot_top_kmer_lollipop(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    rows = top_supported_risk_rows(rows, min_intervals, top_n)
    if not rows:
        return None

    rows = list(reversed(rows))
    labels = [row["kmer"] for row in rows]
    values = [percent(f(row, "error_rate")) for row in rows]
    categories = [dominant_error_label(row) for row in rows]
    colors = {
        "Mismatch-heavy": "#3B78E7",
        "Non-HP indel-heavy": "#E07A5F",
        "HP/other indel-heavy": "#8E5CF3",
    }
    y_positions = list(range(len(rows)))
    fig, ax = plt.subplots(figsize=(8.2, max(4.8, 0.34 * len(rows) + 1.6)))
    for y_pos, value, category in zip(y_positions, values, categories):
        ax.hlines(y_pos, 0, value, color="#D0D5DD", linewidth=1.5)
        ax.scatter(value, y_pos, s=70, color=colors[category], label=category, zorder=3)
    ax.set_yticks(y_positions, labels)
    configure_axes(
        ax,
        title=f"{label}: top high-risk k-mers",
        xlabel="Error rate (%)",
        ylabel="K-mer",
    )
    max_value = max(values) if values else 0.0
    if max_value:
        ax.set_xlim(0, max_value * 1.15)
    for y_pos, row, value in zip(y_positions, rows, values):
        ax.text(
            value + max_value * 0.015,
            y_pos,
            f"risk {f(row, 'risk_score'):.2f}",
            va="center",
            fontsize=8,
            color="#344054",
        )

    handles, labels_seen = ax.get_legend_handles_labels()
    unique = dict(zip(labels_seen, handles))
    ax.legend(unique.values(), unique.keys(), loc="lower right", frameon=False, fontsize=8)
    return (
        "top_kmer_lollipop",
        save_figure(fig, out_dir, figure_stem(file_prefix, "top_kmer_lollipop"), formats),
    )


def plot_top_kmer_clustered_heatmap(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    rows = top_supported_risk_rows(rows, min_intervals, top_n)
    if len(rows) < 2:
        return None

    feature_names = [
        "Error rate",
        "Risk",
        "Mismatch frac",
        "Indel frac",
        "Non-HP indel frac",
        "log10 intervals",
    ]
    feature_matrix = np.array(
        [
            [
                f(row, "error_rate"),
                f(row, "risk_score"),
                f(row, "mismatch_fraction"),
                f(row, "indel_fraction"),
                f(row, "non_hp_indel_fraction"),
                math.log10(i(row, "intervals") + 1),
            ]
            for row in rows
        ],
        dtype=float,
    )
    feature_z = zscore_columns(feature_matrix)
    seq_matrix = []
    for row in rows:
        encoded = []
        for base in row["kmer"]:
            encoded.extend(1.0 if base == expected else 0.0 for expected in BASES)
        seq_matrix.append(encoded)
    cluster_matrix = np.concatenate([feature_z, 0.5 * np.array(seq_matrix, dtype=float)], axis=1)
    order, model = cluster_order(cluster_matrix)
    if model is not None and dendrogram is not None:
        dendro_info = dendrogram(model, no_plot=True)
        order = list(reversed([int(idx) for idx in dendro_info["leaves"]]))

    ordered_rows = [rows[idx] for idx in order]
    ordered_features = feature_z[order, :]
    fig = plt.figure(figsize=(10.2, max(5.0, 0.34 * len(rows) + 1.8)))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.25, 4.8, 0.22], wspace=0.06)
    ax_tree = fig.add_subplot(grid[0, 0])
    ax_heatmap = fig.add_subplot(grid[0, 1])
    ax_cbar = fig.add_subplot(grid[0, 2])

    if model is not None and dendrogram is not None:
        dendrogram(
            model,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#667085",
            ax=ax_tree,
        )
        ax_tree.invert_yaxis()
    ax_tree.axis("off")

    image = ax_heatmap.imshow(ordered_features, aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
    ax_heatmap.set_title(
        f"{label}: top k-mer error-profile clusters",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax_heatmap.set_xticks(list(range(len(feature_names))), feature_names, rotation=35, ha="right")
    ax_heatmap.set_yticks(list(range(len(ordered_rows))), [row["kmer"] for row in ordered_rows])
    ax_heatmap.tick_params(axis="both", length=0)
    cbar = fig.colorbar(image, cax=ax_cbar)
    cbar.set_label("Column z-score")

    return (
        "top_kmer_clustered_heatmap",
        save_figure(fig, out_dir, figure_stem(file_prefix, "top_kmer_clustered_heatmap"), formats),
    )


def plot_top_kmer_sequence_heatmap(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    rows = top_supported_risk_rows(rows, min_intervals, top_n)
    if not rows:
        return None

    kmers = [row["kmer"] for row in rows]
    values = [percent(f(row, "error_rate")) for row in rows]
    fig = plt.figure(figsize=(9.0, max(4.8, 0.34 * len(rows) + 1.6)))
    grid = fig.add_gridspec(1, 2, width_ratios=[4.8, 1.5], wspace=0.06)
    ax_seq = fig.add_subplot(grid[0, 0])
    ax_bar = fig.add_subplot(grid[0, 1], sharey=ax_seq)
    draw_base_sequence_heatmap(
        ax_seq,
        kmers,
        title=f"{label}: top k-mer sequence composition",
    )
    ax_bar.barh(list(range(len(rows))), values, color="#3B78E7", height=0.72)
    ax_bar.set_xlabel("Error rate (%)")
    ax_bar.grid(True, axis="x", color="#D8DEE9", linewidth=0.6, alpha=0.8)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.tick_params(axis="y", left=False, labelleft=False)
    ax_bar.set_ylim(len(rows) - 0.5, -0.5)
    return (
        "top_kmer_sequence_heatmap",
        save_figure(fig, out_dir, figure_stem(file_prefix, "top_kmer_sequence_heatmap"), formats),
    )


def plot_motif_aligned_contexts(
    risk_rows: list[dict[str, str]],
    motif_rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    motif_info = top_motif(motif_rows)
    if motif_info is None:
        return None
    motif, fixed_offset = motif_info
    rows = []
    for row in top_supported_risk_rows(risk_rows, min_intervals, len(risk_rows)):
        kmer = row["kmer"]
        if fixed_offset is not None:
            if kmer[fixed_offset : fixed_offset + len(motif)] == motif:
                rows.append((row, fixed_offset))
        else:
            pos = kmer.find(motif)
            if pos >= 0:
                rows.append((row, pos))
        if len(rows) >= top_n:
            break
    if not rows:
        return None

    max_left = max(pos for _, pos in rows)
    max_right = max(len(row["kmer"]) - (pos + len(motif)) for row, pos in rows)
    width = max_left + len(motif) + max_right
    matrix = np.full((len(rows), width), np.nan)
    labels = []
    placed: list[tuple[int, int, str]] = []
    for y_idx, (row, pos) in enumerate(rows):
        start = max_left - pos
        labels.append(f"{row['kmer']}  {percent(f(row, 'error_rate')):.1f}%")
        for x_idx, base in enumerate(row["kmer"]):
            out_x = start + x_idx
            matrix[y_idx, out_x] = BASE_TO_INDEX.get(base, np.nan)
            placed.append((y_idx, out_x, base))

    cmap = ListedColormap([BASE_COLORS[base] for base in BASES])
    cmap.set_bad("#FFFFFF")
    x_labels = [str(idx - max_left) for idx in range(width)]
    fig, ax = plt.subplots(figsize=(8.8, max(4.8, 0.34 * len(rows) + 1.6)))
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    ax.axvspan(max_left - 0.5, max_left + len(motif) - 0.5, color="#F2CC8F", alpha=0.28)
    ax.set_title(
        f"{label}: top contexts aligned on `{motif}`",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_xticks(list(range(width)), x_labels)
    ax.set_yticks(list(range(len(rows))), labels)
    ax.set_xlabel("Position relative to motif start")
    ax.set_ylabel("K-mer and error rate")
    ax.tick_params(axis="both", length=0)
    for y_idx, x_idx, base in placed:
        if base in BASE_TO_INDEX:
            ax.text(x_idx, y_idx, base, ha="center", va="center", fontsize=8, color="white")

    return (
        "motif_aligned_contexts",
        save_figure(fig, out_dir, figure_stem(file_prefix, "motif_aligned_contexts"), formats),
    )


def plot_rc_pair_asymmetry(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    annotate_top: int,
) -> tuple[str, list[Path]] | None:
    rows = [
        row
        for row in rows
        if row.get("pair_status") == "complete"
        and passes_support(row, min_intervals, "min_member_intervals")
    ]
    if not rows:
        return None

    x_values = [percent(f(row, "member_a_error_rate")) for row in rows]
    y_values = [percent(f(row, "member_b_error_rate")) for row in rows]
    deltas = [percent(f(row, "delta_b_minus_a")) for row in rows]
    sizes = [max(12.0, 10.0 * math.log10(i(row, "min_member_intervals") + 1)) for row in rows]
    limit = max(x_values + y_values) * 1.06

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    norm = Normalize(vmin=-max(abs(value) for value in deltas), vmax=max(abs(value) for value in deltas))
    scatter = ax.scatter(
        x_values,
        y_values,
        c=deltas,
        s=sizes,
        cmap="coolwarm",
        norm=norm,
        alpha=0.55,
        edgecolors="none",
        rasterized=True,
    )
    ax.plot([0, limit], [0, limit], color="#475467", linestyle="--", linewidth=1)
    configure_axes(
        ax,
        title=f"{label}: k-mer vs reverse-complement error rates",
        xlabel="Canonical member error rate (%)",
        ylabel="Reverse-complement member error rate (%)",
    )
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.82)
    cbar.set_label("Reverse-complement minus canonical error rate (percentage points)")

    for row in rows[:annotate_top]:
        x_value = percent(f(row, "member_a_error_rate"))
        y_value = percent(f(row, "member_b_error_rate"))
        ax.annotate(
            f"{row['member_a']}/{row['member_b']}",
            (x_value, y_value),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
            color="#111827",
        )

    return (
        "rc_pair_asymmetry",
        save_figure(fig, out_dir, figure_stem(file_prefix, "rc_pair_asymmetry"), formats),
    )


def plot_rc_strand_mirror(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    rows = [
        row
        for row in rows
        if passes_support(row, min_intervals, "min_member_strand_intervals")
    ][:top_n]
    if not rows:
        return None

    rows = list(reversed(rows))
    labels = [f"{row['member_a']}/{row['member_b']}" for row in rows]
    values = [percent(f(row, "mirror_discordance")) for row in rows]
    within = [percent(f(row, "max_within_member_delta")) for row in rows]
    y_positions = list(range(len(rows)))

    fig, ax = plt.subplots(figsize=(9.0, max(4.5, 0.35 * len(rows) + 1.5)))
    bars = ax.barh(y_positions, values, color="#3B78E7", height=0.72, alpha=0.86)
    ax.scatter(within, y_positions, color="#D92D20", s=28, label="Max strand delta")
    ax.set_yticks(y_positions, labels)
    configure_axes(
        ax,
        title=f"{label}: reverse-complement mirror diagnostics",
        xlabel="Error-rate delta or mirror discordance (percentage points)",
        ylabel="K-mer / reverse complement",
    )
    ax.legend(loc="lower right", frameon=False)
    for bar, row in zip(bars, rows):
        ax.text(
            bar.get_width() + max(values) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"n>={row['min_member_strand_intervals']}",
            va="center",
            fontsize=8,
            color="#344054",
        )

    return (
        "rc_strand_mirror",
        save_figure(fig, out_dir, figure_stem(file_prefix, "rc_strand_mirror"), formats),
    )


def plot_quality_calibration(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    min_intervals: int,
    annotate_top: int,
) -> tuple[str, list[Path]] | None:
    rows = [row for row in rows if i(row, "intervals") >= min_intervals]
    if not rows:
        return None

    x_values = [f(row, "mean_qual") for row in rows]
    y_values = [f(row, "empirical_qv") for row in rows]
    overconfidence = [f(row, "overconfidence_qv") for row in rows]
    sizes = [max(10.0, 7.0 * math.log10(i(row, "intervals") + 1)) for row in rows]
    axis_min = max(0.0, min(x_values + y_values) - 1.0)
    axis_max = max(x_values + y_values) + 1.0

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    scatter = ax.scatter(
        x_values,
        y_values,
        c=overconfidence,
        s=sizes,
        cmap="viridis",
        alpha=0.55,
        edgecolors="none",
        rasterized=True,
    )
    ax.plot([axis_min, axis_max], [axis_min, axis_max], color="#475467", linestyle="--")
    configure_axes(
        ax,
        title=f"{label}: context-specific quality calibration",
        xlabel="Mean predicted QV",
        ylabel="Empirical QV",
    )
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.82)
    cbar.set_label("Predicted minus empirical QV")

    rows_by_overconfidence = sorted(rows, key=lambda row: -f(row, "overconfidence_qv"))
    for row in rows_by_overconfidence[:annotate_top]:
        ax.annotate(
            row["kmer"],
            (f(row, "mean_qual"), f(row, "empirical_qv")),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
            color="#111827",
        )

    return (
        "quality_calibration",
        save_figure(fig, out_dir, figure_stem(file_prefix, "quality_calibration"), formats),
    )


def plot_homopolymer_profile(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    max_run_length: int,
) -> tuple[str, list[Path]] | None:
    if not rows:
        return None

    grouped: dict[tuple[str, int], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in rows:
        run_length = i(row, "run_length_within_kmer")
        if run_length > max_run_length:
            continue
        key = (row["base"], run_length)
        grouped[key]["matches"] += i(row, "matches")
        for column in ["mismatches", "non_hp_ins", "non_hp_del", "hp_ins", "hp_del"]:
            grouped[key]["errors"] += i(row, column)

    if not grouped:
        return None

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    colors = {"A": "#3B78E7", "C": "#6A994E", "G": "#8E5CF3", "T": "#E07A5F"}
    for base in BASES:
        base_points = []
        for (base_name, run_length), values in grouped.items():
            if base_name != base:
                continue
            observations = values["matches"] + values["errors"]
            if observations:
                base_points.append((run_length, percent(values["errors"] / observations)))
        if not base_points:
            continue
        base_points.sort()
        ax.plot(
            [point[0] for point in base_points],
            [point[1] for point in base_points],
            marker="o",
            linewidth=2,
            label=base,
            color=colors[base],
        )

    configure_axes(
        ax,
        title=f"{label}: homopolymer run-length error profile",
        xlabel="Run length within k-mer",
        ylabel="Error rate (%)",
    )
    ax.legend(title="Base", frameon=False)
    return (
        "homopolymer_profile",
        save_figure(fig, out_dir, figure_stem(file_prefix, "homopolymer_run_profile"), formats),
    )


def plot_motif_enrichment(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
    top_n: int,
) -> tuple[str, list[Path]] | None:
    rows = [row for row in rows if i(row, "bad_with_motif") > 0][:top_n]
    if not rows:
        return None

    rows = list(reversed(rows))
    labels = [
        f"{row['motif_type'].replace('_', ' ')}: {row['motif']}"
        for row in rows
    ]
    values = [f(row, "odds_ratio") for row in rows]
    y_positions = list(range(len(rows)))

    fig, ax = plt.subplots(figsize=(8.4, max(4.5, 0.34 * len(rows) + 1.5)))
    ax.barh(y_positions, values, color="#6A994E", height=0.72)
    ax.axvline(1.0, color="#475467", linestyle="--", linewidth=1)
    ax.set_yticks(y_positions, labels)
    configure_axes(
        ax,
        title=f"{label}: model-free motif enrichment among high-error k-mers",
        xlabel="Odds ratio among high-error k-mers",
        ylabel="Motif",
    )
    return (
        "motif_enrichment",
        save_figure(fig, out_dir, figure_stem(file_prefix, "motif_enrichment"), formats),
    )


def plot_offset_base_heatmap(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
) -> tuple[str, list[Path]] | None:
    rows = [row for row in rows if row.get("kmer_base") in BASES]
    if not rows:
        return None

    offsets = sorted({i(row, "offset") for row in rows})
    values: dict[tuple[str, int], float] = defaultdict(float)
    for row in rows:
        values[(row["kmer_base"], i(row, "offset"))] = percent(f(row, "error_rate"))

    matrix = [[values[(base, offset)] for offset in offsets] for base in BASES]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    image = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_xticks(list(range(len(offsets))), offsets)
    ax.set_yticks(list(range(len(BASES))), BASES)
    ax.set_title(
        f"{label}: base- and offset-specific error rate",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_xlabel("Offset within read-oriented k-mer")
    ax.set_ylabel("Reference base")
    cbar = fig.colorbar(image, ax=ax, shrink=0.82)
    cbar.set_label("Error rate (%)")
    return (
        "offset_base_heatmap",
        save_figure(fig, out_dir, figure_stem(file_prefix, "offset_base_error_heatmap"), formats),
    )


def plot_substitution_spectrum(
    rows: list[dict[str, str]],
    out_dir: Path,
    formats: list[str],
    label: str,
    file_prefix: str,
) -> tuple[str, list[Path]] | None:
    if not rows:
        return None

    counts: dict[str, float] = defaultdict(float)
    intervals: dict[str, float] = defaultdict(float)
    for row in rows:
        substitution = row["substitution"]
        counts[substitution] += i(row, "count")
        intervals[substitution] += i(row, "intervals")

    spectrum = [
        (
            substitution,
            counts[substitution] / intervals[substitution] * 1_000_000
            if intervals[substitution]
            else 0.0,
        )
        for substitution in SUBSTITUTIONS
    ]
    spectrum.sort(key=lambda item: -item[1])
    if not any(value > 0 for _, value in spectrum):
        return None

    labels = [item[0] for item in reversed(spectrum)]
    values = [item[1] for item in reversed(spectrum)]
    y_positions = list(range(len(labels)))
    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.barh(y_positions, values, color="#8E5CF3", height=0.72)
    ax.set_yticks(y_positions, labels)
    configure_axes(
        ax,
        title=f"{label}: substitution spectrum",
        xlabel="Substitutions per million k-mer intervals",
        ylabel="Substitution",
    )
    return (
        "substitution_spectrum",
        save_figure(fig, out_dir, figure_stem(file_prefix, "substitution_spectrum"), formats),
    )


def load_inputs(explore_dir: Path) -> dict[str, list[dict[str, str]]]:
    return {
        "risk_ranked_kmers": read_csv(explore_dir / "risk_ranked_kmers.csv"),
        "rc_collapsed_kmers": read_csv(explore_dir / "rc_collapsed_kmers.csv"),
        "rc_pair_diagnostics": read_csv(explore_dir / "rc_pair_diagnostics.csv"),
        "rc_strand_mirror": read_csv(explore_dir / "rc_strand_mirror.csv"),
        "quality_calibration": read_csv(explore_dir / "quality_calibration.csv"),
        "hp_phase_profile": read_csv(explore_dir / "hp_phase_profile.csv"),
        "position_heatmap": read_csv(explore_dir / "position_heatmap.csv"),
        "motif_enrichment": read_csv(explore_dir / "motif_enrichment.csv"),
        "substitution_spectrum": read_csv(explore_dir / "substitution_spectrum.csv"),
    }


def first_supported(
    rows: list[dict[str, str]], min_intervals: int, interval_key: str = "intervals"
) -> dict[str, str] | None:
    for row in rows:
        if passes_support(row, min_intervals, interval_key):
            return row
    return rows[0] if rows else None


def write_manifest(
    path: Path,
    label: str,
    platform: str,
    figure_prefix: str,
    explore_dir: Path,
    generated: list[tuple[str, list[Path]]],
    inputs: dict[str, list[dict[str, str]]],
    min_intervals: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    captions = {
        "context_error_classes": "Reverse-complement-collapsed k-mers ranked by error enrichment, with error classes stacked per k-mer interval.",
        "top_kmer_lollipop": "High-risk raw k-mers ranked by error rate, colored by their dominant mismatch/indel behavior.",
        "top_kmer_clustered_heatmap": "Top raw k-mers clustered by sequence and error-profile features, with a dendrogram when SciPy is available.",
        "top_kmer_sequence_heatmap": "Base-by-position view of top raw k-mers with an adjacent error-rate bar.",
        "motif_aligned_contexts": "Top high-risk k-mers containing the leading enriched motif, aligned to the motif start.",
        "rc_pair_asymmetry": "Direct comparison of each canonical k-mer with its reverse complement; points away from the diagonal mark orientation-specific context effects.",
        "rc_strand_mirror": "Reverse-complement pair differences split by alignment strand, useful for separating true pair asymmetry from ordinary strand bias.",
        "quality_calibration": "Predicted mean QV against empirical k-mer QV; points below the diagonal are overconfident contexts.",
        "homopolymer_profile": "Error rate as a function of homopolymer run length and base within the k-mer.",
        "offset_base_heatmap": "Compact heatmap of error rate by offset and reference base, useful as a bridge from global context enrichment to positional/base-specific effects.",
        "motif_enrichment": "Short motifs enriched among high-error k-mers without specifying a motif model ahead of time.",
        "substitution_spectrum": "Substitution spectrum summarized across offsets from the offset-aware substitution table.",
    }
    lines = [
        f"# {platform} K-mer Error Context Figures",
        "",
        f"Dataset: `{label}`",
        "",
        f"Platform: `{platform}`",
        "",
        f"Exploration directory: `{explore_dir}`",
        "",
        f"Minimum support used for figures: `{min_intervals}` intervals",
        "",
        "## Highlights",
        "",
    ]
    top_rc = first_supported(inputs["rc_collapsed_kmers"], min_intervals)
    if top_rc:
        lines.append(
            "- Top reverse-complement-collapsed context: "
            f"`{top_rc['canonical_kmer']}` ({top_rc['members']}); "
            f"enrichment `{top_rc['enrichment']}x`, error rate `{percent(f(top_rc, 'error_rate')):.2f}%`."
        )
    top_pair = first_supported(
        inputs["rc_pair_diagnostics"], min_intervals, "min_member_intervals"
    )
    if top_pair:
        lines.append(
            "- Largest supported reverse-complement pair difference: "
            f"`{top_pair['member_a']}/{top_pair['member_b']}`; "
            f"absolute delta `{percent(f(top_pair, 'abs_delta_error_rate')):.2f}` percentage points."
        )
    top_quality = first_supported(inputs["quality_calibration"], min_intervals)
    if top_quality:
        lines.append(
            "- Most overconfident supported context: "
            f"`{top_quality['kmer']}`; predicted minus empirical QV "
            f"`{f(top_quality, 'overconfidence_qv'):.2f}`."
        )

    lines.extend(["", "## Figures", ""])
    for key, paths in generated:
        display_path = paths[0]
        lines.append(f"### `{display_path.name}`")
        lines.append("")
        lines.append(captions[key])
        lines.append("")
        lines.append(f"Files: {', '.join(f'`{path.name}`' for path in paths)}")
        lines.append("")

    lines.extend(
        [
            "## Presentation Builder Notes",
            "",
            "Recommended slide order:",
            "",
            f"1. Start with `{figure_stem(figure_prefix, 'context_error_classes')}`: it is the clearest anchor for the claim that sequencing errors are sequence-context dependent.",
            f"2. Use `{figure_stem(figure_prefix, 'top_kmer_lollipop')}` and `{figure_stem(figure_prefix, 'top_kmer_sequence_heatmap')}` to make the actual top k-mer sequences visible.",
            f"3. Use `{figure_stem(figure_prefix, 'top_kmer_clustered_heatmap')}` when you want to show that the top contexts form feature/sequence groups rather than isolated anecdotes.",
            f"4. Follow with `{figure_stem(figure_prefix, 'homopolymer_run_profile')}` to connect top contexts to platform-specific homopolymer behavior.",
            f"5. Use `{figure_stem(figure_prefix, 'motif_aligned_contexts')}` and `{figure_stem(figure_prefix, 'motif_enrichment')}` as discovery slides for compact sequence signatures.",
            f"6. Use `{figure_stem(figure_prefix, 'rc_pair_asymmetry')}` and `{figure_stem(figure_prefix, 'rc_strand_mirror')}` together: the pair scatter shows orientation asymmetry, while the mirror plot checks whether strand alone explains it.",
            f"7. Use `{figure_stem(figure_prefix, 'quality_calibration')}` to show that context effects also affect predicted-vs-empirical quality.",
            f"8. Use `{figure_stem(figure_prefix, 'substitution_spectrum')}` as supporting detail unless substitution bias is the dominant platform signal.",
            f"9. Use `{figure_stem(figure_prefix, 'offset_base_error_heatmap')}` as an optional bridge slide when explaining how BEST can move from k-mer-level summaries to offset/base-level summaries.",
            "",
            "Slide design notes:",
            "",
            "- Prefer PNGs for dense scatter plots; their SVGs intentionally rasterize point clouds to keep file sizes manageable.",
            "- Keep the dataset label in a subtitle if the deck has its own title hierarchy; the generated titles can be shortened in the slide tool.",
            "- Avoid overclaiming from one dataset alone. Phrase single-run results as a demonstration of the feature set and a pilot error-context profile.",
            "- Good first proof statement: BEST recovers platform-specific homopolymer/context signal, then adds reverse-complement and motif diagnostics that expose additional context asymmetries.",
        ]
    )

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    explore_dir = Path(args.explore_dir)
    out_dir = Path(args.out_dir) if args.out_dir else explore_dir / "figures"
    inputs = load_inputs(explore_dir)

    plotters = [
        plot_context_error_classes(
            inputs["rc_collapsed_kmers"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.top_n,
        ),
        plot_top_kmer_lollipop(
            inputs["risk_ranked_kmers"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.top_n,
        ),
        plot_top_kmer_clustered_heatmap(
            inputs["risk_ranked_kmers"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.top_n,
        ),
        plot_top_kmer_sequence_heatmap(
            inputs["risk_ranked_kmers"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.top_n,
        ),
        plot_motif_aligned_contexts(
            inputs["risk_ranked_kmers"],
            inputs["motif_enrichment"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.top_n,
        ),
        plot_rc_pair_asymmetry(
            inputs["rc_pair_diagnostics"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.annotate_top,
        ),
        plot_rc_strand_mirror(
            inputs["rc_strand_mirror"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.top_n,
        ),
        plot_quality_calibration(
            inputs["quality_calibration"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.min_intervals,
            args.annotate_top,
        ),
        plot_homopolymer_profile(
            inputs["hp_phase_profile"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.max_run_length,
        ),
        plot_offset_base_heatmap(
            inputs["position_heatmap"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
        ),
        plot_motif_enrichment(
            inputs["motif_enrichment"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
            args.top_n,
        ),
        plot_substitution_spectrum(
            inputs["substitution_spectrum"],
            out_dir,
            args.formats,
            args.label,
            args.figure_prefix,
        ),
    ]
    generated = [item for item in plotters if item is not None]
    manifest = out_dir / "figure_manifest.md"
    write_manifest(
        manifest,
        args.label,
        args.platform,
        args.figure_prefix,
        explore_dir,
        generated,
        inputs,
        args.min_intervals,
    )

    for _, paths in generated:
        for path in paths:
            print(f"wrote {path}")
    print(f"wrote {manifest}")
    if not generated:
        print("warning: no figures generated; expected kmer_explore.py CSVs were missing or empty")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
