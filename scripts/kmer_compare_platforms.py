#!/usr/bin/env python3
"""Compare k-mer error profiles between two BEST exploration outputs."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ERROR_CLASS_COLUMNS = [
    ("mismatches_per_interval", "Mismatch"),
    ("non_hp_ins_per_interval", "Non-HP insertion"),
    ("non_hp_del_per_interval", "Non-HP deletion"),
    ("hp_ins_per_interval", "HP insertion"),
    ("hp_del_per_interval", "HP deletion"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--left-dir",
        required=True,
        help="First kmer_explore.py output directory",
    )
    parser.add_argument(
        "--right-dir",
        required=True,
        help="Second kmer_explore.py output directory",
    )
    parser.add_argument(
        "--left-label",
        default="Left",
        help="Label for the first platform",
    )
    parser.add_argument(
        "--right-label",
        default="Right",
        help="Label for the second platform",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Directory for comparison figures and tables",
    )
    parser.add_argument(
        "--min-intervals",
        type=int,
        default=1000,
        help="Minimum intervals in both platforms for a k-mer to be included",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=40,
        help="Number of platform-different k-mers to write to the outlier table",
    )
    parser.add_argument(
        "--formats",
        default="png,svg",
        help="Comma-separated output formats supported by Matplotlib",
    )
    parser.add_argument(
        "--gridsize",
        type=int,
        default=55,
        help="Hexbin grid size for density plots",
    )
    args = parser.parse_args()
    if args.min_intervals < 0:
        parser.error("--min-intervals must be nonnegative")
    if args.top_n <= 0:
        parser.error("--top-n must be positive")
    if args.gridsize <= 0:
        parser.error("--gridsize must be positive")
    args.formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    if not args.formats:
        parser.error("--formats must include at least one format")
    return args


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


def error_class_sum(row: dict[str, str]) -> float:
    return sum(f(row, key) for key, _ in ERROR_CLASS_COLUMNS)


def optional_error_class_sum(row: dict[str, str]) -> float | None:
    if not any(row.get(key, "") != "" for key, _ in ERROR_CLASS_COLUMNS):
        return None
    return error_class_sum(row)


def load_joined(
    left_dir: Path,
    right_dir: Path,
    csv_name: str,
    key_field: str,
    min_intervals: int,
) -> list[dict[str, object]]:
    left_rows = {
        row[key_field]: row
        for row in read_csv(left_dir / csv_name)
        if i(row, "intervals") >= min_intervals
    }
    right_rows = {
        row[key_field]: row
        for row in read_csv(right_dir / csv_name)
        if i(row, "intervals") >= min_intervals
    }
    joined = []
    for key in sorted(left_rows.keys() & right_rows.keys()):
        left = left_rows[key]
        right = right_rows[key]
        left_error = f(left, "error_rate")
        right_error = f(right, "error_rate")
        joined.append(
            {
                "key": key,
                "left": left,
                "right": right,
                "left_error_rate": left_error,
                "right_error_rate": right_error,
                "left_error_events_per_interval": optional_error_class_sum(left),
                "right_error_events_per_interval": optional_error_class_sum(right),
                "delta_right_minus_left": right_error - left_error,
                "log2_right_left_error_ratio": log2_ratio(right_error, left_error),
            }
        )
    return joined


def log2_ratio(numerator: float, denominator: float) -> float:
    return math.log2((numerator + 1e-12) / (denominator + 1e-12))


def save_figure(fig, out_dir: Path, stem: str, formats: list[str]) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt}"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def plot_hexbin(
    joined: list[dict[str, object]],
    out_dir: Path,
    formats: list[str],
    stem: str,
    title: str,
    left_label: str,
    right_label: str,
    gridsize: int,
    value_field: str,
    axis_label: str,
) -> list[Path]:
    x_values = [100.0 * float(row[f"left_{value_field}"]) for row in joined]
    y_values = [100.0 * float(row[f"right_{value_field}"]) for row in joined]
    axis_max = max(x_values + y_values) * 1.04 if joined else 1.0

    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    hb = None
    if joined:
        hb = ax.hexbin(
            x_values,
            y_values,
            gridsize=gridsize,
            bins="log",
            mincnt=1,
            cmap="viridis",
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No k-mers passed the support filter",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="#344054",
        )
    ax.plot([0, axis_max], [0, axis_max], linestyle="--", linewidth=1, color="#475467")
    ax.set_xlim(0, axis_max)
    ax.set_ylim(0, axis_max)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.set_xlabel(f"{left_label} {axis_label} (%)")
    ax.set_ylabel(f"{right_label} {axis_label} (%)")
    ax.grid(True, color="#D8DEE9", linewidth=0.6, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if hb is not None:
        cbar = fig.colorbar(hb, ax=ax, shrink=0.84)
        cbar.set_label("K-mers per hexbin, log scale")
    return save_figure(fig, out_dir, stem, formats)


def write_joined_table(
    path: Path,
    joined: list[dict[str, object]],
    left_label: str,
    right_label: str,
    key_column: str,
) -> None:
    fieldnames = [
        key_column,
        f"{left_label}_intervals",
        f"{right_label}_intervals",
        f"{left_label}_error_rate",
        f"{right_label}_error_rate",
        "delta_right_minus_left",
        "log2_right_left_error_ratio",
        f"{left_label}_error_events_per_interval",
        f"{right_label}_error_events_per_interval",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in joined:
            left = row["left"]
            right = row["right"]
            writer.writerow(
                {
                    key_column: row["key"],
                    f"{left_label}_intervals": i(left, "intervals"),
                    f"{right_label}_intervals": i(right, "intervals"),
                    f"{left_label}_error_rate": f"{float(row['left_error_rate']):.8f}",
                    f"{right_label}_error_rate": f"{float(row['right_error_rate']):.8f}",
                    "delta_right_minus_left": f"{float(row['delta_right_minus_left']):.8f}",
                    "log2_right_left_error_ratio": f"{float(row['log2_right_left_error_ratio']):.8f}",
                    f"{left_label}_error_events_per_interval": format_optional_float(
                        row["left_error_events_per_interval"]
                    ),
                    f"{right_label}_error_events_per_interval": format_optional_float(
                        row["right_error_events_per_interval"]
                    ),
                }
            )


def format_optional_float(value: object) -> str:
    if value is None:
        return ""
    return f"{float(value):.8f}"


def write_outliers(
    path: Path,
    joined: list[dict[str, object]],
    top_n: int,
    left_label: str,
    right_label: str,
    key_column: str,
) -> None:
    ranked = sorted(
        joined,
        key=lambda row: abs(float(row["log2_right_left_error_ratio"])),
        reverse=True,
    )[:top_n]
    write_joined_table(path, ranked, left_label, right_label, key_column)


def write_manifest(
    path: Path,
    raw_joined: list[dict[str, object]],
    rc_joined: list[dict[str, object]],
    raw_paths: list[Path],
    rc_paths: list[Path],
    table_paths: list[Path],
    left_label: str,
    right_label: str,
    min_intervals: int,
) -> None:
    left_higher = sum(
        1 for row in raw_joined if float(row["left_error_rate"]) > float(row["right_error_rate"])
    )
    right_higher = sum(
        1 for row in raw_joined if float(row["right_error_rate"]) > float(row["left_error_rate"])
    )
    median_ratio = median([float(row["log2_right_left_error_ratio"]) for row in raw_joined])
    lines = [
        "# K-mer Platform Comparison",
        "",
        f"Left platform: `{left_label}`",
        "",
        f"Right platform: `{right_label}`",
        "",
        f"Minimum support: `{min_intervals}` intervals in both platforms",
        "",
        "## Summary",
        "",
        f"- Joined raw k-mers: `{len(raw_joined)}`",
        f"- Joined reverse-complement-collapsed contexts: `{len(rc_joined)}`",
        f"- Raw k-mers with higher `{left_label}` error rate: `{left_higher}`",
        f"- Raw k-mers with higher `{right_label}` error rate: `{right_higher}`",
        f"- Median log2 `{right_label}`/`{left_label}` error-rate ratio: `{median_ratio:.3f}`",
        "",
        "## Figures",
        "",
        f"- Raw k-mer error-rate hexbin: {', '.join(f'`{path.name}`' for path in raw_paths)}",
        f"- Reverse-complement-collapsed error-rate hexbin: {', '.join(f'`{path.name}`' for path in rc_paths)}",
        "",
        "## Tables",
        "",
        *[f"- `{path.name}`" for path in table_paths],
        "",
        "## Presentation Builder Notes",
        "",
        f"- Use the raw k-mer hexbin as the main {left_label}-vs-{right_label} comparison slide.",
        "- The diagonal is equal error rate; points or density above it have higher error in the right platform, below it higher error in the left platform.",
        "- Use the reverse-complement-collapsed hexbin to show whether the platform contrast survives orientation collapsing.",
        "- Pair the hexbin slide with the outlier CSV to annotate a few concrete k-mers only after the density story is clear.",
    ]
    path.write_text("\n".join(lines) + "\n")


def median(values: list[float]) -> float:
    if not values:
        return 0.0
    values = sorted(values)
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2.0


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    left_dir = Path(args.left_dir)
    right_dir = Path(args.right_dir)
    raw_joined = load_joined(
        left_dir, right_dir, "kmer_enrichment.csv", "kmer", args.min_intervals
    )
    rc_joined = load_joined(
        left_dir,
        right_dir,
        "rc_collapsed_kmers.csv",
        "canonical_kmer",
        args.min_intervals,
    )

    raw_paths = plot_hexbin(
        raw_joined,
        out_dir,
        args.formats,
        "platform_kmer_error_hexbin",
        f"{args.left_label} vs {args.right_label}: every k-mer",
        args.left_label,
        args.right_label,
        args.gridsize,
        "error_rate",
        "error rate",
    )
    rc_paths = plot_hexbin(
        rc_joined,
        out_dir,
        args.formats,
        "platform_rc_collapsed_error_hexbin",
        f"{args.left_label} vs {args.right_label}: RC-collapsed contexts",
        args.left_label,
        args.right_label,
        args.gridsize,
        "error_rate",
        "RC-collapsed error rate",
    )

    joined_path = out_dir / "platform_kmer_joined.csv"
    outlier_path = out_dir / "platform_kmer_outliers.csv"
    rc_joined_path = out_dir / "platform_rc_collapsed_joined.csv"
    rc_outlier_path = out_dir / "platform_rc_collapsed_outliers.csv"
    table_paths = [joined_path, outlier_path, rc_joined_path, rc_outlier_path]
    write_joined_table(joined_path, raw_joined, args.left_label, args.right_label, "kmer")
    write_outliers(outlier_path, raw_joined, args.top_n, args.left_label, args.right_label, "kmer")
    write_joined_table(
        rc_joined_path,
        rc_joined,
        args.left_label,
        args.right_label,
        "canonical_kmer",
    )
    write_outliers(
        rc_outlier_path,
        rc_joined,
        args.top_n,
        args.left_label,
        args.right_label,
        "canonical_kmer",
    )
    manifest_path = out_dir / "platform_comparison_manifest.md"
    write_manifest(
        manifest_path,
        raw_joined,
        rc_joined,
        raw_paths,
        rc_paths,
        table_paths,
        args.left_label,
        args.right_label,
        args.min_intervals,
    )

    for path in raw_paths + rc_paths + table_paths + [manifest_path]:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
