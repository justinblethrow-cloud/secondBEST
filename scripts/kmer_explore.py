#!/usr/bin/env python3
"""Explore BEST k-mer error summaries.

This is a downstream scientific reporting layer for BEST's k-mer CSV outputs.
It intentionally avoids re-reading BAMs: all analyses are derived from one or
more BEST output prefixes. Optional analyses are enabled automatically when the
corresponding summary files exist.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


ERROR_RATE_COLUMNS = [
    "mismatches_per_interval",
    "non_hp_ins_per_interval",
    "non_hp_del_per_interval",
    "hp_ins_per_interval",
    "hp_del_per_interval",
]

BASES = ("A", "C", "G", "T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prefix",
        action="append",
        required=True,
        help="BEST output prefix. Can be provided multiple times for replicate analysis.",
    )
    parser.add_argument(
        "--label",
        action="append",
        help="Label for each --prefix. Defaults to the prefix basename.",
    )
    parser.add_argument("--out-dir", required=True, help="Directory for exploration outputs")
    parser.add_argument(
        "--min-intervals",
        type=int,
        default=100,
        help="Minimum k-mer intervals for ranked k-mer and motif analyses",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of rows to include in Markdown top tables",
    )
    parser.add_argument(
        "--bad-quantile",
        type=float,
        default=0.95,
        help="Error-rate quantile used to define high-error k-mers for motif discovery",
    )
    parser.add_argument(
        "--motif-lengths",
        default="2,3,4",
        help="Comma-separated motif lengths for model-free motif discovery",
    )
    args = parser.parse_args()
    if args.label and len(args.label) != len(args.prefix):
        parser.error("--label must be provided once per --prefix")
    if args.top_n <= 0:
        parser.error("--top-n must be positive")
    if args.min_intervals < 0:
        parser.error("--min-intervals must be nonnegative")
    if not (0.0 < args.bad_quantile < 1.0):
        parser.error("--bad-quantile must be between 0 and 1")
    return args


def parse_int_list(raw: str) -> list[int]:
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    if not values:
        raise ValueError("expected at least one motif length")
    if any(value <= 0 for value in values):
        raise ValueError("motif lengths must be positive")
    return values


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def i(row: dict[str, str], key: str, default: int = 0) -> int:
    value = row.get(key, "")
    if value == "":
        return default
    return int(float(value))


def error_per_interval(row: dict[str, str]) -> float:
    return sum(f(row, column) for column in ERROR_RATE_COLUMNS)


def observed_errors(row: dict[str, str]) -> float:
    return i(row, "intervals") * error_per_interval(row)


def error_rate(row: dict[str, str]) -> float:
    return max(0.0, 1.0 - f(row, "identity", 1.0))


def load_sample(prefix: Path, label: str) -> dict[str, object]:
    def summary(name: str) -> list[dict[str, str]]:
        return read_csv(Path(f"{prefix}.{name}"))

    return {
        "label": label,
        "prefix": str(prefix),
        "kmer": summary("summary_kmer_stats.csv"),
        "position": summary("summary_kmer_position_stats.csv"),
        "position_profile": summary("summary_kmer_position_profile_stats.csv"),
        "context": summary("summary_kmer_context_stats.csv"),
        "strand": summary("summary_kmer_strand_stats.csv"),
        "substitution": summary("summary_kmer_substitution_stats.csv"),
        "substitution_profile": summary("summary_kmer_substitution_profile_stats.csv"),
    }


def kmer_enrichment(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        kmer_rows = sample["kmer"]
        by_k: dict[int, list[dict[str, str]]] = defaultdict(list)
        for row in kmer_rows:
            by_k[i(row, "kmer_len")].append(row)

        for kmer_len, group in by_k.items():
            total_intervals = sum(i(row, "intervals") for row in group)
            total_errors = sum(observed_errors(row) for row in group)
            background = total_errors / total_intervals if total_intervals else 0.0
            for row in group:
                intervals = i(row, "intervals")
                obs = observed_errors(row)
                expected = intervals * background
                enrichment = obs / expected if expected > 0 else 0.0
                z_score = (obs - expected) / math.sqrt(expected) if expected > 0 else 0.0
                rows.append(
                    {
                        "sample": sample["label"],
                        "kmer_len": kmer_len,
                        "kmer": row["kmer"],
                        "intervals": intervals,
                        "observed_errors": f"{obs:.3f}",
                        "expected_errors": f"{expected:.3f}",
                        "enrichment": f"{enrichment:.6f}",
                        "z_score": f"{z_score:.6f}",
                        "identity": row["identity"],
                        "error_rate": f"{error_rate(row):.8f}",
                    }
                )

    rows.sort(key=lambda row: (row["sample"], -float(row["z_score"]), -float(row["enrichment"])))
    return rows


def quality_calibration(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        for row in sample["kmer"]:
            empirical_qv = f(row, "identity_qv")
            mean_qual = f(row, "mean_qual")
            rows.append(
                {
                    "sample": sample["label"],
                    "kmer_len": row["kmer_len"],
                    "kmer": row["kmer"],
                    "intervals": row["intervals"],
                    "mean_qual": f"{mean_qual:.6f}",
                    "empirical_qv": f"{empirical_qv:.6f}",
                    "qv_delta_empirical_minus_predicted": f"{empirical_qv - mean_qual:.6f}",
                    "overconfidence_qv": f"{mean_qual - empirical_qv:.6f}",
                    "identity": row["identity"],
                }
            )
    rows.sort(key=lambda row: (row["sample"], -float(row["overconfidence_qv"])))
    return rows


def hp_phase_profile(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int, str, str, int], dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    for sample in samples:
        for row in sample["position"]:
            kmer = row["kmer"]
            offset = i(row, "offset")
            run_start, run_end = run_bounds(kmer, offset)
            run_len = run_end - run_start + 1
            phase = phase_label(run_start, run_end, offset, len(kmer))
            key = (sample["label"], i(row, "kmer_len"), row["kmer_base"], phase, run_len)
            dest = grouped[key]
            dest["kmers"] += 1
            dest["intervals"] += i(row, "intervals")
            for column in [
                "matches",
                "mismatches",
                "non_hp_ins",
                "non_hp_del",
                "hp_ins",
                "hp_del",
                "skips",
            ]:
                dest[column] += i(row, column)

    rows = []
    for (sample, kmer_len, base, phase, run_len), values in grouped.items():
        errors = sum(values[column] for column in [
            "mismatches",
            "non_hp_ins",
            "non_hp_del",
            "hp_ins",
            "hp_del",
        ])
        observations = values["matches"] + errors
        rows.append(
            {
                "sample": sample,
                "kmer_len": kmer_len,
                "base": base,
                "phase": phase,
                "run_length_within_kmer": run_len,
                "position_rows": int(values["kmers"]),
                "intervals": int(values["intervals"]),
                "matches": int(values["matches"]),
                "mismatches": int(values["mismatches"]),
                "non_hp_ins": int(values["non_hp_ins"]),
                "non_hp_del": int(values["non_hp_del"]),
                "hp_ins": int(values["hp_ins"]),
                "hp_del": int(values["hp_del"]),
                "error_rate": f"{errors / observations:.8f}" if observations else "0.00000000",
            }
        )
    rows.sort(key=lambda row: (row["sample"], row["kmer_len"], row["base"], row["run_length_within_kmer"], row["phase"]))
    return rows


def position_heatmap(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        for row in sample["position_profile"]:
            intervals = i(row, "intervals")
            errors = sum(i(row, column) for column in [
                "mismatches",
                "non_hp_ins",
                "non_hp_del",
                "hp_ins",
                "hp_del",
            ])
            rows.append(
                {
                    "sample": sample["label"],
                    "kmer_len": row["kmer_len"],
                    "offset": row["offset"],
                    "kmer_base": row["kmer_base"],
                    "intervals": intervals,
                    "error_rate": row["error_rate"],
                    "mismatch_rate": rate(i(row, "mismatches"), intervals),
                    "non_hp_ins_rate": rate(i(row, "non_hp_ins"), intervals),
                    "non_hp_del_rate": rate(i(row, "non_hp_del"), intervals),
                    "hp_ins_rate": rate(i(row, "hp_ins"), intervals),
                    "hp_del_rate": rate(i(row, "hp_del"), intervals),
                }
            )
    rows.sort(key=lambda row: (row["sample"], int(row["kmer_len"]), int(row["offset"]), row["kmer_base"]))
    return rows


def motif_enrichment(
    samples: list[dict[str, object]],
    motif_lengths: list[int],
    min_intervals: int,
    bad_quantile: float,
) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        eligible = [
            row for row in sample["kmer"] if i(row, "intervals") >= min_intervals
        ]
        if not eligible:
            continue
        threshold = quantile([error_rate(row) for row in eligible], bad_quantile)
        bad = {row_key(row) for row in eligible if error_rate(row) >= threshold}
        total_bad = len(bad)
        total_good = len(eligible) - total_bad
        motif_to_keys: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
        for row in eligible:
            key = row_key(row)
            kmer = row["kmer"]
            for length in motif_lengths:
                if length > len(kmer):
                    continue
                for start in range(0, len(kmer) - length + 1):
                    motif = kmer[start : start + length]
                    motif_to_keys[("anywhere", motif)].add(key)
                    motif_to_keys[(f"offset_{start}", motif)].add(key)

        for (motif_type, motif), keys in motif_to_keys.items():
            bad_with = len(keys & bad)
            total_with = len(keys)
            good_with = total_with - bad_with
            bad_without = total_bad - bad_with
            good_without = total_good - good_with
            odds_ratio = ((bad_with + 0.5) * (good_without + 0.5)) / (
                (good_with + 0.5) * (bad_without + 0.5)
            )
            rows.append(
                {
                    "sample": sample["label"],
                    "motif_type": motif_type,
                    "motif": motif,
                    "motif_length": len(motif),
                    "bad_quantile_threshold": f"{threshold:.8f}",
                    "bad_with_motif": bad_with,
                    "kmers_with_motif": total_with,
                    "bad_fraction_with_motif": f"{bad_with / total_with:.8f}",
                    "bad_fraction_background": f"{total_bad / len(eligible):.8f}",
                    "odds_ratio": f"{odds_ratio:.6f}",
                }
            )
    rows.sort(key=lambda row: (row["sample"], -float(row["odds_ratio"]), -int(row["bad_with_motif"])))
    return rows


def risk_ranked_kmers(samples: list[dict[str, object]], min_intervals: int) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        for row in sample["kmer"]:
            intervals = i(row, "intervals")
            if intervals < min_intervals:
                continue
            errors_pi = error_per_interval(row)
            indel_pi = sum(f(row, column) for column in [
                "non_hp_ins_per_interval",
                "non_hp_del_per_interval",
                "hp_ins_per_interval",
                "hp_del_per_interval",
            ])
            non_hp_indel_pi = f(row, "non_hp_ins_per_interval") + f(row, "non_hp_del_per_interval")
            mismatch_pi = f(row, "mismatches_per_interval")
            indel_fraction = indel_pi / errors_pi if errors_pi > 0 else 0.0
            non_hp_indel_fraction = non_hp_indel_pi / errors_pi if errors_pi > 0 else 0.0
            mismatch_fraction = mismatch_pi / errors_pi if errors_pi > 0 else 0.0
            support = math.log10(intervals + 1)
            risk = error_rate(row) * support * (
                1.0 + non_hp_indel_fraction + 0.5 * mismatch_fraction
            )
            rows.append(
                {
                    "sample": sample["label"],
                    "kmer_len": row["kmer_len"],
                    "kmer": row["kmer"],
                    "intervals": intervals,
                    "identity": row["identity"],
                    "error_rate": f"{error_rate(row):.8f}",
                    "errors_per_interval": f"{errors_pi:.8f}",
                    "mismatch_fraction": f"{mismatch_fraction:.8f}",
                    "indel_fraction": f"{indel_fraction:.8f}",
                    "non_hp_indel_fraction": f"{non_hp_indel_fraction:.8f}",
                    "risk_score": f"{risk:.8f}",
                }
            )
    rows.sort(key=lambda row: (row["sample"], -float(row["risk_score"])))
    return rows


def replicate_reliability(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    if len(samples) < 2:
        return []
    by_key: dict[tuple[str, str], list[tuple[str, dict[str, str]]]] = defaultdict(list)
    for sample in samples:
        for row in sample["kmer"]:
            by_key[row_key(row)].append((sample["label"], row))

    rows = []
    for (kmer_len, kmer), values in by_key.items():
        if len(values) < 2:
            continue
        rates = [error_rate(row) for _, row in values]
        intervals = [i(row, "intervals") for _, row in values]
        mean_rate = statistics.mean(rates)
        sd_rate = statistics.stdev(rates) if len(rates) > 1 else 0.0
        rows.append(
            {
                "kmer_len": kmer_len,
                "kmer": kmer,
                "samples": len(values),
                "mean_error_rate": f"{mean_rate:.8f}",
                "sd_error_rate": f"{sd_rate:.8f}",
                "cv_error_rate": f"{sd_rate / mean_rate:.8f}" if mean_rate else "0.00000000",
                "min_error_rate": f"{min(rates):.8f}",
                "max_error_rate": f"{max(rates):.8f}",
                "max_delta_error_rate": f"{max(rates) - min(rates):.8f}",
                "min_intervals": min(intervals),
                "labels": ";".join(label for label, _ in values),
            }
        )
    rows.sort(key=lambda row: (-float(row["mean_error_rate"]), float(row["cv_error_rate"])))
    return rows


def strand_asymmetry(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        grouped: dict[tuple[str, str], dict[str, dict[str, str]]] = defaultdict(dict)
        for row in sample["strand"]:
            grouped[row_key(row)][row["strand"]] = row
        for (kmer_len, kmer), strands in grouped.items():
            if "forward" not in strands or "reverse" not in strands:
                continue
            fwd = 1.0 - f(strands["forward"], "identity", 1.0)
            rev = 1.0 - f(strands["reverse"], "identity", 1.0)
            rows.append(
                {
                    "sample": sample["label"],
                    "kmer_len": kmer_len,
                    "kmer": kmer,
                    "forward_intervals": strands["forward"]["intervals"],
                    "reverse_intervals": strands["reverse"]["intervals"],
                    "forward_error_rate": f"{fwd:.8f}",
                    "reverse_error_rate": f"{rev:.8f}",
                    "delta_reverse_minus_forward": f"{rev - fwd:.8f}",
                    "log2_reverse_forward_ratio": f"{log2_ratio(rev, fwd):.8f}",
                }
            )
    rows.sort(key=lambda row: (row["sample"], -abs(float(row["delta_reverse_minus_forward"]))))
    return rows


def substitution_spectrum(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sample in samples:
        for row in sample["substitution_profile"]:
            count = i(row, "count")
            intervals = i(row, "intervals")
            rows.append(
                {
                    "sample": sample["label"],
                    "kmer_len": row["kmer_len"],
                    "offset": row["offset"],
                    "ref_base": row["ref_base"],
                    "read_base": row["read_base"],
                    "substitution": f"{row['ref_base']}>{row['read_base']}",
                    "kmers": row["kmers"],
                    "intervals": intervals,
                    "count": count,
                    "rate_per_interval": row["rate_per_interval"],
                    "count_per_million_intervals": f"{count / intervals * 1_000_000:.6f}"
                    if intervals
                    else "0.000000",
                }
            )
    rows.sort(key=lambda row: (row["sample"], -int(row["count"]), row["substitution"]))
    return rows


def run_bounds(kmer: str, offset: int) -> tuple[int, int]:
    base = kmer[offset]
    start = offset
    while start > 0 and kmer[start - 1] == base:
        start -= 1
    end = offset
    while end + 1 < len(kmer) and kmer[end + 1] == base:
        end += 1
    return start, end


def phase_label(start: int, end: int, offset: int, kmer_len: int) -> str:
    if start == end:
        return "singleton"
    if offset == start and start == 0:
        return "left_edge"
    if offset == end and end == kmer_len - 1:
        return "right_edge"
    if offset == start:
        return "left_boundary"
    if offset == end:
        return "right_boundary"
    return "interior"


def rate(count: int, intervals: int) -> str:
    return f"{count / intervals:.8f}" if intervals else "0.00000000"


def row_key(row: dict[str, str]) -> tuple[str, str]:
    return (row["kmer_len"], row["kmer"])


def quantile(values: list[float], q: float) -> float:
    values = sorted(values)
    if not values:
        return 0.0
    idx = min(len(values) - 1, max(0, math.ceil(q * len(values)) - 1))
    return values[idx]


def log2_ratio(numerator: float, denominator: float) -> float:
    return math.log2((numerator + 1e-12) / (denominator + 1e-12))


def markdown_table(rows: list[dict[str, object]], columns: list[str], limit: int) -> str:
    if not rows:
        return "No rows."
    shown = rows[:limit]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in shown:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def write_report(
    path: Path,
    output_files: dict[str, Path],
    outputs: dict[str, list[dict[str, object]]],
    args: argparse.Namespace,
) -> None:
    lines = [
        "# BEST k-mer exploration report",
        "",
        "## Inputs",
        "",
    ]
    for idx, prefix in enumerate(args.prefix):
        label = args.label[idx] if args.label else Path(prefix).name
        lines.append(f"- `{label}`: `{prefix}`")

    lines.extend(
        [
            "",
            "## Generated Files",
            "",
        ]
    )
    for name, file_path in output_files.items():
        lines.append(f"- `{name}`: `{file_path}`")

    sections = [
        (
            "Top Enriched Error Contexts",
            "kmer_enrichment",
            ["sample", "kmer_len", "kmer", "intervals", "enrichment", "z_score", "error_rate"],
        ),
        (
            "Most Overconfident Contexts",
            "quality_calibration",
            [
                "sample",
                "kmer_len",
                "kmer",
                "intervals",
                "mean_qual",
                "empirical_qv",
                "overconfidence_qv",
            ],
        ),
        (
            "Top Motif Enrichments",
            "motif_enrichment",
            [
                "sample",
                "motif_type",
                "motif",
                "bad_with_motif",
                "kmers_with_motif",
                "odds_ratio",
            ],
        ),
        (
            "Consensus-Risk Ranked K-mers",
            "risk_ranked_kmers",
            [
                "sample",
                "kmer_len",
                "kmer",
                "intervals",
                "error_rate",
                "non_hp_indel_fraction",
                "risk_score",
            ],
        ),
        (
            "Largest Strand Asymmetries",
            "strand_asymmetry",
            [
                "sample",
                "kmer_len",
                "kmer",
                "forward_error_rate",
                "reverse_error_rate",
                "delta_reverse_minus_forward",
            ],
        ),
        (
            "Substitution Spectrum",
            "substitution_spectrum",
            [
                "sample",
                "kmer_len",
                "offset",
                "substitution",
                "count",
                "count_per_million_intervals",
            ],
        ),
        (
            "Replicate Reliability",
            "replicate_reliability",
            [
                "kmer_len",
                "kmer",
                "samples",
                "mean_error_rate",
                "sd_error_rate",
                "cv_error_rate",
            ],
        ),
    ]

    for title, key, columns in sections:
        if key not in outputs or not outputs[key]:
            continue
        lines.extend(["", f"## {title}", "", markdown_table(outputs[key], columns, args.top_n)])

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    motif_lengths = parse_int_list(args.motif_lengths)
    labels = args.label or [Path(prefix).name for prefix in args.prefix]
    samples = [
        load_sample(Path(prefix), label)
        for prefix, label in zip(args.prefix, labels)
    ]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "kmer_enrichment": kmer_enrichment(samples),
        "quality_calibration": quality_calibration(samples),
        "hp_phase_profile": hp_phase_profile(samples),
        "position_heatmap": position_heatmap(samples),
        "motif_enrichment": motif_enrichment(
            samples, motif_lengths, args.min_intervals, args.bad_quantile
        ),
        "risk_ranked_kmers": risk_ranked_kmers(samples, args.min_intervals),
        "replicate_reliability": replicate_reliability(samples),
        "strand_asymmetry": strand_asymmetry(samples),
        "substitution_spectrum": substitution_spectrum(samples),
    }

    fieldnames = {
        "kmer_enrichment": [
            "sample",
            "kmer_len",
            "kmer",
            "intervals",
            "observed_errors",
            "expected_errors",
            "enrichment",
            "z_score",
            "identity",
            "error_rate",
        ],
        "quality_calibration": [
            "sample",
            "kmer_len",
            "kmer",
            "intervals",
            "mean_qual",
            "empirical_qv",
            "qv_delta_empirical_minus_predicted",
            "overconfidence_qv",
            "identity",
        ],
        "hp_phase_profile": [
            "sample",
            "kmer_len",
            "base",
            "phase",
            "run_length_within_kmer",
            "position_rows",
            "intervals",
            "matches",
            "mismatches",
            "non_hp_ins",
            "non_hp_del",
            "hp_ins",
            "hp_del",
            "error_rate",
        ],
        "position_heatmap": [
            "sample",
            "kmer_len",
            "offset",
            "kmer_base",
            "intervals",
            "error_rate",
            "mismatch_rate",
            "non_hp_ins_rate",
            "non_hp_del_rate",
            "hp_ins_rate",
            "hp_del_rate",
        ],
        "motif_enrichment": [
            "sample",
            "motif_type",
            "motif",
            "motif_length",
            "bad_quantile_threshold",
            "bad_with_motif",
            "kmers_with_motif",
            "bad_fraction_with_motif",
            "bad_fraction_background",
            "odds_ratio",
        ],
        "risk_ranked_kmers": [
            "sample",
            "kmer_len",
            "kmer",
            "intervals",
            "identity",
            "error_rate",
            "errors_per_interval",
            "mismatch_fraction",
            "indel_fraction",
            "non_hp_indel_fraction",
            "risk_score",
        ],
        "replicate_reliability": [
            "kmer_len",
            "kmer",
            "samples",
            "mean_error_rate",
            "sd_error_rate",
            "cv_error_rate",
            "min_error_rate",
            "max_error_rate",
            "max_delta_error_rate",
            "min_intervals",
            "labels",
        ],
        "strand_asymmetry": [
            "sample",
            "kmer_len",
            "kmer",
            "forward_intervals",
            "reverse_intervals",
            "forward_error_rate",
            "reverse_error_rate",
            "delta_reverse_minus_forward",
            "log2_reverse_forward_ratio",
        ],
        "substitution_spectrum": [
            "sample",
            "kmer_len",
            "offset",
            "ref_base",
            "read_base",
            "substitution",
            "kmers",
            "intervals",
            "count",
            "rate_per_interval",
            "count_per_million_intervals",
        ],
    }

    output_files = {}
    for key, rows in outputs.items():
        if not rows:
            continue
        path = out_dir / f"{key}.csv"
        write_csv(path, rows, fieldnames[key])
        output_files[key] = path

    report_path = out_dir / "kmer_exploration_report.md"
    write_report(report_path, output_files, outputs, args)
    print(f"wrote {report_path}")
    for path in output_files.values():
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
