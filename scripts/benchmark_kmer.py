#!/usr/bin/env python3
"""Run reproducible BEST k-mer benchmarks and write proof-ready reports.

Each benchmark executes BEST under `/usr/bin/time -v`, writes normal BEST
outputs under the requested output directory, and appends one machine-readable
CSV row with timing, memory, output size, row-count, and output-signature
metadata. After the matrix finishes, the script also writes a JSON summary and a
Markdown report suitable for sharing with colleagues.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shlex
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


SUMMARY_FILES = [
    "summary_yield_stats.csv",
    "summary_identity_stats.csv",
    "summary_feature_stats.csv",
    "summary_kmer_stats.csv",
    "summary_kmer_length_stats.csv",
    "summary_kmer_context_stats.csv",
    "summary_kmer_position_stats.csv",
    "summary_kmer_position_profile_stats.csv",
    "summary_cigar_stats.csv",
    "summary_bin_stats.csv",
    "summary_qual_score_stats.csv",
]


CSV_FIELDS = [
    "timestamp_utc",
    "run_id",
    "repeat",
    "mode",
    "kmer_len",
    "position_stats",
    "threads",
    "record_batch_size",
    "bam_reader_threads",
    "exit_status",
    "best_runtime_seconds",
    "elapsed_seconds",
    "user_seconds",
    "system_seconds",
    "cpu_percent",
    "max_rss_kb",
    "bam",
    "bam_bytes",
    "bam_sha256",
    "reference",
    "reference_bytes",
    "reference_sha256",
    "output_prefix",
    "output_signature",
    "output_byte_signature",
    "command",
]

for name in SUMMARY_FILES:
    stem = name.removesuffix(".csv")
    CSV_FIELDS.append(f"{stem}_bytes")
    CSV_FIELDS.append(f"{stem}_rows")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best", default="target/release/best", help="BEST binary path")
    parser.add_argument("--bam", required=True, help="Input BAM")
    parser.add_argument("--reference", required=True, help="Input reference FASTA")
    parser.add_argument("--output-dir", required=True, help="Directory for outputs and logs")
    parser.add_argument(
        "--threads",
        default="1,2,4,8,16,32,61,128",
        help="Comma-separated thread counts",
    )
    parser.add_argument("--kmers", default="5,7", help="Comma-separated k-mer lengths")
    parser.add_argument(
        "--record-batch-size",
        default="64",
        help="Comma-separated records-per-Rayon-task values",
    )
    parser.add_argument(
        "--bam-reader-threads",
        default="8",
        help="Comma-separated BGZF decompression worker counts for the BAM reader",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Number of repeated runs for every matrix cell",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Benchmark CSV path. Defaults to <output-dir>/benchmark_results.csv",
    )
    parser.add_argument(
        "--report-md",
        default=None,
        help="Markdown report path. Defaults to <output-dir>/benchmark_report.md",
    )
    parser.add_argument(
        "--summary-json",
        default=None,
        help="JSON summary path. Defaults to <output-dir>/benchmark_summary.json",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Human-readable label included in the report title",
    )
    parser.add_argument(
        "--checksum-inputs",
        action="store_true",
        help="SHA-256 checksum the full BAM and reference before running",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Do not run BEST; regenerate report and JSON from the CSV",
    )
    parser.add_argument(
        "--no-position-modes",
        action="store_true",
        help="Skip k-mer runs with --kmer-position-stats",
    )
    parser.add_argument(
        "--no-baseline",
        action="store_true",
        help="Skip the non-k-mer baseline mode",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running them",
    )
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("--repeats must be positive")
    return args


def parse_int_list(raw: str) -> list[int]:
    values = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        value = int(item)
        if value <= 0:
            raise ValueError(f"expected positive integer, got {value}")
        values.append(value)
    if not values:
        raise ValueError("expected at least one integer")
    return values


def elapsed_to_seconds(raw: str) -> float | None:
    raw = raw.strip()
    pieces = raw.split(":")
    try:
        if len(pieces) == 2:
            minutes, seconds = pieces
            return int(minutes) * 60 + float(seconds)
        if len(pieces) == 3:
            hours, minutes, seconds = pieces
            return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
        return float(raw)
    except ValueError:
        return None


def parse_time_output(stderr: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for line in stderr.splitlines():
        line = line.strip()
        if line.startswith("User time (seconds):"):
            parsed["user_seconds"] = line.rsplit(":", 1)[1].strip()
        elif line.startswith("System time (seconds):"):
            parsed["system_seconds"] = line.rsplit(":", 1)[1].strip()
        elif line.startswith("Percent of CPU this job got:"):
            parsed["cpu_percent"] = line.rsplit(":", 1)[1].strip().rstrip("%")
        elif line.startswith("Elapsed (wall clock) time"):
            elapsed = line.split("):", 1)[1].strip()
            seconds = elapsed_to_seconds(elapsed)
            if seconds is not None:
                parsed["elapsed_seconds"] = f"{seconds:.2f}"
        elif line.startswith("Maximum resident set size (kbytes):"):
            parsed["max_rss_kb"] = line.rsplit(":", 1)[1].strip()
    return parsed


def parse_best_runtime(stdout: str) -> str:
    match = re.search(r"Run time \(s\):\s+([0-9.]+)", stdout)
    return match.group(1) if match else ""


def file_rows(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def output_metrics(prefix: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    byte_signature_parts = []
    signature_parts = []
    for name in SUMMARY_FILES:
        stem = name.removesuffix(".csv")
        path = Path(f"{prefix}.{name}")
        if path.exists():
            size = path.stat().st_size
            rows = file_rows(path)
            file_hash = hash_file(path)
            normalized_hash = normalized_csv_hash(path)
            metrics[f"{stem}_bytes"] = str(size)
            metrics[f"{stem}_rows"] = str(rows)
            byte_signature_parts.append(f"{name}:{size}:{rows}:{file_hash}")
            signature_parts.append(f"{name}:{size}:{rows}:{normalized_hash}")
        else:
            metrics[f"{stem}_bytes"] = ""
            metrics[f"{stem}_rows"] = ""

    if byte_signature_parts:
        digest = hashlib.sha256()
        digest.update("\n".join(byte_signature_parts).encode())
        metrics["output_byte_signature"] = digest.hexdigest()
    else:
        metrics["output_byte_signature"] = ""

    if signature_parts:
        digest = hashlib.sha256()
        digest.update("\n".join(signature_parts).encode())
        metrics["output_signature"] = digest.hexdigest()
    else:
        metrics["output_signature"] = ""
    return metrics


def normalized_csv_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            digest.update("\t".join(normalize_csv_cell(cell) for cell in row).encode())
            digest.update(b"\n")
    return digest.hexdigest()


def normalize_csv_cell(cell: str) -> str:
    if "." not in cell and "e" not in cell.lower():
        return cell
    try:
        return f"{float(cell):.5f}"
    except ValueError:
        return cell


def benchmark_modes(
    kmers: list[int],
    include_position_modes: bool,
    include_baseline: bool,
) -> list[tuple[str, int | None, bool]]:
    modes: list[tuple[str, int | None, bool]] = []
    if include_baseline:
        modes.append(("baseline", None, False))
    for kmer in kmers:
        modes.append((f"k{kmer}", kmer, False))
    if include_position_modes:
        for kmer in kmers:
            modes.append((f"k{kmer}_position", kmer, True))
    return modes


def build_command(
    best: Path,
    threads: int,
    record_batch_size: int,
    bam_reader_threads: int,
    bam: Path,
    reference: Path,
    output_prefix: Path,
    kmer_len: int | None,
    position_stats: bool,
) -> list[str]:
    command = [
        str(best),
        "-t",
        str(threads),
        "--record-batch-size",
        str(record_batch_size),
        "--bam-reader-threads",
        str(bam_reader_threads),
        "--no-per-aln-stats",
        "--no-feature-qual-score-stats",
    ]
    if position_stats:
        command.append("--kmer-position-stats")
    if kmer_len is not None:
        command.extend(["--intervals-kmer", str(kmer_len)])
    command.extend(["--", str(bam), str(reference), str(output_prefix)])
    return command


def run_one(
    best: Path,
    bam: Path,
    reference: Path,
    bam_sha256: str,
    reference_sha256: str,
    output_dir: Path,
    run_id: str,
    repeat: int,
    mode: str,
    kmer_len: int | None,
    position_stats: bool,
    threads: int,
    record_batch_size: int,
    bam_reader_threads: int,
    dry_run: bool,
) -> dict[str, str]:
    run_stem = f"{mode}.t{threads}.b{record_batch_size}.bgzf{bam_reader_threads}.r{repeat}"
    output_prefix = output_dir / run_stem
    command = build_command(
        best,
        threads,
        record_batch_size,
        bam_reader_threads,
        bam,
        reference,
        output_prefix,
        kmer_len,
        position_stats,
    )
    timed_command = ["/usr/bin/time", "-v", *command]
    row = {
        "timestamp_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "run_id": run_id,
        "repeat": str(repeat),
        "mode": mode,
        "kmer_len": "" if kmer_len is None else str(kmer_len),
        "position_stats": str(position_stats).lower(),
        "threads": str(threads),
        "record_batch_size": str(record_batch_size),
        "bam_reader_threads": str(bam_reader_threads),
        "bam": str(bam),
        "bam_bytes": str(bam.stat().st_size),
        "bam_sha256": bam_sha256,
        "reference": str(reference),
        "reference_bytes": str(reference.stat().st_size),
        "reference_sha256": reference_sha256,
        "output_prefix": str(output_prefix),
        "command": shlex.join(command),
    }
    if dry_run:
        print(shlex.join(timed_command))
        row["exit_status"] = "dry-run"
        return row

    completed = subprocess.run(timed_command, text=True, capture_output=True)
    (output_dir / f"{run_stem}.stdout.log").write_text(completed.stdout)
    (output_dir / f"{run_stem}.stderr.log").write_text(completed.stderr)
    row["exit_status"] = str(completed.returncode)
    row["best_runtime_seconds"] = parse_best_runtime(completed.stdout)
    row.update(parse_time_output(completed.stderr))
    row.update(output_metrics(output_prefix))
    return row


def command_output(command: list[str], cwd: Path | None = None) -> str:
    try:
        completed = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return ""
    if completed.returncode != 0:
        return ""
    return completed.stdout.strip()


def load_rows(csv_path: Path) -> list[dict[str, str]]:
    if not csv_path.exists():
        return []
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def refresh_output_metrics(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    refreshed = []
    for row in rows:
        output_prefix = row.get("output_prefix")
        if output_prefix:
            metrics = output_metrics(Path(output_prefix))
            row = {**row, **metrics}
        refreshed.append(row)
    return refreshed


def safe_float(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def safe_int(value: str | None) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def fmt_float(value: float | None, digits: int = 2) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def group_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str]:
    return (
        row.get("mode", ""),
        row.get("kmer_len", ""),
        row.get("position_stats", ""),
        row.get("threads", ""),
        row.get("record_batch_size", ""),
        row.get("bam_reader_threads", ""),
    )


def mode_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("mode", ""),
        row.get("kmer_len", ""),
        row.get("position_stats", ""),
    )


def summarize_rows(rows: list[dict[str, str]]) -> dict[str, object]:
    successful = [row for row in rows if row.get("exit_status") == "0"]
    grouped: dict[tuple[str, str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in successful:
        grouped[group_key(row)].append(row)

    groups = []
    for key, members in grouped.items():
        elapsed = [v for row in members if (v := safe_float(row.get("elapsed_seconds"))) is not None]
        runtime = [
            v for row in members if (v := safe_float(row.get("best_runtime_seconds"))) is not None
        ]
        cpu = [v for row in members if (v := safe_float(row.get("cpu_percent"))) is not None]
        rss = [v for row in members if (v := safe_int(row.get("max_rss_kb"))) is not None]
        signatures = sorted(
            {row.get("output_signature", "") for row in members if row.get("output_signature")}
        )
        groups.append(
            {
                "mode": key[0],
                "kmer_len": key[1],
                "position_stats": key[2],
                "threads": key[3],
                "record_batch_size": key[4],
                "bam_reader_threads": key[5],
                "runs": len(members),
                "elapsed_median_seconds": median(elapsed),
                "elapsed_min_seconds": min(elapsed) if elapsed else None,
                "elapsed_max_seconds": max(elapsed) if elapsed else None,
                "best_runtime_median_seconds": median(runtime),
                "cpu_percent_median": median(cpu),
                "max_rss_kb_max": max(rss) if rss else None,
                "output_signature_count": len(signatures),
                "output_signatures": signatures,
            }
        )

    groups.sort(
        key=lambda row: (
            row["mode"],
            safe_int(str(row["kmer_len"])) or -1,
            row["position_stats"],
            safe_int(str(row["threads"])) or -1,
            safe_int(str(row["record_batch_size"])) or -1,
            safe_int(str(row["bam_reader_threads"])) or -1,
        )
    )

    best_by_mode = []
    single_thread_best: dict[tuple[str, str, str], dict[str, object]] = {}
    groups_by_mode: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for group in groups:
        key = (str(group["mode"]), str(group["kmer_len"]), str(group["position_stats"]))
        groups_by_mode[key].append(group)
        if group["threads"] == "1":
            current = single_thread_best.get(key)
            if current is None or elapsed_less(group, current):
                single_thread_best[key] = group

    for key, members in groups_by_mode.items():
        best = min(members, key=lambda row: elapsed_sort_value(row))
        speedup = None
        if key in single_thread_best:
            single = single_thread_best[key].get("elapsed_median_seconds")
            best_elapsed = best.get("elapsed_median_seconds")
            if isinstance(single, float) and isinstance(best_elapsed, float) and best_elapsed > 0:
                speedup = single / best_elapsed
        best_by_mode.append({**best, "speedup_vs_best_1_thread": speedup})

    best_by_mode.sort(key=lambda row: (str(row["mode"]), str(row["kmer_len"])))

    signatures_by_mode: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for row in successful:
        if row.get("output_signature"):
            signatures_by_mode[mode_key(row)].add(row["output_signature"])

    consistency = [
        {
            "mode": key[0],
            "kmer_len": key[1],
            "position_stats": key[2],
            "output_signature_count": len(signatures),
            "consistent": len(signatures) <= 1,
        }
        for key, signatures in sorted(signatures_by_mode.items())
    ]

    return {
        "rows": len(rows),
        "successful_rows": len(successful),
        "failed_rows": len(
            [row for row in rows if row.get("exit_status") not in {"0", "", "dry-run"}]
        ),
        "groups": groups,
        "best_by_mode": best_by_mode,
        "output_consistency_by_mode": consistency,
    }


def elapsed_sort_value(row: dict[str, object]) -> float:
    elapsed = row.get("elapsed_median_seconds")
    return elapsed if isinstance(elapsed, float) else float("inf")


def elapsed_less(lhs: dict[str, object], rhs: dict[str, object]) -> bool:
    return elapsed_sort_value(lhs) < elapsed_sort_value(rhs)


def metadata(
    args: argparse.Namespace,
    best: Path,
    bam: Path,
    reference: Path,
    bam_sha256: str,
    reference_sha256: str,
) -> dict[str, object]:
    cwd = Path.cwd()
    return {
        "generated_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "label": args.label or "",
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "git_commit": command_output(["git", "rev-parse", "HEAD"], cwd),
        "git_branch": command_output(["git", "branch", "--show-current"], cwd),
        "best_path": str(best),
        "best_version": command_output([str(best), "--version"]),
        "bam": str(bam),
        "bam_bytes": bam.stat().st_size if bam.exists() else None,
        "bam_sha256": bam_sha256,
        "reference": str(reference),
        "reference_bytes": reference.stat().st_size if reference.exists() else None,
        "reference_sha256": reference_sha256,
    }


def write_json_summary(path: Path, meta: dict[str, object], summary: dict[str, object]) -> None:
    path.write_text(json.dumps({"metadata": meta, "summary": summary}, indent=2) + "\n")


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def render_report(
    csv_path: Path,
    json_path: Path,
    meta: dict[str, object],
    summary: dict[str, object],
) -> str:
    title = "BEST k-mer benchmark report"
    if meta.get("label"):
        title += f": {meta['label']}"

    lines = [
        f"# {title}",
        "",
        f"Generated: `{meta['generated_utc']}`",
        "",
        "## Inputs and environment",
        "",
        markdown_table(
            ["Field", "Value"],
            [
                ["Platform", str(meta.get("platform", ""))],
                ["CPU count", str(meta.get("cpu_count", ""))],
                ["Python", str(meta.get("python", ""))],
                ["Git branch", str(meta.get("git_branch", ""))],
                ["Git commit", str(meta.get("git_commit", ""))],
                ["BEST binary", str(meta.get("best_path", ""))],
                ["BEST version", str(meta.get("best_version", ""))],
                ["BAM", str(meta.get("bam", ""))],
                ["BAM bytes", str(meta.get("bam_bytes", ""))],
                ["BAM SHA-256", str(meta.get("bam_sha256", ""))],
                ["Reference", str(meta.get("reference", ""))],
                ["Reference bytes", str(meta.get("reference_bytes", ""))],
                ["Reference SHA-256", str(meta.get("reference_sha256", ""))],
                ["Raw CSV", str(csv_path)],
                ["JSON summary", str(json_path)],
            ],
        ),
        "",
        "## Matrix summary",
        "",
        markdown_table(
            ["Rows", "Successful", "Failed"],
            [
                [
                    str(summary["rows"]),
                    str(summary["successful_rows"]),
                    str(summary["failed_rows"]),
                ]
            ],
        ),
        "",
        "## Best observed setting per mode",
        "",
    ]

    best_rows = []
    for row in summary["best_by_mode"]:
        rss_mib = None
        if isinstance(row.get("max_rss_kb_max"), int):
            rss_mib = row["max_rss_kb_max"] / 1024.0
        best_rows.append(
            [
                str(row["mode"]),
                str(row["kmer_len"]),
                str(row["position_stats"]),
                str(row["threads"]),
                str(row["record_batch_size"]),
                str(row["bam_reader_threads"]),
                fmt_float(row.get("elapsed_median_seconds")),
                fmt_float(row.get("elapsed_min_seconds")),
                fmt_float(row.get("elapsed_max_seconds")),
                fmt_float(row.get("cpu_percent_median")),
                fmt_float(rss_mib),
                fmt_float(row.get("speedup_vs_best_1_thread")),
                str(row["runs"]),
            ]
        )

    lines.append(
        markdown_table(
            [
                "Mode",
                "k",
                "Position",
                "Threads",
                "Batch",
                "BGZF",
                "Median s",
                "Min s",
                "Max s",
                "CPU %",
                "Max RSS MiB",
                "Speedup vs 1t",
                "Runs",
            ],
            best_rows,
        )
        if best_rows
        else "No successful benchmark rows were found."
    )

    consistency_rows = []
    for row in summary["output_consistency_by_mode"]:
        consistency_rows.append(
            [
                str(row["mode"]),
                str(row["kmer_len"]),
                str(row["position_stats"]),
                str(row["output_signature_count"]),
                "yes" if row["consistent"] else "no",
            ]
        )
    lines.extend(
        [
            "",
            "## Output consistency",
            "",
            "Output signatures hash all generated summary CSVs for a run after normalizing floating-point fields to 5 decimal places. For a fixed mode, signatures should match across thread, batch, BGZF, and repeat settings while preserving integer/count differences.",
            "",
            markdown_table(
                ["Mode", "k", "Position", "Unique signatures", "Consistent"],
                consistency_rows,
            )
            if consistency_rows
            else "No output signatures were collected.",
            "",
            "## Full grouped timings",
            "",
        ]
    )

    grouped_rows = []
    for row in summary["groups"]:
        rss_mib = None
        if isinstance(row.get("max_rss_kb_max"), int):
            rss_mib = row["max_rss_kb_max"] / 1024.0
        grouped_rows.append(
            [
                str(row["mode"]),
                str(row["kmer_len"]),
                str(row["position_stats"]),
                str(row["threads"]),
                str(row["record_batch_size"]),
                str(row["bam_reader_threads"]),
                fmt_float(row.get("elapsed_median_seconds")),
                fmt_float(row.get("elapsed_min_seconds")),
                fmt_float(row.get("elapsed_max_seconds")),
                fmt_float(row.get("cpu_percent_median")),
                fmt_float(rss_mib),
                str(row["output_signature_count"]),
                str(row["runs"]),
            ]
        )

    lines.append(
        markdown_table(
            [
                "Mode",
                "k",
                "Position",
                "Threads",
                "Batch",
                "BGZF",
                "Median s",
                "Min s",
                "Max s",
                "CPU %",
                "Max RSS MiB",
                "Signatures",
                "Runs",
            ],
            grouped_rows,
        )
        if grouped_rows
        else "No grouped timings were available."
    )

    proof_lines = proof_statements(meta, summary)
    if proof_lines:
        lines.extend(["", "## Candidate proof statements", ""])
        lines.extend(f"- {line}" for line in proof_lines)

    return "\n".join(lines) + "\n"


def proof_statements(meta: dict[str, object], summary: dict[str, object]) -> list[str]:
    statements = []
    bam = Path(str(meta.get("bam", ""))).name
    for row in summary["best_by_mode"]:
        if row["mode"] == "baseline":
            continue
        elapsed = row.get("elapsed_median_seconds")
        if not isinstance(elapsed, float):
            continue
        statement = (
            f"{row['mode']} completed on {bam} with median wall time "
            f"{elapsed:.2f}s at --threads {row['threads']}, "
            f"--record-batch-size {row['record_batch_size']}, and "
            f"--bam-reader-threads {row['bam_reader_threads']}."
        )
        speedup = row.get("speedup_vs_best_1_thread")
        if isinstance(speedup, float):
            statement += f" This was {speedup:.2f}x faster than the best single-thread run."
        statements.append(statement)

    inconsistent = [
        row for row in summary["output_consistency_by_mode"] if not row.get("consistent")
    ]
    if summary["output_consistency_by_mode"] and not inconsistent:
        statements.append(
            "All successful runs produced matching output signatures within each benchmark mode."
        )
    return statements


def main() -> int:
    args = parse_args()
    best = Path(args.best)
    bam = Path(args.bam)
    reference = Path(args.reference)
    output_dir = Path(args.output_dir)
    csv_path = Path(args.csv) if args.csv else output_dir / "benchmark_results.csv"
    report_path = Path(args.report_md) if args.report_md else output_dir / "benchmark_report.md"
    json_path = (
        Path(args.summary_json) if args.summary_json else output_dir / "benchmark_summary.json"
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    bam_sha256 = hash_file(bam) if args.checksum_inputs and bam.exists() else ""
    reference_sha256 = (
        hash_file(reference) if args.checksum_inputs and reference.exists() else ""
    )

    if not args.report_only:
        threads = parse_int_list(args.threads)
        kmers = parse_int_list(args.kmers)
        record_batch_sizes = parse_int_list(args.record_batch_size)
        bam_reader_thread_counts = parse_int_list(args.bam_reader_threads)
        modes = benchmark_modes(kmers, not args.no_position_modes, not args.no_baseline)
        run_id = dt.datetime.now(dt.UTC).strftime("%Y%m%dT%H%M%SZ")

        write_header = not csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
            if write_header:
                writer.writeheader()
            for repeat in range(1, args.repeats + 1):
                for mode, kmer_len, position_stats in modes:
                    for thread_count in threads:
                        for record_batch_size in record_batch_sizes:
                            for bam_reader_threads in bam_reader_thread_counts:
                                print(
                                    f"running repeat={repeat} mode={mode} "
                                    f"threads={thread_count} batch={record_batch_size} "
                                    f"bgzf={bam_reader_threads}",
                                    file=sys.stderr,
                                    flush=True,
                                )
                                row = run_one(
                                    best,
                                    bam,
                                    reference,
                                    bam_sha256,
                                    reference_sha256,
                                    output_dir,
                                    run_id,
                                    repeat,
                                    mode,
                                    kmer_len,
                                    position_stats,
                                    thread_count,
                                    record_batch_size,
                                    bam_reader_threads,
                                    args.dry_run,
                                )
                                writer.writerow(row)
                                handle.flush()
                                if row.get("exit_status") not in {"0", "dry-run"}:
                                    print(
                                        f"warning: repeat={repeat} mode={mode} "
                                        f"threads={thread_count} batch={record_batch_size} "
                                        f"bgzf={bam_reader_threads} exited with "
                                        f"{row['exit_status']}",
                                        file=sys.stderr,
                                    )

    rows = load_rows(csv_path)
    if args.report_only:
        rows = refresh_output_metrics(rows)
    summary = summarize_rows(rows)
    meta = metadata(args, best, bam, reference, bam_sha256, reference_sha256)
    write_json_summary(json_path, meta, summary)
    report_path.write_text(render_report(csv_path, json_path, meta, summary))
    print(f"{'read' if args.report_only else 'wrote'} {csv_path}", file=sys.stderr)
    print(f"wrote {json_path}", file=sys.stderr)
    print(f"wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
