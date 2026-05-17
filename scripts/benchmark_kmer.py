#!/usr/bin/env python3
"""Run reproducible BEST k-mer benchmarks.

The benchmark matrix is intentionally simple: each run executes BEST under
`/usr/bin/time -v`, writes normal BEST outputs under the requested output
directory, and appends one machine-readable CSV row with timing, memory, output
size, and row-count metadata.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import subprocess
import sys
from pathlib import Path


SUMMARY_FILES = [
    "summary_yield_stats.csv",
    "summary_identity_stats.csv",
    "summary_feature_stats.csv",
    "summary_kmer_stats.csv",
    "summary_kmer_position_stats.csv",
    "summary_cigar_stats.csv",
    "summary_bin_stats.csv",
    "summary_qual_score_stats.csv",
]


CSV_FIELDS = [
    "timestamp_utc",
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
    "reference",
    "reference_bytes",
    "output_prefix",
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
        "--csv",
        default=None,
        help="Benchmark CSV path. Defaults to <output-dir>/benchmark_results.csv",
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
    return parser.parse_args()


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


def output_metrics(prefix: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    for name in SUMMARY_FILES:
        stem = name.removesuffix(".csv")
        path = Path(f"{prefix}.{name}")
        if path.exists():
            metrics[f"{stem}_bytes"] = str(path.stat().st_size)
            metrics[f"{stem}_rows"] = str(file_rows(path))
        else:
            metrics[f"{stem}_bytes"] = ""
            metrics[f"{stem}_rows"] = ""
    return metrics


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
    output_dir: Path,
    mode: str,
    kmer_len: int | None,
    position_stats: bool,
    threads: int,
    record_batch_size: int,
    bam_reader_threads: int,
    dry_run: bool,
) -> dict[str, str]:
    output_prefix = output_dir / (
        f"{mode}.t{threads}.b{record_batch_size}.bgzf{bam_reader_threads}"
    )
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
        "mode": mode,
        "kmer_len": "" if kmer_len is None else str(kmer_len),
        "position_stats": str(position_stats).lower(),
        "threads": str(threads),
        "record_batch_size": str(record_batch_size),
        "bam_reader_threads": str(bam_reader_threads),
        "bam": str(bam),
        "bam_bytes": str(bam.stat().st_size),
        "reference": str(reference),
        "reference_bytes": str(reference.stat().st_size),
        "output_prefix": str(output_prefix),
        "command": " ".join(command),
    }
    if dry_run:
        print(" ".join(timed_command))
        row["exit_status"] = "dry-run"
        return row

    log_stem = f"{mode}.t{threads}.b{record_batch_size}.bgzf{bam_reader_threads}"
    completed = subprocess.run(timed_command, text=True, capture_output=True)
    (output_dir / f"{log_stem}.stdout.log").write_text(completed.stdout)
    (output_dir / f"{log_stem}.stderr.log").write_text(completed.stderr)
    row["exit_status"] = str(completed.returncode)
    row["best_runtime_seconds"] = parse_best_runtime(completed.stdout)
    row.update(parse_time_output(completed.stderr))
    row.update(output_metrics(output_prefix))
    return row


def main() -> int:
    args = parse_args()
    best = Path(args.best)
    bam = Path(args.bam)
    reference = Path(args.reference)
    output_dir = Path(args.output_dir)
    csv_path = Path(args.csv) if args.csv else output_dir / "benchmark_results.csv"
    threads = parse_int_list(args.threads)
    kmers = parse_int_list(args.kmers)
    record_batch_sizes = parse_int_list(args.record_batch_size)
    bam_reader_thread_counts = parse_int_list(args.bam_reader_threads)
    modes = benchmark_modes(kmers, not args.no_position_modes, not args.no_baseline)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_header = not csv_path.exists()
    with csv_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for mode, kmer_len, position_stats in modes:
            for thread_count in threads:
                for record_batch_size in record_batch_sizes:
                    for bam_reader_threads in bam_reader_thread_counts:
                        print(
                            f"running mode={mode} threads={thread_count} "
                            f"batch={record_batch_size} bgzf={bam_reader_threads}",
                            file=sys.stderr,
                            flush=True,
                        )
                        row = run_one(
                            best,
                            bam,
                            reference,
                            output_dir,
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
                                f"warning: mode={mode} threads={thread_count} "
                                f"batch={record_batch_size} bgzf={bam_reader_threads} "
                                f"exited with {row['exit_status']}",
                                file=sys.stderr,
                            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
