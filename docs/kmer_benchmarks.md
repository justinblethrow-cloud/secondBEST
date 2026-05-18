# K-mer Benchmark Notes

These notes summarize local proof runs for the optimized k-mer summary path.
They are intended to make performance claims auditable without committing the
large generated benchmark outputs.

## Environment

- Date: 2026-05-18 UTC
- Commit: `c9e230a4e9a84fe8a4060f148feb095547ec3c86`
- Host platform: Linux 6.8.0-107-generic x86_64
- Reported CPU count: 256
- BEST binary: `target/release/best`
- BEST version: `best 0.1.0`
- BAM: `simplex_all.hprc_hg002.downsampled.bam`
- BAM size: 5,301,770,242 bytes
- BAM SHA-256: `6f596cd4a04066774f99bc04680273d9fe01b89c6b20f10af395e1aa96ad2807`
- Reference: `HPRC-HG002.fasta`
- Reference size: 6,096,357,567 bytes
- Reference SHA-256: `bd29273d71833636a7f3c6efc252176fb84e6a907bcd2cdcdefd24ef5dbd8211`

The benchmark harness records both raw per-run CSV rows and generated Markdown
and JSON reports. Output consistency uses a normalized signature that rounds
floating-point fields to 5 decimal places while preserving integer/count
differences. This avoids treating harmless last-digit differences from parallel
floating-point reduction order as biological or counting differences.

## Scaling Benchmark

Command:

```bash
python3 -B scripts/benchmark_kmer.py \
    --best target/release/best \
    --bam /mnt/datavault/Agentic/BEST/benchmarking/simplex_all.hprc_hg002.downsampled.bam \
    --reference /mnt/datavault/Agentic/BEST/HPRC-HG002/HPRC-HG002.fasta \
    --output-dir target/tmp/kmer_proof_scaling_20260518 \
    --threads 1,8,16,32,61,128 \
    --kmers 7 \
    --record-batch-size 64 \
    --bam-reader-threads 8 \
    --repeats 2 \
    --no-position-modes \
    --no-baseline \
    --checksum-inputs \
    --label full-downsample-k7-scaling
```

Results for k7 aggregate mode:

| Threads | Batch | BGZF | Median s | Min s | Max s | CPU % | Max RSS MiB |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 64 | 8 | 231.84 | 231.08 | 232.60 | 150.00 | 5757.82 |
| 8 | 64 | 8 | 38.02 | 37.94 | 38.10 | 925.50 | 6187.25 |
| 16 | 64 | 8 | 26.85 | 26.64 | 27.06 | 1429.50 | 5947.30 |
| 32 | 64 | 8 | 27.38 | 25.23 | 29.52 | 2481.00 | 5918.05 |
| 61 | 64 | 8 | 29.92 | 29.81 | 30.03 | 4571.50 | 6007.66 |
| 128 | 64 | 8 | 29.59 | 29.30 | 29.87 | 8997.50 | 6272.71 |

Best observed setting: `--threads 16 --record-batch-size 64
--bam-reader-threads 8`, with median wall time 26.85 seconds. This was 8.63x
faster than the best single-thread run in the same matrix. All successful runs
produced matching normalized output signatures within the k7 aggregate mode.

## Mode-Cost Benchmark

Command:

```bash
python3 -B scripts/benchmark_kmer.py \
    --best target/release/best \
    --bam /mnt/datavault/Agentic/BEST/benchmarking/simplex_all.hprc_hg002.downsampled.bam \
    --reference /mnt/datavault/Agentic/BEST/HPRC-HG002/HPRC-HG002.fasta \
    --output-dir target/tmp/kmer_proof_modes_20260518 \
    --threads 16 \
    --kmers 5,7 \
    --record-batch-size 64 \
    --bam-reader-threads 8 \
    --repeats 2 \
    --no-baseline \
    --checksum-inputs \
    --label full-downsample-kmer-mode-cost
```

Results at the tuned 16-thread setting:

| Mode | k | Position stats | Median s | Min s | Max s | CPU % | Max RSS MiB |
| --- | --- | --- | --- | --- | --- | --- | --- |
| k5 | 5 | false | 24.62 | 23.19 | 26.05 | 1522.50 | 5860.21 |
| k5_position | 5 | true | 25.16 | 23.77 | 26.56 | 1505.00 | 5862.72 |
| k7 | 7 | false | 26.55 | 25.95 | 27.16 | 1411.50 | 5973.13 |
| k7_position | 7 | true | 38.00 | 33.80 | 42.20 | 1123.50 | 6812.93 |

All successful runs produced matching normalized output signatures within each
mode. In this matrix, k5 position summaries were nearly free relative to k5
aggregate summaries, while k7 position summaries added measurable cost and
memory because they maintain a larger per-k-mer offset profile.

## Practical Recommendation

For high-core machines similar to this benchmark host, use:

```bash
best -t 16 \
    --record-batch-size 64 \
    --bam-reader-threads 8 \
    --no-feature-qual-score-stats \
    --intervals-kmer 7 \
    -- input.bam reference.fasta output
```

The 32-, 61-, and 128-thread runs did not improve median wall time for this
input, and they consumed substantially more total CPU. On smaller or shared
machines, reduce `--bam-reader-threads` to 1 or 2.
