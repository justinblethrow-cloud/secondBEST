# `best` Usage Guide

This guide will give a general overview of how to use `best`.

## Alignment Settings
The alignment `bam` file must contain the read sequences and quality scores.
The alignment CIGAR string can either use `M` or `=`/`X` for matches/mismatches.

## Example Analysis
Let's say we have aligned reads in `aln.bam` and the reference assembly
`ref.fasta.gz`. We want to collect statistics on the types of
errors that occur in the alignments. Then, we can run `best` like
```
best -t 4 aln.bam ref.fasta.gz output
```
This will use 4 threads to collect stats and generate the following files:
```
output.per_aln_stats.csv.gz
output.summary_cigar_stats.csv
output.summary_identity_stats.csv
output.summary_qual_score_stats.csv
output.summary_yield_stats.csv
```
The per-alignment stats file is gzipped to save space. The
`output.summary_cigar_stats.csv` file contains the distribution of lengths of
consecutive insertions and deletions. The `output.summary_identity_stats.csv`
file contains some general stats on the error rates across all alignments.
The `output.summary_qual_score_stats.csv` file contains the empirical Q-value
calculated from matches and mismatches for each corresponding quality score.
The `all_alignments` feature indicates the quality score stats across all
alignments. If intervals are specified, then this is computed per interval feature.
The `output.summary_yield_stats.csv` contains the yield (number of reads/bases)
above certain quality thresholds.

We can collect even more data in one run of `best`. If we run
```
best -t 4 --intervals-hp --bin-types q_len:1000 gc_content:0.05 -- aln.bam ref.fasta.gz output
```
then we will get two more files:
```
output.summary_feature_stats.csv
output.summary_bin_stats.csv
```
The `output.summary_feature_stats.csv` file contains stats stratified by
intervals. In this case, we use `--intervals-hp` to indicate that the intervals
are homopolymer regions, so this will produce the error types at homopolymer
regions of different lengths. The `output.summary_bin_stats.csv` file will
contain the error types for reads binned by both the read (query) length (`q_len`,
bin by increments of 1000bp) and the GC content (`gc_content`, bin by increments
of 0.05).

It is also possible to use bed files as custom intervals. These files can have
three columns (all intervals will have the same feature) or four columns, where
the last column indicates the feature. The feature stats are aggregated across
all bed intervals with the same feature.

To summarize error profiles by reference sequence context, use
`--intervals-kmer` with one or more k-mer lengths. For example,
```
best -t 4 --intervals-kmer 7 -- aln.bam ref.fasta.gz output
```
This writes k-mer rows to `output.summary_feature_stats.csv`, where the
`feature` column is the A/C/G/T k-mer sequence. When BEST can use the optimized
k-mer-only path, it also writes `output.summary_kmer_stats.csv`, which has the
same aggregate statistics plus explicit `kmer_len` and `kmer` columns. K-mers
are generated from the aligned reference span, k-mers containing ambiguous bases
are skipped, and reverse-strand alignments are labeled using the
reverse-complemented k-mer so features are reported in read orientation.
Because k-mer intervals overlap, k-mer-stratified base counts are context counts
and are not expected to sum to the genome-wide base counts.

Optimized k-mer runs also write two compact rollup files for interpretation:
`output.summary_kmer_length_stats.csv` aggregates all observed k-mers by k-mer
length, and `output.summary_kmer_context_stats.csv` aggregates rows by
interpretable sequence-context features. The context file currently includes GC
count, central base or central dinucleotide, maximum homopolymer run length,
edge bases, and whether the k-mer is a reverse-complement palindrome. Context
rows are separate stratifications; the same k-mer contributes to multiple
context types.

Use `--kmer-position-stats` with the optimized k-mer path to also write
`output.summary_kmer_position_stats.csv`, which breaks errors down by offset
within each reported k-mer. In this file, `offset` is zero-based within the
reported read-oriented k-mer and `kmer_base` is the base at that offset.
Position-level matches are derived from k-mer occurrences after subtracting
observed mismatches, deletions, and skipped reference positions; insertions are
attributed to the current reference offset using the same convention as the
interval feature summary. Position runs also write
`output.summary_kmer_position_profile_stats.csv`, a compact profile aggregated
by k-mer length, offset, and base. Rows with `kmer_base` set to `*` aggregate all
bases at that offset.

For k-mer runs with many features, `summary_qual_score_stats.csv` can become
large because it reports empirical quality values for each feature by default.
Use `--no-feature-qual-score-stats` to keep the `all_alignments` quality-score
summary while skipping feature-specific quality-score rows and their associated
in-memory counters:
```
best -t 4 --intervals-kmer 7 --no-feature-qual-score-stats -- aln.bam ref.fasta.gz output
```
This form also enables the optimized encoded k-mer backend for k-mer-only runs
with lengths up to 32 bases.

To include the per-offset k-mer summary:
```
best -t 4 --intervals-kmer 7 --no-feature-qual-score-stats --kmer-position-stats -- aln.bam ref.fasta.gz output
```

Use `--kmer-advanced-stats` when substitution identity or strand asymmetry is
scientifically important. This writes `output.summary_kmer_strand_stats.csv`,
`output.summary_kmer_substitution_stats.csv`, and
`output.summary_kmer_substitution_profile_stats.csv`. Substitution bases are
reported in the same read-oriented frame as the k-mer labels. The detailed
substitution table can be large for high k values, so this mode is opt-in.

BEST processes BAM records in parallel batches. The default batch size is 64
records, which reduces scheduler overhead for high-thread k-mer runs while
keeping memory use bounded. Use `--record-batch-size` to tune this when
benchmarking unusual inputs. BEST also uses BGZF decompression workers for BAM
input; these are additional threads beyond `--threads`. The default is
`--bam-reader-threads 8`, which was fastest in the included k-mer benchmark on a
high-core machine. Decrease this value, for example to 1 or 2, on low-core or
shared machines.

## K-mer Benchmarking

The repository includes `scripts/benchmark_kmer.py` to run reproducible k-mer
benchmark matrices under `/usr/bin/time -v`. The script writes raw per-run CSV
rows, a JSON summary, and a Markdown report with environment details, median
timings, peak memory, output-signature checks, and candidate proof statements.
For example:
```
python3 scripts/benchmark_kmer.py \
    --best target/release/best \
    --bam benchmarking/simplex_all.hprc_hg002.downsampled.bam \
    --reference HPRC-HG002/HPRC-HG002.fasta \
    --output-dir target/tmp/kmer_benchmark \
    --threads 16 \
    --kmers 7 \
    --record-batch-size 64,128,256,512 \
    --bam-reader-threads 1,2,4,8 \
    --repeats 3 \
    --no-position-modes \
    --no-baseline
```

Use `--checksum-inputs` when preparing final numbers for publication or shared
reports; it hashes the full BAM and reference so the report identifies the exact
input bytes. Use `--report-only` to regenerate the Markdown and JSON reports
from an existing `benchmark_results.csv` without rerunning BEST.

On the available 5.0 GiB HPRC HG002 downsampled BAM, the fastest k7 aggregate
setting in this matrix was `--record-batch-size 64 --bam-reader-threads 8`.
The current full-downsample benchmark notes are in `docs/kmer_benchmarks.md`.

## K-mer Exploration Reports

Use `scripts/kmer_explore.py` to turn one or more BEST k-mer output prefixes
into scientific summary tables and a Markdown report. The script consumes
existing CSVs; it does not reread the BAM.

```
python3 scripts/kmer_explore.py \
    --prefix output \
    --label ont_duplex_v1 \
    --out-dir target/tmp/kmer_explore \
    --min-intervals 1000 \
    --top-n 50
```

The report includes k-mer error enrichment, quality calibration by context,
reverse-complement-collapsed context summaries, direct reverse-complement pair
diagnostics, homopolymer phase profiles, per-offset heatmap-ready tables,
model-free motif enrichment, consensus-risk ranked k-mers, replicate reliability
when multiple prefixes are supplied, and strand/substitution summaries when
`--kmer-advanced-stats` outputs are present. See `docs/kmer_exploration.md` for
output descriptions and interpretation notes.

To make presentation figures from the exploration output:
```
python3 scripts/kmer_explore_figures.py \
    --explore-dir target/tmp/kmer_explore \
    --out-dir target/tmp/kmer_explore/figures \
    --label "ONT simplex chr20 5x k=7" \
    --platform Nanopore \
    --figure-prefix nanopore \
    --min-intervals 1000 \
    --top-n 20
```

This writes PNG/SVG figures plus `figure_manifest.md` with highlights,
captions, and notes for building an error-context slide deck. Use `--platform`
and `--figure-prefix` to label PacBio, Revio, or other datasets.

To compare two platforms or chemistries directly:
```
python3 scripts/kmer_compare_platforms.py \
    --left-dir target/tmp/kmer_simplex/explore \
    --right-dir target/tmp/kmer_revio/explore \
    --left-label Simplex \
    --right-label Revio \
    --out-dir target/tmp/kmer_platform_compare \
    --min-intervals 1000
```

This writes hexbin plots for every shared k-mer and for
reverse-complement-collapsed contexts, plus joined and outlier CSVs for
annotating platform-specific contexts.

## Help Message:
```
best 0.1.0
Daniel Liu, Daniel E. Cook
Bam Error Stats Tool (best): analysis of error types in aligned reads.

USAGE:
    best [OPTIONS] <INPUT> <REFERENCE> <STATS_PREFIX>

ARGS:
    <INPUT>
            Input BAM file

    <REFERENCE>
            Input reference FASTA file. Can be gzipped

    <STATS_PREFIX>
            Prefix for output files that contain statistics

OPTIONS:
    -b, --bin-types <BIN_TYPES>...
            Types of bins to use for per alignment stats.

            Each bin should be of the format <bin_type>:<step_size>.

            Supported bin types: q_len (read sequence length), subread_passes, mapq, mean_qual,
            gc_content, concordance_qv (phred scale Q-value)

    -h, --help
            Print help information

        --intervals-bed <INTERVALS_BED>...
            Use intervals from a BED file.

            The BED file should have the columns chrom, start, stop, and feature. The feature column
            is optional.

            This allows stats to be gathered separately for different types of intervals. Note that
            all intervals are on the reference, not the reads.

        --intervals-border <INTERVALS_BORDER>...
            Use fixed-width nonoverlapping window border regions as intervals.

            This is used to specify the window widths.

        --intervals-hp
            Use homopolymer regions in the reference as intervals

        --intervals-kmer <INTERVALS_KMER>...
            Use all A/C/G/T k-mers of the specified lengths as intervals.

            K-mers are taken from the aligned reference span and reported in read orientation.

        --intervals-match <INTERVALS_MATCH>...
            Use regions in the reference that match any of the specified subsequences as intervals

        --intervals-window <INTERVALS_WINDOW>...
            Use fixed-width nonoverlapping windows as intervals.

            This is used to specify the window widths.

        --intervals-window-pos <INTERVALS_WINDOW_POS>...
            Use fixed-width nonoverlapping windows with positions as intervals.

            This is used to specify the window widths.

    -n, --name-column <NAME_COLUMN>
            Add column with a specific name in CSV outputs

        --no-per-aln-stats
            Turn off outputting per alignment stats

        --no-feature-qual-score-stats
            Turn off interval feature rows in summary_qual_score_stats.csv.

            The all_alignments quality score summary is still reported.

        --kmer-position-stats
            Write summary_kmer_position_stats.csv for optimized k-mer runs.

            This adds per-offset error counts within each read-oriented k-mer.

        --kmer-advanced-stats
            Write strand and substitution spectrum summaries for optimized k-mer runs.

            This adds summary_kmer_strand_stats.csv,
            summary_kmer_substitution_stats.csv, and
            summary_kmer_substitution_profile_stats.csv.

        --record-batch-size <RECORD_BATCH_SIZE>
            Number of records to hand to each Rayon task.

            Larger batches reduce scheduler overhead; smaller batches can improve load balancing.

            [default: 64]

        --bam-reader-threads <BAM_READER_THREADS>
            Number of BGZF decompression workers for the BAM reader.

            These workers are in addition to the Rayon compute threads set by --threads.
            Decrease this on low-core machines.

            [default: 8]

    -t, --threads <THREADS>
            Number of threads. Will be automatically determined if this is set to 0

            [default: 0]

    -V, --version
            Print version information
```
