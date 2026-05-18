// Copyright (c) 2022 Google LLC
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
// the Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use clap::Parser;

use rayon::prelude::*;

use noodles::{bam, bgzf, fasta, sam};

use fxhash::FxHashMap;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::num::NonZeroUsize;
use std::str::FromStr;
use std::sync::Mutex;
use std::time::Instant;

mod stats;
use stats::*;
mod summary;
use summary::*;
mod bed;
use bed::*;
mod intervals;
use intervals::*;
mod kmer;
use kmer::*;

const PER_ALN_STATS_NAME: &str = "per_aln_stats.csv.gz";
const YIELD_STATS_NAME: &str = "summary_yield_stats.csv";
const IDENTITY_STATS_NAME: &str = "summary_identity_stats.csv";
const FEATURE_STATS_NAME: &str = "summary_feature_stats.csv";
const KMER_STATS_NAME: &str = "summary_kmer_stats.csv";
const KMER_LENGTH_STATS_NAME: &str = "summary_kmer_length_stats.csv";
const KMER_CONTEXT_STATS_NAME: &str = "summary_kmer_context_stats.csv";
const KMER_STRAND_STATS_NAME: &str = "summary_kmer_strand_stats.csv";
const KMER_SUBSTITUTION_STATS_NAME: &str = "summary_kmer_substitution_stats.csv";
const KMER_SUBSTITUTION_PROFILE_STATS_NAME: &str = "summary_kmer_substitution_profile_stats.csv";
const KMER_POSITION_STATS_NAME: &str = "summary_kmer_position_stats.csv";
const KMER_POSITION_PROFILE_STATS_NAME: &str = "summary_kmer_position_profile_stats.csv";
const CIGAR_STATS_NAME: &str = "summary_cigar_stats.csv";
const BIN_STATS_NAME: &str = "summary_bin_stats.csv";
const QUAL_SCORE_STATS_NAME: &str = "summary_qual_score_stats.csv";
const DEFAULT_RECORD_BATCH_SIZE: usize = 64;
const DEFAULT_BAM_READER_THREADS: usize = 8;

fn run(
    input_path: String,
    reference_path: String,
    stats_prefix: String,
    bin_types: Option<Vec<BinType>>,
    intervals_types: Vec<IntervalsType>,
    name_column: Option<String>,
    output_per_aln_stats: bool,
    feature_qual_score_stats: bool,
    kmer_position_stats: bool,
    kmer_advanced_stats: bool,
    record_batch_size: usize,
    bam_reader_threads: usize,
) {
    // read reference sequences from fasta file
    let mut ref_reader = {
        let f = File::open(&reference_path).unwrap();
        let r: Box<dyn Read> = if reference_path.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        };
        fasta::Reader::new(BufReader::new(r))
    };
    let reference_seqs: FxHashMap<String, fasta::Record> = ref_reader
        .records()
        .map(|r| r.unwrap())
        .map(|r| (r.name().to_string(), r))
        .collect();

    // read bam file
    let bgzf_reader = bgzf::reader::Builder::default()
        .set_worker_count(NonZeroUsize::new(bam_reader_threads).unwrap())
        .build_from_reader(File::open(input_path).unwrap());
    let mut reader = bam::Reader::from(bgzf_reader);
    reader.read_header().unwrap();
    let references = reader.read_reference_sequences().unwrap();

    // create per alignment stats writer that is shared between threads
    let aln_stats_path = format!("{}.{}", stats_prefix, PER_ALN_STATS_NAME);
    let aln_stats_writer = if output_per_aln_stats {
        let mut w = GzEncoder::new(
            BufWriter::new(File::create(&aln_stats_path).unwrap()),
            flate2::Compression::default(),
        );
        write!(
            w,
            "{}{}\n",
            if name_column.is_some() { "name," } else { "" },
            AlnStats::header()
        )
        .unwrap();
        Some(Mutex::new(w))
    } else {
        None
    };

    let optimized_kmer_lens = optimized_kmer_lens(&intervals_types, feature_qual_score_stats);
    let use_optimized_kmer = optimized_kmer_lens.is_some();
    let optimized_kmer_lens = optimized_kmer_lens.unwrap_or_else(Vec::new);
    let feature_summary_enabled = !intervals_types.is_empty() && !use_optimized_kmer;
    let bin_types = bin_types.unwrap_or_else(Vec::new);
    let record_processor = RecordProcessor {
        references: &references,
        reference_seqs: &reference_seqs,
        intervals_types: &intervals_types,
        feature_qual_score_stats,
        use_optimized_kmer,
        name_column: name_column.as_deref(),
    };

    // Feed Rayon record batches rather than individual records to reduce scheduling overhead.
    let mut lazy_records = reader.lazy_records();
    let record_batches = std::iter::from_fn(move || {
        let mut batch = Vec::with_capacity(record_batch_size);
        for _ in 0..record_batch_size {
            match lazy_records.next() {
                Some(r) => batch.push(r.unwrap()),
                None => break,
            }
        }

        if batch.is_empty() {
            None
        } else {
            Some(batch)
        }
    });

    let mut summaries = record_batches
        .par_bridge()
        .fold(
            || {
                SummaryAccumulator::new(
                    name_column.clone(),
                    feature_summary_enabled,
                    &bin_types,
                    feature_qual_score_stats,
                    kmer_position_stats,
                    kmer_advanced_stats,
                    if use_optimized_kmer {
                        Some(&optimized_kmer_lens)
                    } else {
                        None
                    },
                )
            },
            |mut summaries, records| {
                let mut aln_stats_rows = if aln_stats_writer.is_some() {
                    Some(String::new())
                } else {
                    None
                };

                for record in records {
                    record_processor.update_summaries(
                        &mut summaries,
                        &record,
                        aln_stats_rows.as_mut(),
                    );
                }

                if let (Some(ref w), Some(rows)) = (&aln_stats_writer, aln_stats_rows) {
                    let mut w = w.lock().unwrap();
                    w.write_all(rows.as_bytes()).unwrap();
                }

                summaries
            },
        )
        .reduce(
            || {
                SummaryAccumulator::new(
                    name_column.clone(),
                    feature_summary_enabled,
                    &bin_types,
                    feature_qual_score_stats,
                    kmer_position_stats,
                    kmer_advanced_stats,
                    if use_optimized_kmer {
                        Some(&optimized_kmer_lens)
                    } else {
                        None
                    },
                )
            },
            |mut a, b| {
                a.assign_add(b);
                a
            },
        );

    write_summary(summaries.yield_summary, &stats_prefix, YIELD_STATS_NAME);

    summaries.identity_summary.total_alns = summaries.total_alns;
    write_summary(
        summaries.identity_summary,
        &stats_prefix,
        IDENTITY_STATS_NAME,
    );

    if let Some(f) = summaries.feature_summary {
        write_summary(f, &stats_prefix, FEATURE_STATS_NAME);
    }
    if let Some(ref k) = summaries.kmer_summary {
        write_summary(k.feature_summary(), &stats_prefix, FEATURE_STATS_NAME);
        write_summary(k.kmer_summary(), &stats_prefix, KMER_STATS_NAME);
        write_summary(k.length_summary(), &stats_prefix, KMER_LENGTH_STATS_NAME);
        write_summary(k.context_summary(), &stats_prefix, KMER_CONTEXT_STATS_NAME);
        if kmer_advanced_stats {
            write_summary(k.strand_summary(), &stats_prefix, KMER_STRAND_STATS_NAME);
            write_summary(
                k.substitution_summary(),
                &stats_prefix,
                KMER_SUBSTITUTION_STATS_NAME,
            );
            write_summary(
                k.substitution_profile_summary(),
                &stats_prefix,
                KMER_SUBSTITUTION_PROFILE_STATS_NAME,
            );
        }
        if kmer_position_stats {
            write_summary(
                k.position_summary(),
                &stats_prefix,
                KMER_POSITION_STATS_NAME,
            );
            write_summary(
                k.position_profile_summary(),
                &stats_prefix,
                KMER_POSITION_PROFILE_STATS_NAME,
            );
        }
    }

    write_summary(summaries.cigar_summary, &stats_prefix, CIGAR_STATS_NAME);

    if let Some(b) = summaries.bin_summary {
        write_summary(b, &stats_prefix, BIN_STATS_NAME);
    }

    write_summary(
        summaries.qual_score_summary,
        &stats_prefix,
        QUAL_SCORE_STATS_NAME,
    );
}

struct RecordProcessor<'a> {
    references: &'a sam::header::ReferenceSequences,
    reference_seqs: &'a FxHashMap<String, fasta::Record>,
    intervals_types: &'a [IntervalsType],
    feature_qual_score_stats: bool,
    use_optimized_kmer: bool,
    name_column: Option<&'a str>,
}

impl<'a> RecordProcessor<'a> {
    fn update_summaries(
        &self,
        summaries: &mut SummaryAccumulator,
        record: &bam::lazy::Record,
        aln_stats_rows: Option<&mut String>,
    ) {
        summaries.total_alns += 1;

        let flags = record.flags().unwrap();
        if flags.is_unmapped() || flags.is_secondary() {
            return;
        }

        let strand_rev = flags.is_reverse_complemented();
        let aln_ref = self.references[record.reference_sequence_id().unwrap().unwrap()]
            .name()
            .as_str();
        if !self.reference_seqs.contains_key(aln_ref) {
            panic!(
                "{} is not found in the input reference sequence names!",
                aln_ref
            );
        }

        let aln_start = usize::from(record.alignment_start().unwrap().unwrap());
        let aln_end = aln_start
            + sam::record::Cigar::try_from(record.cigar())
                .unwrap()
                .alignment_span();

        let generated_intervals = if self.use_optimized_kmer {
            Vec::new()
        } else {
            self.generated_intervals(aln_ref, aln_start, aln_end, strand_rev)
        };
        let mut overlap_intervals = if self.use_optimized_kmer {
            Vec::new()
        } else {
            self.bed_intervals(aln_ref, aln_start, aln_end)
        };
        overlap_intervals.extend(&generated_intervals);
        overlap_intervals.sort();

        let stats = AlnStats::from_record(
            self.references,
            self.reference_seqs,
            record,
            &overlap_intervals,
            self.feature_qual_score_stats,
        );

        summaries.update(&stats);
        if let Some(ref mut kmer_summary) = summaries.kmer_summary {
            kmer_summary.update_record(self.references, self.reference_seqs, record);
        }

        if let Some(rows) = aln_stats_rows {
            self.write_aln_stats_row(rows, &stats);
        }
    }

    fn generated_intervals(
        &self,
        aln_ref: &str,
        aln_start: usize,
        aln_end: usize,
        strand_rev: bool,
    ) -> Vec<FeatureInterval> {
        let mut intervals_vec = Vec::new();
        self.intervals_types
            .iter()
            .for_each(|intervals_type| match intervals_type {
                IntervalsType::Homopolymer => intervals_vec.extend(find_homopolymers(
                    self.reference_seqs[aln_ref].sequence(),
                    aln_start,
                    aln_end,
                    strand_rev,
                )),
                IntervalsType::Window(win_len) => intervals_vec
                    .extend(get_windows(aln_start, aln_end, *win_len, false, strand_rev)),
                IntervalsType::WindowPos(win_len) => intervals_vec
                    .extend(get_windows(aln_start, aln_end, *win_len, true, strand_rev)),
                IntervalsType::Border(win_len) => {
                    intervals_vec.extend(get_borders(aln_start, aln_end, *win_len, strand_rev))
                }
                IntervalsType::Match(seq) => intervals_vec.extend(get_matches(
                    self.reference_seqs[aln_ref].sequence(),
                    aln_start,
                    aln_end,
                    seq,
                    strand_rev,
                )),
                IntervalsType::Kmer(kmer_len) => intervals_vec.extend(get_kmers(
                    self.reference_seqs[aln_ref].sequence(),
                    aln_start,
                    aln_end,
                    *kmer_len,
                    strand_rev,
                )),
                IntervalsType::Bed(_) => {}
            });
        intervals_vec
    }

    fn bed_intervals(
        &self,
        aln_ref: &str,
        aln_start: usize,
        aln_end: usize,
    ) -> Vec<&FeatureInterval> {
        let mut overlap_intervals = Vec::new();
        self.intervals_types
            .iter()
            .for_each(|intervals_type| match intervals_type {
                IntervalsType::Bed(intervals) => {
                    overlap_intervals.extend(intervals.find(aln_ref, aln_start, aln_end))
                }
                _ => {}
            });
        overlap_intervals
    }

    fn write_aln_stats_row(&self, rows: &mut String, stats: &AlnStats) {
        if let Some(name) = self.name_column {
            writeln!(rows, "{},{}", name, stats.to_csv()).unwrap();
        } else {
            writeln!(rows, "{}", stats.to_csv()).unwrap();
        }
    }
}

struct SummaryAccumulator {
    yield_summary: YieldSummary,
    identity_summary: IdentitySummary,
    feature_summary: Option<FeatureSummary>,
    kmer_summary: Option<KmerSummary>,
    cigar_summary: CigarLenSummary,
    bin_summary: Option<BinSummary>,
    qual_score_summary: QualScoreSummary,
    total_alns: usize,
}

impl SummaryAccumulator {
    fn new(
        name_column: Option<String>,
        feature_summary_enabled: bool,
        bin_types: &[BinType],
        feature_qual_score_stats: bool,
        kmer_position_stats: bool,
        kmer_advanced_stats: bool,
        optimized_kmer_lens: Option<&[usize]>,
    ) -> Self {
        Self {
            yield_summary: YieldSummary::new(name_column.clone()),
            identity_summary: IdentitySummary::new(name_column.clone()),
            feature_summary: if feature_summary_enabled {
                Some(FeatureSummary::new_with_q_scores(
                    name_column.clone(),
                    feature_qual_score_stats,
                ))
            } else {
                None
            },
            kmer_summary: optimized_kmer_lens.map(|lens| {
                KmerSummary::new(
                    name_column.clone(),
                    lens.to_vec(),
                    kmer_position_stats,
                    kmer_advanced_stats,
                )
            }),
            cigar_summary: CigarLenSummary::new(name_column.clone()),
            bin_summary: if bin_types.is_empty() {
                None
            } else {
                Some(BinSummary::new(name_column.clone(), bin_types.to_vec()))
            },
            qual_score_summary: QualScoreSummary::new_with_feature_scores(
                name_column,
                feature_qual_score_stats,
            ),
            total_alns: 0,
        }
    }

    fn update(&mut self, stats: &AlnStats) {
        self.yield_summary.update(stats);
        self.identity_summary.update(stats);
        if let Some(ref mut f) = self.feature_summary {
            f.update(stats);
        }
        self.cigar_summary.update(stats);
        if let Some(ref mut b) = self.bin_summary {
            b.update(stats);
        }
        self.qual_score_summary.update(stats);
    }

    fn assign_add(&mut self, o: Self) {
        self.total_alns += o.total_alns;
        self.yield_summary.assign_add(&o.yield_summary);
        self.identity_summary.assign_add(&o.identity_summary);
        if let (Some(f), Some(other_f)) = (&mut self.feature_summary, o.feature_summary) {
            f.assign_add(&other_f);
        }
        if let (Some(k), Some(other_k)) = (&mut self.kmer_summary, o.kmer_summary) {
            k.assign_add(other_k);
        }
        self.cigar_summary.assign_add(&o.cigar_summary);
        if let (Some(b), Some(other_b)) = (&mut self.bin_summary, o.bin_summary) {
            b.assign_add(&other_b);
        }
        self.qual_score_summary.assign_add(&o.qual_score_summary);
    }
}

fn write_summary<D: std::fmt::Display>(s: D, prefix: &str, name: &str) {
    let summary_path = format!("{}.{}", prefix, name);
    let mut summary_writer = File::create(&summary_path).unwrap();
    write!(summary_writer, "{}", s).unwrap();
}

fn main() {
    let start_time = Instant::now();
    let args = Args::parse();

    let bin_types = args
        .bin_types
        .map(|b| b.iter().map(|s| BinType::from_str(s).unwrap()).collect());

    let mut intervals_types = Vec::new();
    if args.intervals_hp {
        intervals_types.push(IntervalsType::Homopolymer);
    }
    if let Some(paths) = args.intervals_bed {
        intervals_types.extend(paths.iter().map(|p| IntervalsType::Bed(Intervals::new(p))));
    }
    if let Some(win_lens) = args.intervals_window {
        intervals_types.extend(win_lens.into_iter().map(|l| IntervalsType::Window(l)));
    }
    if let Some(win_lens) = args.intervals_window_pos {
        intervals_types.extend(win_lens.into_iter().map(|l| IntervalsType::WindowPos(l)));
    }
    if let Some(win_lens) = args.intervals_border {
        intervals_types.extend(win_lens.into_iter().map(|l| IntervalsType::Border(l)));
    }
    if let Some(seqs) = args.intervals_match {
        intervals_types.extend(seqs.into_iter().map(|mut s| {
            s.make_ascii_uppercase();
            IntervalsType::Match(s)
        }));
    }
    if let Some(kmer_lens) = args.intervals_kmer {
        intervals_types.extend(kmer_lens.into_iter().map(IntervalsType::Kmer));
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    run(
        args.input,
        args.reference,
        args.stats_prefix,
        bin_types,
        intervals_types,
        args.name_column,
        !args.no_per_aln_stats,
        !args.no_feature_qual_score_stats,
        args.kmer_position_stats,
        args.kmer_advanced_stats,
        args.record_batch_size,
        args.bam_reader_threads,
    );

    let duration = start_time.elapsed();
    println!("Run time (s): {}", duration.as_secs());
}

enum IntervalsType {
    Bed(Intervals),
    Homopolymer,
    Window(usize),
    WindowPos(usize),
    Border(usize),
    Match(String),
    Kmer(usize),
}

fn optimized_kmer_lens(
    intervals_types: &[IntervalsType],
    feature_qual_score_stats: bool,
) -> Option<Vec<usize>> {
    if feature_qual_score_stats || intervals_types.is_empty() {
        return None;
    }

    let mut kmer_lens = Vec::new();
    for intervals_type in intervals_types {
        match intervals_type {
            IntervalsType::Kmer(kmer_len) if *kmer_len <= MAX_ENCODED_KMER_LEN => {
                kmer_lens.push(*kmer_len)
            }
            _ => return None,
        }
    }

    Some(kmer_lens)
}

fn parse_kmer_len(s: &str) -> Result<usize, String> {
    let len = s
        .parse::<usize>()
        .map_err(|e| format!("invalid k-mer length: {}", e))?;

    if len == 0 {
        Err("k-mer length must be greater than 0".to_string())
    } else {
        Ok(len)
    }
}

fn parse_record_batch_size(s: &str) -> Result<usize, String> {
    let len = s
        .parse::<usize>()
        .map_err(|e| format!("invalid record batch size: {}", e))?;

    if len == 0 {
        Err("record batch size must be greater than 0".to_string())
    } else {
        Ok(len)
    }
}

fn parse_bam_reader_threads(s: &str) -> Result<usize, String> {
    let threads = s
        .parse::<usize>()
        .map_err(|e| format!("invalid BAM reader thread count: {}", e))?;

    if threads == 0 {
        Err("BAM reader thread count must be greater than 0".to_string())
    } else {
        Ok(threads)
    }
}

#[derive(Parser)]
#[clap(author, version, about)]
struct Args {
    /// Input BAM file.
    input: String,

    /// Input reference FASTA file. Can be gzipped.
    reference: String,

    /// Prefix for output files that contain statistics.
    stats_prefix: String,

    /// Add column with a specific name in CSV outputs.
    #[clap(short, long)]
    name_column: Option<String>,

    /// Turn off outputting per alignment stats.
    #[clap(long)]
    no_per_aln_stats: bool,

    /// Turn off interval feature rows in summary_qual_score_stats.csv.
    ///
    /// The all_alignments quality score summary is still reported.
    #[clap(long)]
    no_feature_qual_score_stats: bool,

    /// Write summary_kmer_position_stats.csv for optimized k-mer runs.
    ///
    /// This adds per-offset error counts within each read-oriented k-mer.
    #[clap(long)]
    kmer_position_stats: bool,

    /// Write strand and substitution spectrum summaries for optimized k-mer runs.
    ///
    /// This adds summary_kmer_strand_stats.csv,
    /// summary_kmer_substitution_stats.csv, and
    /// summary_kmer_substitution_profile_stats.csv.
    #[clap(long)]
    kmer_advanced_stats: bool,

    /// Number of records to hand to each Rayon task.
    ///
    /// Larger batches reduce scheduler overhead; smaller batches can improve load balancing.
    #[clap(long, default_value_t = DEFAULT_RECORD_BATCH_SIZE, parse(try_from_str = parse_record_batch_size))]
    record_batch_size: usize,

    /// Number of BGZF decompression workers for the BAM reader.
    ///
    /// These workers are in addition to the Rayon compute threads set by --threads.
    /// Decrease this on low-core machines.
    #[clap(long, default_value_t = DEFAULT_BAM_READER_THREADS, parse(try_from_str = parse_bam_reader_threads))]
    bam_reader_threads: usize,

    /// Types of bins to use for per alignment stats.
    ///
    /// Each bin should be of the format <bin_type>:<step_size>.
    ///
    /// Supported bin types:
    /// q_len (read sequence length),
    /// subread_passes,
    /// mapq,
    /// mean_qual,
    /// gc_content,
    /// concordance_qv (phred scale Q-value)
    #[clap(short, long, min_values = 1)]
    bin_types: Option<Vec<String>>,

    /// Use intervals from a BED file.
    ///
    /// The BED file should have the columns chrom, start, stop, and feature.
    /// The feature column is optional.
    ///
    /// This allows stats to be gathered separately for different types of intervals.
    /// Note that all intervals are on the reference, not the reads.
    #[clap(long, min_values = 1)]
    intervals_bed: Option<Vec<String>>,

    /// Use homopolymer regions in the reference as intervals.
    #[clap(long)]
    intervals_hp: bool,

    /// Use fixed-width nonoverlapping windows as intervals.
    ///
    /// This is used to specify the window widths.
    #[clap(long, min_values = 1)]
    intervals_window: Option<Vec<usize>>,

    /// Use fixed-width nonoverlapping windows with positions as intervals.
    ///
    /// This is used to specify the window widths.
    #[clap(long, min_values = 1)]
    intervals_window_pos: Option<Vec<usize>>,

    /// Use fixed-width nonoverlapping window border regions as intervals.
    ///
    /// This is used to specify the window widths.
    #[clap(long, min_values = 1)]
    intervals_border: Option<Vec<usize>>,

    /// Use regions in the reference that match any of the specified subsequences as intervals.
    #[clap(long, min_values = 1)]
    intervals_match: Option<Vec<String>>,

    /// Use all A/C/G/T k-mers of the specified lengths as intervals.
    ///
    /// K-mers are taken from the aligned reference span and reported in read orientation.
    #[clap(long, min_values = 1, parse(try_from_str = parse_kmer_len))]
    intervals_kmer: Option<Vec<usize>>,

    /// Number of threads. Will be automatically determined if this is set to 0.
    #[clap(short, long, default_value_t = 0usize)]
    threads: usize,
}
