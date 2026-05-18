//! Optimized k-mer context summaries.
//!
//! This module backs the k-mer-only path used when `--intervals-kmer` is run
//! with `--no-feature-qual-score-stats`. K-mers are encoded as two bits per
//! A/C/G/T base and stored in read orientation, so reverse-strand alignments use
//! the rolling reverse-complement encoding. Small k values use a dense direct
//! index from encoded k-mer bits to stats rows; larger k values fall back to a
//! sparse hash map while still avoiding string allocation in the hot path.
//!
//! Position-level summaries are optional because they multiply the number of
//! counters by k and add work to every overlapping k-mer update.

use std::{collections::BTreeMap, fmt};

use fxhash::FxHashMap;
use noodles::{bam, core::Position, fasta, sam};

use sam::record::cigar::op::Kind;

use crate::stats::{concordance_qv, qual_to_error, FeatureStats};

pub const MAX_ENCODED_KMER_LEN: usize = 32;

type KmerKey = (usize, u64);
const INVALID_STAT_IDX: usize = usize::MAX;
const MAX_DENSE_KMER_STATES: usize = 1 << 18;
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

struct KmerWindowContext {
    kmer_len: usize,
    start: usize,
    stat_idxs: Vec<usize>,
    has_error: Vec<bool>,
}

impl KmerWindowContext {
    fn new(kmer_len: usize, start: usize, num_windows: usize) -> Self {
        Self {
            kmer_len,
            start,
            stat_idxs: vec![INVALID_STAT_IDX; num_windows],
            has_error: vec![false; num_windows],
        }
    }

    fn last_start(&self) -> usize {
        self.start + self.stat_idxs.len() - 1
    }

    fn overlapping_idx_bounds(&self, ref_pos: usize) -> Option<(usize, usize)> {
        if self.stat_idxs.is_empty() {
            return None;
        }

        let lo = ref_pos.saturating_sub(self.kmer_len - 1).max(self.start);
        let hi = ref_pos.min(self.last_start());
        if lo > hi {
            None
        } else {
            Some((lo - self.start, hi - self.start))
        }
    }

    fn read_oriented_offset(&self, idx: usize, ref_pos: usize, strand_rev: bool) -> usize {
        let offset = ref_pos - (self.start + idx);
        if strand_rev {
            self.kmer_len - 1 - offset
        } else {
            offset
        }
    }
}

#[derive(Default)]
struct KmerPositionStats {
    mismatches: usize,
    non_hp_ins: usize,
    non_hp_del: usize,
    hp_ins: usize,
    hp_del: usize,
    skips: usize,
}

impl KmerPositionStats {
    fn assign_add_ref(&mut self, o: &Self) {
        self.mismatches += o.mismatches;
        self.non_hp_ins += o.non_hp_ins;
        self.non_hp_del += o.non_hp_del;
        self.hp_ins += o.hp_ins;
        self.hp_del += o.hp_del;
        self.skips += o.skips;
    }

    fn assign_add(&mut self, o: Self) {
        self.assign_add_ref(&o);
    }

    fn num_errors(&self) -> usize {
        self.mismatches + self.non_hp_ins + self.hp_ins + self.non_hp_del + self.hp_del
    }

    fn matches(&self, intervals: usize) -> usize {
        intervals.saturating_sub(self.mismatches + self.non_hp_del + self.hp_del + self.skips)
    }

    fn num_observations(&self, intervals: usize) -> usize {
        self.matches(intervals) + self.num_errors()
    }

    fn identity(&self, intervals: usize) -> f64 {
        let observations = self.num_observations(intervals);
        if observations == 0 {
            0.0
        } else {
            (self.matches(intervals) as f64) / (observations as f64)
        }
    }

    fn error_rate(&self, intervals: usize) -> f64 {
        let observations = self.num_observations(intervals);
        if observations == 0 {
            0.0
        } else {
            (self.num_errors() as f64) / (observations as f64)
        }
    }
}

struct KmerStats {
    aggregate: FeatureStats,
    positions: Vec<KmerPositionStats>,
    strands: Vec<FeatureStats>,
    substitutions: Vec<[usize; 16]>,
}

impl KmerStats {
    fn new(kmer_len: usize, track_position_stats: bool, track_advanced_stats: bool) -> Self {
        let positions = if track_position_stats {
            (0..kmer_len)
                .map(|_| KmerPositionStats::default())
                .collect()
        } else {
            Vec::new()
        };
        let strands = if track_advanced_stats {
            vec![FeatureStats::new(false), FeatureStats::new(false)]
        } else {
            Vec::new()
        };
        let substitutions = if track_advanced_stats {
            vec![[0usize; 16]; kmer_len]
        } else {
            Vec::new()
        };
        Self {
            aggregate: FeatureStats::new(false),
            positions,
            strands,
            substitutions,
        }
    }

    fn assign_add(&mut self, o: Self) {
        self.aggregate.assign_add(&o.aggregate);
        self.positions
            .iter_mut()
            .zip(o.positions)
            .for_each(|(stats, other_stats)| stats.assign_add(other_stats));
        self.strands
            .iter_mut()
            .zip(o.strands)
            .for_each(|(stats, other_stats)| stats.assign_add(&other_stats));
        self.substitutions
            .iter_mut()
            .zip(o.substitutions)
            .for_each(|(counts, other_counts)| {
                counts
                    .iter_mut()
                    .zip(other_counts)
                    .for_each(|(count, other_count)| *count += other_count);
            });
    }
}

pub struct KmerSummary {
    name_column: Option<String>,
    kmer_lens: Vec<usize>,
    kmer_index: FxHashMap<KmerKey, usize>,
    dense_kmer_index: FxHashMap<usize, Vec<usize>>,
    kmer_stats: Vec<(KmerKey, KmerStats)>,
    track_position_stats: bool,
    track_advanced_stats: bool,
}

struct KmerContextAggregate {
    kmers: usize,
    stats: FeatureStats,
}

impl KmerContextAggregate {
    fn new() -> Self {
        Self {
            kmers: 0,
            stats: FeatureStats::new(false),
        }
    }
}

#[derive(Default)]
struct KmerPositionProfileAggregate {
    kmers: usize,
    intervals: usize,
    stats: KmerPositionStats,
}

#[derive(Default)]
struct KmerSubstitutionProfileAggregate {
    kmers: usize,
    intervals: usize,
    count: usize,
}

impl KmerSummary {
    pub fn new(
        mut name_column: Option<String>,
        kmer_lens: Vec<usize>,
        track_position_stats: bool,
        track_advanced_stats: bool,
    ) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        Self {
            name_column,
            dense_kmer_index: dense_kmer_index(&kmer_lens),
            kmer_lens,
            kmer_index: FxHashMap::default(),
            kmer_stats: Vec::new(),
            track_position_stats,
            track_advanced_stats,
        }
    }

    pub fn update_record(
        &mut self,
        references: &sam::header::ReferenceSequences,
        reference_seqs: &FxHashMap<String, fasta::Record>,
        r: &bam::lazy::Record,
    ) {
        let flags = r.flags().unwrap();
        if flags.is_supplementary() {
            return;
        }

        let chr = references[r.reference_sequence_id().unwrap().unwrap()]
            .name()
            .to_string();
        let curr_ref_seq = reference_seqs[&chr].sequence().as_ref();
        let strand_rev = flags.is_reverse_complemented();
        let aln_start = usize::from(r.alignment_start().unwrap().unwrap());
        let cigar = sam::record::Cigar::try_from(r.cigar()).unwrap();
        let aln_end = aln_start + cigar.alignment_span();
        let mut contexts = self.get_contexts(curr_ref_seq, aln_start, aln_end, strand_rev);
        if contexts.is_empty() {
            return;
        }

        let sequence = sam::record::Sequence::try_from(r.sequence()).unwrap();
        let q_scores = sam::record::QualityScores::try_from(r.quality_scores()).unwrap();
        if self.track_position_stats {
            self.update_record_with_positions(
                curr_ref_seq,
                &sequence,
                &q_scores,
                &cigar,
                aln_start,
                strand_rev,
                &mut contexts,
            );
        } else {
            self.update_record_aggregate(
                curr_ref_seq,
                &sequence,
                &q_scores,
                &cigar,
                aln_start,
                strand_rev,
                &mut contexts,
            );
        }
    }

    fn update_record_aggregate(
        &mut self,
        curr_ref_seq: &[u8],
        sequence: &sam::record::Sequence,
        q_scores: &sam::record::QualityScores,
        cigar: &sam::record::Cigar,
        mut ref_pos: usize,
        strand_rev: bool,
        contexts: &mut [KmerWindowContext],
    ) {
        let mut query_pos = 1;
        let strand_idx = strand_idx(strand_rev);

        for op in cigar.iter() {
            for _i in 0..op.len() {
                match op.kind() {
                    Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match => {
                        let c = ref_base(curr_ref_seq, ref_pos);
                        let query_base = u8::from(sequence[Position::new(query_pos).unwrap()])
                            .to_ascii_uppercase();
                        let is_match = op.kind() == Kind::SequenceMatch
                            || (op.kind() == Kind::Match && c == query_base);
                        let q_score = u8::from(q_scores[Position::new(query_pos).unwrap()]);
                        let qual_error = qual_to_error(q_score);

                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        stats.aggregate.total_qual_error += qual_error;
                                        if self.track_advanced_stats {
                                            stats.strands[strand_idx].total_qual_error +=
                                                qual_error;
                                        }
                                        if is_match {
                                            stats.aggregate.matches += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].matches += 1;
                                            }
                                        } else {
                                            stats.aggregate.mismatches += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].mismatches += 1;
                                                let offset = context
                                                    .read_oriented_offset(idx, ref_pos, strand_rev);
                                                if let Some(sub_idx) = substitution_idx(
                                                    read_oriented_base(c, strand_rev),
                                                    read_oriented_base(query_base, strand_rev),
                                                ) {
                                                    stats.substitutions[offset][sub_idx] += 1;
                                                }
                                            }
                                            context.has_error[idx] = true;
                                        }
                                    }
                                }
                            }
                        }

                        query_pos += 1;
                        ref_pos += 1;
                    }
                    Kind::Insertion => {
                        // BEST attributes insertions to the current reference position.
                        // This matches the existing interval-feature convention.
                        let before_ins = ref_base(curr_ref_seq, ref_pos);
                        let after_ins = ref_base_opt(curr_ref_seq, ref_pos + 1).unwrap_or(b'?');
                        let query_ins = &sequence[Position::new(query_pos).unwrap()
                            ..Position::new(query_pos + op.len()).unwrap()];
                        let hp_before = query_ins
                            .iter()
                            .map(|&c| u8::from(c).to_ascii_uppercase())
                            .all(|c| c == before_ins);
                        let hp_after = query_ins
                            .iter()
                            .map(|&c| u8::from(c).to_ascii_uppercase())
                            .all(|c| c == after_ins);

                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        if hp_before || hp_after {
                                            stats.aggregate.hp_ins += op.len();
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].hp_ins += op.len();
                                            }
                                        } else {
                                            stats.aggregate.non_hp_ins += op.len();
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].non_hp_ins += op.len();
                                            }
                                        }
                                        context.has_error[idx] = true;
                                    }
                                }
                            }
                        }

                        query_pos += op.len();
                        break;
                    }
                    Kind::Deletion => {
                        // Deletions are counted once per deleted reference base and per
                        // overlapping k-mer window.
                        let before_curr = ref_base_opt(curr_ref_seq, ref_pos - 1).unwrap_or(b'?');
                        let after_curr = ref_base_opt(curr_ref_seq, ref_pos + 1).unwrap_or(b'?');
                        let curr = ref_base(curr_ref_seq, ref_pos);
                        let hp = curr == before_curr || curr == after_curr;

                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        if hp {
                                            stats.aggregate.hp_del += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].hp_del += 1;
                                            }
                                        } else {
                                            stats.aggregate.non_hp_del += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].non_hp_del += 1;
                                            }
                                        }
                                        context.has_error[idx] = true;
                                    }
                                }
                            }
                        }

                        ref_pos += 1;
                    }
                    Kind::SoftClip => {
                        query_pos += op.len();
                        break;
                    }
                    Kind::HardClip => {
                        break;
                    }
                    Kind::Skip => {
                        ref_pos += 1;
                    }
                    _ => panic!("Unexpected CIGAR operation: {}", op),
                }
            }
        }

        self.add_identical_overlaps(contexts, strand_idx);
    }

    fn update_record_with_positions(
        &mut self,
        curr_ref_seq: &[u8],
        sequence: &sam::record::Sequence,
        q_scores: &sam::record::QualityScores,
        cigar: &sam::record::Cigar,
        mut ref_pos: usize,
        strand_rev: bool,
        contexts: &mut [KmerWindowContext],
    ) {
        let mut query_pos = 1;
        let strand_idx = strand_idx(strand_rev);

        for op in cigar.iter() {
            for _i in 0..op.len() {
                match op.kind() {
                    Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match => {
                        let c = ref_base(curr_ref_seq, ref_pos);
                        let query_base = u8::from(sequence[Position::new(query_pos).unwrap()])
                            .to_ascii_uppercase();
                        let is_match = op.kind() == Kind::SequenceMatch
                            || (op.kind() == Kind::Match && c == query_base);
                        let q_score = u8::from(q_scores[Position::new(query_pos).unwrap()]);
                        let qual_error = qual_to_error(q_score);

                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        stats.aggregate.total_qual_error += qual_error;
                                        if self.track_advanced_stats {
                                            stats.strands[strand_idx].total_qual_error +=
                                                qual_error;
                                        }
                                        if is_match {
                                            stats.aggregate.matches += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].matches += 1;
                                            }
                                        } else {
                                            stats.aggregate.mismatches += 1;
                                            let offset = context
                                                .read_oriented_offset(idx, ref_pos, strand_rev);
                                            stats.positions[offset].mismatches += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].mismatches += 1;
                                                if let Some(sub_idx) = substitution_idx(
                                                    read_oriented_base(c, strand_rev),
                                                    read_oriented_base(query_base, strand_rev),
                                                ) {
                                                    stats.substitutions[offset][sub_idx] += 1;
                                                }
                                            }
                                            context.has_error[idx] = true;
                                        }
                                    }
                                }
                            }
                        }

                        query_pos += 1;
                        ref_pos += 1;
                    }
                    Kind::Insertion => {
                        // Position stats use the same current-reference-position
                        // insertion attribution as the aggregate interval stats.
                        let before_ins = ref_base(curr_ref_seq, ref_pos);
                        let after_ins = ref_base_opt(curr_ref_seq, ref_pos + 1).unwrap_or(b'?');
                        let query_ins = &sequence[Position::new(query_pos).unwrap()
                            ..Position::new(query_pos + op.len()).unwrap()];
                        let hp_before = query_ins
                            .iter()
                            .map(|&c| u8::from(c).to_ascii_uppercase())
                            .all(|c| c == before_ins);
                        let hp_after = query_ins
                            .iter()
                            .map(|&c| u8::from(c).to_ascii_uppercase())
                            .all(|c| c == after_ins);

                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        let offset =
                                            context.read_oriented_offset(idx, ref_pos, strand_rev);
                                        if hp_before || hp_after {
                                            stats.aggregate.hp_ins += op.len();
                                            stats.positions[offset].hp_ins += op.len();
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].hp_ins += op.len();
                                            }
                                        } else {
                                            stats.aggregate.non_hp_ins += op.len();
                                            stats.positions[offset].non_hp_ins += op.len();
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].non_hp_ins += op.len();
                                            }
                                        }
                                        context.has_error[idx] = true;
                                    }
                                }
                            }
                        }

                        query_pos += op.len();
                        break;
                    }
                    Kind::Deletion => {
                        // The reported offset is in the read-oriented k-mer, so reverse
                        // alignments flip the reference offset within the window.
                        let before_curr = ref_base_opt(curr_ref_seq, ref_pos - 1).unwrap_or(b'?');
                        let after_curr = ref_base_opt(curr_ref_seq, ref_pos + 1).unwrap_or(b'?');
                        let curr = ref_base(curr_ref_seq, ref_pos);
                        let hp = curr == before_curr || curr == after_curr;

                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        let offset =
                                            context.read_oriented_offset(idx, ref_pos, strand_rev);
                                        if hp {
                                            stats.aggregate.hp_del += 1;
                                            stats.positions[offset].hp_del += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].hp_del += 1;
                                            }
                                        } else {
                                            stats.aggregate.non_hp_del += 1;
                                            stats.positions[offset].non_hp_del += 1;
                                            if self.track_advanced_stats {
                                                stats.strands[strand_idx].non_hp_del += 1;
                                            }
                                        }
                                        context.has_error[idx] = true;
                                    }
                                }
                            }
                        }

                        ref_pos += 1;
                    }
                    Kind::SoftClip => {
                        query_pos += op.len();
                        break;
                    }
                    Kind::HardClip => {
                        break;
                    }
                    Kind::Skip => {
                        for context in contexts.iter_mut() {
                            if let Some((lo, hi)) = context.overlapping_idx_bounds(ref_pos) {
                                for idx in lo..=hi {
                                    let stat_idx = context.stat_idxs[idx];
                                    if stat_idx != INVALID_STAT_IDX {
                                        let stats = &mut self.kmer_stats[stat_idx].1;
                                        let offset =
                                            context.read_oriented_offset(idx, ref_pos, strand_rev);
                                        stats.positions[offset].skips += 1;
                                    }
                                }
                            }
                        }

                        ref_pos += 1;
                    }
                    _ => panic!("Unexpected CIGAR operation: {}", op),
                }
            }
        }

        self.add_identical_overlaps(contexts, strand_idx);
    }

    fn add_identical_overlaps(&mut self, contexts: &[KmerWindowContext], strand_idx: usize) {
        for context in contexts {
            for (idx, stat_idx) in context.stat_idxs.iter().enumerate() {
                if *stat_idx != INVALID_STAT_IDX {
                    if !context.has_error[idx] {
                        let stats = &mut self.kmer_stats[*stat_idx].1;
                        stats.aggregate.identical_overlaps += 1;
                        if self.track_advanced_stats {
                            stats.strands[strand_idx].identical_overlaps += 1;
                        }
                    }
                }
            }
        }
    }

    fn get_or_insert_stat_idx(&mut self, key: KmerKey) -> usize {
        if self.dense_kmer_index.contains_key(&key.0) {
            let dense_idx = key.1 as usize;
            let stat_idx = self.dense_kmer_index[&key.0][dense_idx];
            if stat_idx != INVALID_STAT_IDX {
                return stat_idx;
            }

            let stat_idx = self.push_stat(key);
            self.dense_kmer_index.get_mut(&key.0).unwrap()[dense_idx] = stat_idx;
            return stat_idx;
        }

        self.get_or_insert_sparse_stat_idx(key)
    }

    fn get_or_insert_sparse_stat_idx(&mut self, key: KmerKey) -> usize {
        if let Some(&stat_idx) = self.kmer_index.get(&key) {
            stat_idx
        } else {
            let stat_idx = self.push_stat(key);
            self.kmer_index.insert(key, stat_idx);
            stat_idx
        }
    }

    fn push_stat(&mut self, key: KmerKey) -> usize {
        let stat_idx = self.kmer_stats.len();
        self.kmer_stats.push((
            key,
            KmerStats::new(key.0, self.track_position_stats, self.track_advanced_stats),
        ));
        stat_idx
    }

    pub fn assign_add(&mut self, o: Self) {
        for (k, v) in o.kmer_stats {
            let stat_idx = self.get_or_insert_stat_idx(k);
            self.kmer_stats[stat_idx].1.assign_add(v);
        }
    }

    pub fn feature_summary(&self) -> KmerFeatureSummary {
        KmerFeatureSummary { summary: self }
    }

    pub fn kmer_summary(&self) -> KmerStatsSummary {
        KmerStatsSummary { summary: self }
    }

    pub fn length_summary(&self) -> KmerLengthSummary {
        KmerLengthSummary { summary: self }
    }

    pub fn context_summary(&self) -> KmerContextSummary {
        KmerContextSummary { summary: self }
    }

    pub fn strand_summary(&self) -> KmerStrandSummary {
        KmerStrandSummary { summary: self }
    }

    pub fn substitution_summary(&self) -> KmerSubstitutionSummary {
        KmerSubstitutionSummary { summary: self }
    }

    pub fn substitution_profile_summary(&self) -> KmerSubstitutionProfileSummary {
        KmerSubstitutionProfileSummary { summary: self }
    }

    pub fn position_summary(&self) -> KmerPositionSummary {
        KmerPositionSummary { summary: self }
    }

    pub fn position_profile_summary(&self) -> KmerPositionProfileSummary {
        KmerPositionProfileSummary { summary: self }
    }

    fn get_contexts(
        &mut self,
        seq: &[u8],
        start: usize,
        end: usize,
        strand_rev: bool,
    ) -> Vec<KmerWindowContext> {
        let mut contexts = Vec::new();

        if end <= start {
            return contexts;
        }

        let bounded_end = end.min(seq.len() + 1);
        if bounded_end <= start {
            return contexts;
        }

        let kmer_lens = self.kmer_lens.clone();
        for kmer_len in kmer_lens {
            if kmer_len == 0 || kmer_len > MAX_ENCODED_KMER_LEN {
                continue;
            }
            if bounded_end - start < kmer_len {
                continue;
            }

            let last_start = bounded_end - kmer_len;
            let mut context = KmerWindowContext::new(kmer_len, start, last_start - start + 1);
            let mut fwd_bits = 0u64;
            let mut rev_bits = 0u64;
            let mut valid_len = 0usize;
            let mask = if kmer_len == MAX_ENCODED_KMER_LEN {
                u64::MAX
            } else {
                (1u64 << (2 * kmer_len)) - 1
            };
            let rev_shift = 2 * (kmer_len - 1);
            let mut has_valid_kmers = false;
            let mut dense_index = self.dense_kmer_index.remove(&kmer_len);

            for ref_pos in start..bounded_end {
                if let Some(bits) = encode_base(seq[ref_pos - 1]) {
                    // Maintain both forward and reverse-complement rolling encodings.
                    // Ambiguous bases reset the window, so no k-mer crossing an N is used.
                    fwd_bits = ((fwd_bits << 2) | bits) & mask;
                    rev_bits = (rev_bits >> 2) | ((bits ^ 0b11) << rev_shift);
                    valid_len += 1;

                    if valid_len >= kmer_len {
                        let window_start = ref_pos + 1 - kmer_len;
                        if window_start >= start {
                            let key_bits = if strand_rev { rev_bits } else { fwd_bits };
                            let idx = window_start - start;
                            let stat_idx = if let Some(ref mut dense_index) = dense_index {
                                // Dense lookup avoids hashing for small k values where
                                // 4^k states fit comfortably in memory.
                                let dense_idx = key_bits as usize;
                                let existing_idx = dense_index[dense_idx];
                                if existing_idx == INVALID_STAT_IDX {
                                    let new_idx = self.push_stat((kmer_len, key_bits));
                                    dense_index[dense_idx] = new_idx;
                                    new_idx
                                } else {
                                    existing_idx
                                }
                            } else {
                                self.get_or_insert_sparse_stat_idx((kmer_len, key_bits))
                            };
                            context.stat_idxs[idx] = stat_idx;
                            let stats = &mut self.kmer_stats[stat_idx].1;
                            stats.aggregate.overlaps += 1;
                            if self.track_advanced_stats {
                                stats.strands[strand_idx(strand_rev)].overlaps += 1;
                            }
                            has_valid_kmers = true;
                        }
                    }
                } else {
                    fwd_bits = 0;
                    rev_bits = 0;
                    valid_len = 0;
                }
            }

            if let Some(dense_index) = dense_index {
                self.dense_kmer_index.insert(kmer_len, dense_index);
            }
            if has_valid_kmers {
                contexts.push(context);
            }
        }

        contexts
    }
}

pub struct KmerFeatureSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerStatsSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerLengthSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerContextSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerStrandSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerSubstitutionSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerSubstitutionProfileSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerPositionSummary<'a> {
    summary: &'a KmerSummary,
}

pub struct KmerPositionProfileSummary<'a> {
    summary: &'a KmerSummary,
}

fn write_aggregate_metrics(f: &mut fmt::Formatter, stats: &FeatureStats) -> fmt::Result {
    let per_interval = |x| (x as f64) / (stats.overlaps as f64);
    let id = stats.identity();
    write!(
        f,
        "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
        stats.overlaps,
        per_interval(stats.identical_overlaps),
        id,
        concordance_qv(id, id != 1.0),
        stats.mean_qual(),
        per_interval(stats.num_bases()),
        per_interval(stats.matches),
        per_interval(stats.mismatches),
        per_interval(stats.non_hp_ins),
        per_interval(stats.non_hp_del),
        per_interval(stats.hp_ins),
        per_interval(stats.hp_del)
    )
}

impl fmt::Display for KmerFeatureSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}feature,intervals,identical_intervals,identity,identity_qv,mean_qual,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut v = self
            .summary
            .kmer_stats
            .iter()
            .map(|((k, bits), stats)| (decode_kmer(*k, *bits), &stats.aggregate))
            .collect::<Vec<_>>();
        v.sort_by(|a, b| a.0.cmp(&b.0));

        for (feature, stats) in v.into_iter() {
            let per_interval = |x| (x as f64) / (stats.overlaps as f64);
            let id = stats.identity();
            writeln!(
                f,
                "{}{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                feature,
                stats.overlaps,
                per_interval(stats.identical_overlaps),
                id,
                concordance_qv(id, id != 1.0),
                stats.mean_qual(),
                per_interval(stats.num_bases()),
                per_interval(stats.matches),
                per_interval(stats.mismatches),
                per_interval(stats.non_hp_ins),
                per_interval(stats.non_hp_del),
                per_interval(stats.hp_ins),
                per_interval(stats.hp_del)
            )?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerStatsSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,kmer,intervals,identical_intervals,identity,identity_qv,mean_qual,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut v = self
            .summary
            .kmer_stats
            .iter()
            .map(|((k, bits), stats)| (*k, decode_kmer(*k, *bits), &stats.aggregate))
            .collect::<Vec<_>>();
        v.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));

        for (kmer_len, kmer, stats) in v.into_iter() {
            let per_interval = |x| (x as f64) / (stats.overlaps as f64);
            let id = stats.identity();
            writeln!(
                f,
                "{}{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                kmer,
                stats.overlaps,
                per_interval(stats.identical_overlaps),
                id,
                concordance_qv(id, id != 1.0),
                stats.mean_qual(),
                per_interval(stats.num_bases()),
                per_interval(stats.matches),
                per_interval(stats.mismatches),
                per_interval(stats.non_hp_ins),
                per_interval(stats.non_hp_del),
                per_interval(stats.hp_ins),
                per_interval(stats.hp_del)
            )?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerLengthSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,kmers,intervals,identical_intervals,identity,identity_qv,mean_qual,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut rows: BTreeMap<usize, KmerContextAggregate> = BTreeMap::new();
        for ((k, _bits), stats) in &self.summary.kmer_stats {
            let row = rows.entry(*k).or_insert_with(KmerContextAggregate::new);
            row.kmers += 1;
            row.stats.assign_add(&stats.aggregate);
        }

        for (kmer_len, row) in rows {
            write!(
                f,
                "{}{},{},",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                row.kmers
            )?;
            write_aggregate_metrics(f, &row.stats)?;
            writeln!(f)?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerContextSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,context_type,context_value,kmers,intervals,identical_intervals,identity,identity_qv,mean_qual,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut rows: BTreeMap<(usize, String, String), KmerContextAggregate> = BTreeMap::new();
        for ((k, bits), stats) in &self.summary.kmer_stats {
            let kmer = decode_kmer(*k, *bits);
            for (context_type, context_value) in kmer_contexts(&kmer) {
                let row = rows
                    .entry((*k, context_type.to_string(), context_value))
                    .or_insert_with(KmerContextAggregate::new);
                row.kmers += 1;
                row.stats.assign_add(&stats.aggregate);
            }
        }

        for ((kmer_len, context_type, context_value), row) in rows {
            write!(
                f,
                "{}{},{},{},{},",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                context_type,
                context_value,
                row.kmers
            )?;
            write_aggregate_metrics(f, &row.stats)?;
            writeln!(f)?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerStrandSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,kmer,strand,intervals,identical_intervals,identity,identity_qv,mean_qual,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut rows = Vec::new();
        for ((k, bits), stats) in &self.summary.kmer_stats {
            let kmer = decode_kmer(*k, *bits);
            for (strand_idx, strand_stats) in stats.strands.iter().enumerate() {
                if strand_stats.overlaps > 0 {
                    rows.push((*k, kmer.clone(), strand_idx, strand_stats));
                }
            }
        }
        rows.sort_by(|a, b| (a.0, &a.1, a.2).cmp(&(b.0, &b.1, b.2)));

        for (kmer_len, kmer, strand_idx, stats) in rows {
            write!(
                f,
                "{}{},{},{},",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                kmer,
                strand_label(strand_idx)
            )?;
            write_aggregate_metrics(f, stats)?;
            writeln!(f)?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerSubstitutionSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,kmer,offset,ref_base,read_base,count,rate_per_interval,fraction_of_kmer_mismatches", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut rows = Vec::new();
        for ((k, bits), stats) in &self.summary.kmer_stats {
            let kmer = decode_kmer(*k, *bits);
            let total_substitutions = stats
                .substitutions
                .iter()
                .flat_map(|counts| counts.iter())
                .sum::<usize>();
            for (offset, counts) in stats.substitutions.iter().enumerate() {
                for (sub_idx, count) in counts.iter().enumerate() {
                    if *count > 0 {
                        let (ref_base, read_base) = substitution_bases(sub_idx);
                        rows.push((
                            *k,
                            kmer.clone(),
                            offset,
                            ref_base,
                            read_base,
                            *count,
                            stats.aggregate.overlaps,
                            total_substitutions,
                        ));
                    }
                }
            }
        }
        rows.sort_by(|a, b| (a.0, &a.1, a.2, a.3, a.4).cmp(&(b.0, &b.1, b.2, b.3, b.4)));

        for (kmer_len, kmer, offset, ref_base, read_base, count, intervals, total_substitutions) in
            rows
        {
            let rate = (count as f64) / (intervals as f64);
            let fraction = if total_substitutions == 0 {
                0.0
            } else {
                (count as f64) / (total_substitutions as f64)
            };
            writeln!(
                f,
                "{}{},{},{},{},{},{},{:.8},{:.8}",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                kmer,
                offset,
                ref_base as char,
                read_base as char,
                count,
                rate,
                fraction
            )?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerSubstitutionProfileSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}kmer_len,offset,ref_base,read_base,kmers,intervals,count,rate_per_interval",
            if self.summary.name_column.is_some() {
                "name,"
            } else {
                ""
            }
        )?;
        let mut rows: BTreeMap<(usize, usize, u8, u8), KmerSubstitutionProfileAggregate> =
            BTreeMap::new();
        for ((k, bits), stats) in &self.summary.kmer_stats {
            let kmer = decode_kmer(*k, *bits);
            for (offset, counts) in stats.substitutions.iter().enumerate() {
                let ref_base = kmer.as_bytes()[offset];
                for read_base in BASES {
                    if read_base == ref_base {
                        continue;
                    }
                    let Some(sub_idx) = substitution_idx(ref_base, read_base) else {
                        continue;
                    };
                    let row = rows.entry((*k, offset, ref_base, read_base)).or_default();
                    row.kmers += 1;
                    row.intervals += stats.aggregate.overlaps;
                    row.count += counts[sub_idx];
                }
            }
        }

        for ((kmer_len, offset, ref_base, read_base), row) in rows {
            let rate = if row.intervals == 0 {
                0.0
            } else {
                (row.count as f64) / (row.intervals as f64)
            };
            writeln!(
                f,
                "{}{},{},{},{},{},{},{},{:.8}",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                offset,
                ref_base as char,
                read_base as char,
                row.kmers,
                row.intervals,
                row.count,
                rate
            )?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerPositionSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,kmer,offset,kmer_base,intervals,identity,identity_qv,matches,mismatches,non_hp_ins,non_hp_del,hp_ins,hp_del,skips,error_rate", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut v = Vec::new();
        for ((k, bits), kmer_stats) in &self.summary.kmer_stats {
            let kmer = decode_kmer(*k, *bits);
            for (offset, stats) in kmer_stats.positions.iter().enumerate() {
                v.push((
                    *k,
                    kmer.clone(),
                    offset,
                    kmer_stats.aggregate.overlaps,
                    stats,
                ));
            }
        }
        v.sort_by(|a, b| (a.0, &a.1, a.2).cmp(&(b.0, &b.1, b.2)));

        for (kmer_len, kmer, offset, intervals, stats) in v.into_iter() {
            let id = stats.identity(intervals);
            let matches = stats.matches(intervals);
            let base = kmer.as_bytes()[offset] as char;
            writeln!(
                f,
                "{}{},{},{},{},{},{:.6},{:.6},{},{},{},{},{},{},{},{:.6}",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                kmer,
                offset,
                base,
                intervals,
                id,
                concordance_qv(id, stats.num_errors() > 0),
                matches,
                stats.mismatches,
                stats.non_hp_ins,
                stats.non_hp_del,
                stats.hp_ins,
                stats.hp_del,
                stats.skips,
                stats.error_rate(intervals)
            )?;
        }

        Ok(())
    }
}

impl fmt::Display for KmerPositionProfileSummary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}kmer_len,offset,kmer_base,kmers,intervals,identity,identity_qv,matches,mismatches,non_hp_ins,non_hp_del,hp_ins,hp_del,skips,error_rate", if self.summary.name_column.is_some() { "name," } else { "" })?;
        let mut rows: BTreeMap<(usize, usize, String), KmerPositionProfileAggregate> =
            BTreeMap::new();

        for ((k, bits), kmer_stats) in &self.summary.kmer_stats {
            let kmer = decode_kmer(*k, *bits);
            for (offset, stats) in kmer_stats.positions.iter().enumerate() {
                let base = (kmer.as_bytes()[offset] as char).to_string();
                for key_base in ["*".to_string(), base] {
                    let row = rows.entry((*k, offset, key_base)).or_default();
                    row.kmers += 1;
                    row.intervals += kmer_stats.aggregate.overlaps;
                    row.stats.assign_add_ref(stats);
                }
            }
        }

        for ((kmer_len, offset, kmer_base), row) in rows {
            let id = row.stats.identity(row.intervals);
            let matches = row.stats.matches(row.intervals);
            writeln!(
                f,
                "{}{},{},{},{},{},{:.6},{:.6},{},{},{},{},{},{},{},{:.6}",
                self.summary
                    .name_column
                    .as_ref()
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                kmer_len,
                offset,
                kmer_base,
                row.kmers,
                row.intervals,
                id,
                concordance_qv(id, row.stats.num_errors() > 0),
                matches,
                row.stats.mismatches,
                row.stats.non_hp_ins,
                row.stats.non_hp_del,
                row.stats.hp_ins,
                row.stats.hp_del,
                row.stats.skips,
                row.stats.error_rate(row.intervals)
            )?;
        }

        Ok(())
    }
}

fn kmer_contexts(kmer: &str) -> Vec<(&'static str, String)> {
    let bases = kmer.as_bytes();
    let gc_count = bases.iter().filter(|&&base| is_gc(base)).count();
    let mut contexts = vec![
        ("gc_count", gc_count.to_string()),
        (
            "max_homopolymer_run",
            max_homopolymer_run(bases).to_string(),
        ),
        (
            "reverse_complement_palindrome",
            is_reverse_complement_palindrome(bases).to_string(),
        ),
        (
            "edge_bases",
            format!("{}{}", bases[0] as char, bases[bases.len() - 1] as char),
        ),
    ];

    if bases.len() % 2 == 1 {
        contexts.push(("central_base", (bases[bases.len() / 2] as char).to_string()));
    } else {
        let mid = bases.len() / 2;
        contexts.push((
            "central_dinucleotide",
            format!("{}{}", bases[mid - 1] as char, bases[mid] as char),
        ));
    }

    contexts
}

fn is_gc(base: u8) -> bool {
    matches!(base, b'C' | b'G')
}

fn max_homopolymer_run(bases: &[u8]) -> usize {
    let mut best = 0;
    let mut current = 0;
    let mut prev = None;
    for &base in bases {
        if Some(base) == prev {
            current += 1;
        } else {
            current = 1;
            prev = Some(base);
        }
        best = best.max(current);
    }
    best
}

fn is_reverse_complement_palindrome(bases: &[u8]) -> bool {
    bases
        .iter()
        .zip(bases.iter().rev())
        .all(|(&lhs, &rhs)| complement_base(lhs) == rhs)
}

fn strand_idx(strand_rev: bool) -> usize {
    usize::from(strand_rev)
}

fn strand_label(strand_idx: usize) -> &'static str {
    match strand_idx {
        0 => "forward",
        1 => "reverse",
        _ => unreachable!(),
    }
}

fn read_oriented_base(base: u8, strand_rev: bool) -> u8 {
    if strand_rev {
        complement_base(base)
    } else {
        base
    }
}

fn substitution_idx(ref_base: u8, read_base: u8) -> Option<usize> {
    let ref_idx = base_idx(ref_base)?;
    let read_idx = base_idx(read_base)?;
    Some(ref_idx * 4 + read_idx)
}

fn substitution_bases(idx: usize) -> (u8, u8) {
    (BASES[idx / 4], BASES[idx % 4])
}

fn base_idx(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => base,
    }
}

fn ref_base(seq: &[u8], ref_pos: usize) -> u8 {
    seq[ref_pos - 1].to_ascii_uppercase()
}

fn dense_kmer_index(kmer_lens: &[usize]) -> FxHashMap<usize, Vec<usize>> {
    let mut index = FxHashMap::default();
    for &kmer_len in kmer_lens {
        if let Some(num_states) = dense_kmer_states(kmer_len) {
            index
                .entry(kmer_len)
                .or_insert_with(|| vec![INVALID_STAT_IDX; num_states]);
        }
    }
    index
}

fn dense_kmer_states(kmer_len: usize) -> Option<usize> {
    if kmer_len == 0 || kmer_len > MAX_ENCODED_KMER_LEN {
        return None;
    }

    let bits = 2 * kmer_len;
    if bits >= usize::BITS as usize {
        return None;
    }

    let num_states = 1usize << bits;
    if num_states <= MAX_DENSE_KMER_STATES {
        Some(num_states)
    } else {
        None
    }
}

fn encode_base(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn ref_base_opt(seq: &[u8], ref_pos: usize) -> Option<u8> {
    if ref_pos == 0 {
        None
    } else {
        seq.get(ref_pos - 1).map(|c| c.to_ascii_uppercase())
    }
}

fn decode_kmer(kmer_len: usize, mut bits: u64) -> String {
    let mut bases = vec![b'A'; kmer_len];
    for base in bases.iter_mut().rev() {
        *base = match bits & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        bits >>= 2;
    }

    String::from_utf8(bases).unwrap()
}
