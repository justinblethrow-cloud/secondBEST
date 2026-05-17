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

use noodles::core::Position;
use noodles::fasta;

use crate::bed::FeatureInterval;

static COMPLEMENT: [u8; 128] = {
    let mut c = [0u8; 128];
    c[b'A' as usize] = b'T';
    c[b'T' as usize] = b'A';
    c[b'C' as usize] = b'G';
    c[b'G' as usize] = b'C';
    c[b'N' as usize] = b'N';
    c
};

/// Find homopolymers in a sequence to use as intervals.
pub fn find_homopolymers(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();
    let mut hp_len = 0;
    let mut prev = b'?';

    for i in start..=end {
        let curr = seq
            .get(Position::new(i).unwrap())
            .unwrap()
            .to_ascii_uppercase();

        if curr == prev && i != end {
            hp_len += 1;
            continue;
        }

        if hp_len > 1 {
            let c = if strand_rev {
                COMPLEMENT[prev as usize]
            } else {
                prev
            };
            res.push(FeatureInterval {
                start: i - hp_len,
                stop: i,
                val: format!("{: >5}{}", hp_len, c as char),
            });
        }
        hp_len = 1;
        prev = curr;
    }

    res
}

/// Get fixed-length windows as intervals.
pub fn get_windows(
    start: usize,
    end: usize,
    win_len: usize,
    pos: bool,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in (start..end).step_by(win_len) {
        let lo;
        let hi;
        if strand_rev {
            hi = end - (i - start);
            lo = hi.saturating_sub(win_len).max(start);
        } else {
            lo = i;
            hi = (i + win_len).min(end);
        };

        if pos {
            res.push(FeatureInterval {
                start: lo,
                stop: hi,
                val: format!("window_{}_pos_{}", win_len, i - start),
            });
        } else {
            res.push(FeatureInterval {
                start: lo,
                stop: hi,
                val: format!("window_{}", win_len),
            });
        }
    }

    res
}

const BORDER_CONTEXT: usize = 1;

/// Get small intervals that represent the region near fixed-width window borders.
pub fn get_borders(
    start: usize,
    end: usize,
    win_len: usize,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in (start..end).step_by(win_len).skip(1) {
        let idx;
        if strand_rev {
            idx = end - (i - start);
        } else {
            idx = i;
        };

        res.push(FeatureInterval {
            start: (idx - 1 - BORDER_CONTEXT).max(start),
            stop: (idx + BORDER_CONTEXT + 1).min(end),
            val: format!("border_{}", win_len),
        });
    }

    res
}

/// Get regions that match a sequence as intervals.
pub fn get_matches(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
    s: &str,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in start..end {
        // convert to zero-indexed
        let seq_iter = seq.as_ref()[i - 1..(i - 1 + s.len()).min(end - 1)]
            .iter()
            .map(|c| c.to_ascii_uppercase());
        let is_match = if strand_rev {
            seq_iter.eq(s.bytes().rev().map(|c| COMPLEMENT[c as usize]))
        } else {
            seq_iter.eq(s.bytes())
        };
        if is_match {
            res.push(FeatureInterval {
                start: i,
                stop: i + s.len(),
                val: s.to_owned(),
            });
        }
    }

    res
}

/// Get all A/C/G/T k-mers fully contained in a one-indexed [start, end) range.
///
/// The feature value is the k-mer sequence in read orientation. For alignments
/// on the reverse strand, this means the reverse complement of the reference
/// k-mer is reported.
pub fn get_kmers(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
    kmer_len: usize,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    if kmer_len == 0 || end <= start {
        return res;
    }

    let bounded_end = end.min(seq.len() + 1);
    if bounded_end <= start || bounded_end - start < kmer_len {
        return res;
    }

    let last_start = bounded_end - kmer_len;
    for i in start..=last_start {
        let kmer_start = i - 1;
        let kmer_end = kmer_start + kmer_len;
        let kmer = &seq.as_ref()[kmer_start..kmer_end];

        if !kmer.iter().all(|&c| is_acgt(c)) {
            continue;
        }

        let val = if strand_rev {
            reverse_complement(kmer)
        } else {
            uppercase_dna(kmer)
        };

        res.push(FeatureInterval {
            start: i,
            stop: i + kmer_len,
            val,
        });
    }

    res
}

fn is_acgt(c: u8) -> bool {
    matches!(c.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn uppercase_dna(seq: &[u8]) -> String {
    let bytes = seq.iter().map(|c| c.to_ascii_uppercase()).collect();
    String::from_utf8(bytes).unwrap()
}

fn reverse_complement(seq: &[u8]) -> String {
    let bytes = seq
        .iter()
        .rev()
        .map(|c| COMPLEMENT[c.to_ascii_uppercase() as usize])
        .collect();
    String::from_utf8(bytes).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(s: &str) -> fasta::record::Sequence {
        fasta::record::Sequence::from(s.as_bytes().to_vec())
    }

    fn interval_tuples(intervals: Vec<FeatureInterval>) -> Vec<(usize, usize, String)> {
        intervals
            .into_iter()
            .map(|i| (i.start, i.stop, i.val))
            .collect()
    }

    #[test]
    fn get_kmers_reports_forward_kmers_in_half_open_range() {
        let actual = interval_tuples(get_kmers(&seq("ACGTAC"), 2, 6, 3, false));

        assert_eq!(
            actual,
            vec![(2, 5, "CGT".to_string()), (3, 6, "GTA".to_string())]
        );
    }

    #[test]
    fn get_kmers_reports_reverse_complement_for_reverse_strand() {
        let actual = interval_tuples(get_kmers(&seq("ACGTAC"), 1, 5, 3, true));

        assert_eq!(
            actual,
            vec![(1, 4, "CGT".to_string()), (2, 5, "ACG".to_string())]
        );
    }

    #[test]
    fn get_kmers_skips_ambiguous_bases_and_uppercases_output() {
        let actual = interval_tuples(get_kmers(&seq("ACNTac"), 1, 7, 3, false));

        assert_eq!(actual, vec![(4, 7, "TAC".to_string())]);
    }

    #[test]
    fn get_kmers_returns_empty_for_zero_or_too_long_kmer_lengths() {
        assert!(get_kmers(&seq("ACGT"), 1, 5, 0, false).is_empty());
        assert!(get_kmers(&seq("ACGT"), 1, 3, 3, false).is_empty());
    }

    #[test]
    fn get_kmers_clamps_ranges_to_reference_length() {
        let actual = interval_tuples(get_kmers(&seq("ACGT"), 3, 20, 2, false));

        assert_eq!(actual, vec![(3, 5, "GT".to_string())]);
    }
}
