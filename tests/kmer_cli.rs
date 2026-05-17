use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use noodles::{bam, sam};
use sam::{
    alignment::Record,
    header::ReferenceSequence,
    record::{Flags, MappingQuality},
};

struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Result<Self, Box<dyn Error>> {
        let id = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
        let path =
            std::env::temp_dir().join(format!("best-{}-{}-{}", name, std::process::id(), id));
        fs::create_dir(&path)?;
        Ok(Self(path))
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn intervals_kmer_cli_writes_feature_summary() -> Result<(), Box<dyn Error>> {
    let tmp = TestDir::new("kmer-cli")?;
    let fasta_path = tmp.path().join("ref.fa");
    let bam_path = tmp.path().join("aln.bam");
    let output_prefix = tmp.path().join("out");

    write_fasta(&fasta_path)?;
    write_bam(&bam_path)?;

    let output = Command::new(env!("CARGO_BIN_EXE_best"))
        .arg("--no-per-aln-stats")
        .arg("--intervals-kmer")
        .arg("3")
        .arg("--")
        .arg(&bam_path)
        .arg(&fasta_path)
        .arg(&output_prefix)
        .output()?;

    assert!(
        output.status.success(),
        "best failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let summary_path = PathBuf::from(format!(
        "{}.summary_feature_stats.csv",
        output_prefix.display()
    ));
    let summary = fs::read_to_string(summary_path)?;
    let rows = parse_feature_summary(&summary);

    assert_eq!(rows.len(), 4);
    assert_eq!(rows["ACG"].intervals, 3);
    assert_eq!(rows["CGT"].intervals, 3);
    assert_eq!(rows["GTA"].intervals, 2);
    assert_eq!(rows["TAC"].intervals, 2);

    for feature in ["ACN", "CNA", "NAC"] {
        assert!(!rows.contains_key(feature));
    }

    for stats in rows.values() {
        assert_eq!(stats.identity, "1.000000");
        assert_eq!(stats.bases_per_interval, "3.000000");
        assert_eq!(stats.matches_per_interval, "3.000000");
    }

    let output_prefix_no_feature_qual = tmp.path().join("out_no_feature_qual");
    let output = Command::new(env!("CARGO_BIN_EXE_best"))
        .arg("--no-per-aln-stats")
        .arg("--no-feature-qual-score-stats")
        .arg("--record-batch-size")
        .arg("1")
        .arg("--intervals-kmer")
        .arg("3")
        .arg("--")
        .arg(&bam_path)
        .arg(&fasta_path)
        .arg(&output_prefix_no_feature_qual)
        .output()?;

    assert!(
        output.status.success(),
        "best failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let qual_summary_path = PathBuf::from(format!(
        "{}.summary_qual_score_stats.csv",
        output_prefix_no_feature_qual.display()
    ));
    let feature_summary_path = PathBuf::from(format!(
        "{}.summary_feature_stats.csv",
        output_prefix_no_feature_qual.display()
    ));
    let kmer_summary_path = PathBuf::from(format!(
        "{}.summary_kmer_stats.csv",
        output_prefix_no_feature_qual.display()
    ));
    let kmer_position_summary_path = PathBuf::from(format!(
        "{}.summary_kmer_position_stats.csv",
        output_prefix_no_feature_qual.display()
    ));
    let feature_summary = fs::read_to_string(feature_summary_path)?;
    assert_eq!(parse_feature_summary(&feature_summary), rows);

    let kmer_summary = fs::read_to_string(kmer_summary_path)?;
    let kmer_rows = parse_kmer_summary(&kmer_summary);
    assert_eq!(kmer_rows.len(), 4);
    assert_eq!(kmer_rows["ACG"].kmer_len, 3);
    assert_eq!(kmer_rows["ACG"].intervals, 3);
    assert_eq!(kmer_rows["CGT"].intervals, 3);
    assert_eq!(kmer_rows["GTA"].intervals, 2);
    assert_eq!(kmer_rows["TAC"].intervals, 2);

    assert!(!kmer_position_summary_path.exists());

    let output_prefix_kmer_position = tmp.path().join("out_kmer_position");
    let output = Command::new(env!("CARGO_BIN_EXE_best"))
        .arg("--no-per-aln-stats")
        .arg("--no-feature-qual-score-stats")
        .arg("--kmer-position-stats")
        .arg("--intervals-kmer")
        .arg("3")
        .arg("--")
        .arg(&bam_path)
        .arg(&fasta_path)
        .arg(&output_prefix_kmer_position)
        .output()?;

    assert!(
        output.status.success(),
        "best failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let kmer_position_summary_path = PathBuf::from(format!(
        "{}.summary_kmer_position_stats.csv",
        output_prefix_kmer_position.display()
    ));
    let kmer_position_summary = fs::read_to_string(kmer_position_summary_path)?;
    let kmer_position_rows = parse_kmer_position_summary(&kmer_position_summary);
    assert_eq!(kmer_position_rows.len(), 12);
    let acg_first_base = kmer_position_rows.get(&("ACG".to_string(), 0)).unwrap();
    assert_eq!(acg_first_base.kmer_len, 3);
    assert_eq!(acg_first_base.kmer_base, "A");
    assert_eq!(acg_first_base.intervals, 3);
    assert_eq!(acg_first_base.identity, "1.000000");
    assert_eq!(acg_first_base.matches, 3);
    assert_eq!(acg_first_base.mismatches, 0);

    let qual_summary = fs::read_to_string(qual_summary_path)?;
    let qual_rows = qual_summary.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(qual_rows, vec!["all_alignments,40,75.00"]);

    Ok(())
}

#[test]
fn intervals_kmer_cli_counts_mismatches_and_indels() -> Result<(), Box<dyn Error>> {
    let tmp = TestDir::new("kmer-errors-cli")?;
    let fasta_path = tmp.path().join("ref.fa");
    let bam_path = tmp.path().join("aln.bam");
    let output_prefix = tmp.path().join("out");

    write_nonperfect_fasta(&fasta_path)?;
    write_nonperfect_bam(&bam_path)?;

    let output = Command::new(env!("CARGO_BIN_EXE_best"))
        .arg("--no-per-aln-stats")
        .arg("--no-feature-qual-score-stats")
        .arg("--kmer-position-stats")
        .arg("--record-batch-size")
        .arg("1")
        .arg("--bam-reader-threads")
        .arg("1")
        .arg("--intervals-kmer")
        .arg("3")
        .arg("--")
        .arg(&bam_path)
        .arg(&fasta_path)
        .arg(&output_prefix)
        .output()?;

    assert!(
        output.status.success(),
        "best failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let kmer_summary_path = PathBuf::from(format!(
        "{}.summary_kmer_stats.csv",
        output_prefix.display()
    ));
    let kmer_summary = fs::read_to_string(kmer_summary_path)?;
    let kmer_rows = parse_kmer_summary(&kmer_summary);

    assert_eq!(kmer_rows.len(), 4);
    assert_eq!(kmer_rows["ACG"].intervals, 2);
    assert_eq!(kmer_rows["ACG"].identical_intervals, "0.000000");
    assert_eq!(kmer_rows["ACG"].identity, "0.666667");
    assert_eq!(kmer_rows["ACG"].matches_per_interval, "2.000000");
    assert_eq!(kmer_rows["ACG"].mismatches_per_interval, "0.500000");
    assert_eq!(kmer_rows["ACG"].non_hp_ins_per_interval, "0.000000");
    assert_eq!(kmer_rows["ACG"].non_hp_del_per_interval, "0.500000");

    assert_eq!(kmer_rows["CGT"].intervals, 2);
    assert_eq!(kmer_rows["CGT"].identity, "0.571429");
    assert_eq!(kmer_rows["CGT"].mismatches_per_interval, "0.500000");
    assert_eq!(kmer_rows["CGT"].non_hp_ins_per_interval, "0.500000");
    assert_eq!(kmer_rows["CGT"].non_hp_del_per_interval, "0.500000");

    assert_eq!(kmer_rows["GTA"].intervals, 2);
    assert_eq!(kmer_rows["GTA"].identity, "0.714286");
    assert_eq!(kmer_rows["GTA"].matches_per_interval, "2.500000");
    assert_eq!(kmer_rows["GTA"].non_hp_ins_per_interval, "0.500000");
    assert_eq!(kmer_rows["GTA"].non_hp_del_per_interval, "0.500000");

    assert_eq!(kmer_rows["TAC"].intervals, 2);
    assert_eq!(kmer_rows["TAC"].identical_intervals, "0.500000");
    assert_eq!(kmer_rows["TAC"].identity, "0.857143");
    assert_eq!(kmer_rows["TAC"].matches_per_interval, "3.000000");
    assert_eq!(kmer_rows["TAC"].non_hp_ins_per_interval, "0.500000");

    let position_summary_path = PathBuf::from(format!(
        "{}.summary_kmer_position_stats.csv",
        output_prefix.display()
    ));
    let position_summary = fs::read_to_string(position_summary_path)?;
    let position_rows = parse_kmer_position_summary(&position_summary);

    let acg_mismatch = position_rows.get(&("ACG".to_string(), 1)).unwrap();
    assert_eq!(acg_mismatch.kmer_base, "C");
    assert_eq!(acg_mismatch.matches, 1);
    assert_eq!(acg_mismatch.mismatches, 1);
    assert_eq!(acg_mismatch.non_hp_del, 0);

    let acg_deletion = position_rows.get(&("ACG".to_string(), 2)).unwrap();
    assert_eq!(acg_deletion.kmer_base, "G");
    assert_eq!(acg_deletion.matches, 1);
    assert_eq!(acg_deletion.mismatches, 0);
    assert_eq!(acg_deletion.non_hp_del, 1);

    let cgt_insertion = position_rows.get(&("CGT".to_string(), 2)).unwrap();
    assert_eq!(cgt_insertion.kmer_base, "T");
    assert_eq!(cgt_insertion.matches, 2);
    assert_eq!(cgt_insertion.non_hp_ins, 1);
    assert_eq!(cgt_insertion.non_hp_del, 0);

    Ok(())
}

fn write_fasta(path: &Path) -> Result<(), Box<dyn Error>> {
    fs::write(path, b">chr1\nACGTACNACG\n")?;
    Ok(())
}

fn write_nonperfect_fasta(path: &Path) -> Result<(), Box<dyn Error>> {
    fs::write(path, b">chr1\nACGTACGTAC\n")?;
    Ok(())
}

fn write_bam(path: &Path) -> Result<(), Box<dyn Error>> {
    let header = sam::Header::builder()
        .add_reference_sequence(ReferenceSequence::new("chr1".parse()?, 10)?)
        .build();

    let mut writer = bam::Writer::new(File::create(path)?);
    writer.write_header(&header)?;
    writer.write_reference_sequences(header.reference_sequences())?;

    for (name, flags) in [
        ("forward", Flags::empty()),
        ("reverse", Flags::REVERSE_COMPLEMENTED),
    ] {
        let record = Record::builder()
            .set_read_name(name.parse()?)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(1.try_into()?)
            .set_mapping_quality(MappingQuality::try_from(60)?)
            .set_cigar("10M".parse()?)
            .set_sequence("ACGTACNACG".parse()?)
            .set_quality_scores("IIIIIIIIII".parse()?)
            .build();

        writer.write_record(&header, &record)?;
    }

    writer.try_finish()?;
    Ok(())
}

fn write_nonperfect_bam(path: &Path) -> Result<(), Box<dyn Error>> {
    let header = sam::Header::builder()
        .add_reference_sequence(ReferenceSequence::new("chr1".parse()?, 10)?)
        .build();

    let mut writer = bam::Writer::new(File::create(path)?);
    writer.write_header(&header)?;
    writer.write_reference_sequences(header.reference_sequences())?;

    let record = Record::builder()
        .set_read_name("mismatch_ins_del".parse()?)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(1.try_into()?)
        .set_mapping_quality(MappingQuality::try_from(60)?)
        .set_cigar("3M1I3M1D3M".parse()?)
        .set_sequence("AAGGTACTAC".parse()?)
        .set_quality_scores("IIIIIIIIII".parse()?)
        .build();

    writer.write_record(&header, &record)?;
    writer.try_finish()?;
    Ok(())
}

#[derive(Debug, Eq, PartialEq)]
struct FeatureRow {
    intervals: usize,
    identity: String,
    bases_per_interval: String,
    matches_per_interval: String,
}

#[derive(Debug, Eq, PartialEq)]
struct KmerRow {
    kmer_len: usize,
    intervals: usize,
    identical_intervals: String,
    identity: String,
    matches_per_interval: String,
    mismatches_per_interval: String,
    non_hp_ins_per_interval: String,
    non_hp_del_per_interval: String,
}

#[derive(Debug, Eq, PartialEq)]
struct KmerPositionRow {
    kmer_len: usize,
    kmer_base: String,
    intervals: usize,
    identity: String,
    matches: usize,
    mismatches: usize,
    non_hp_ins: usize,
    non_hp_del: usize,
}

fn parse_feature_summary(summary: &str) -> HashMap<String, FeatureRow> {
    summary
        .lines()
        .skip(1)
        .map(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            assert_eq!(fields.len(), 13, "unexpected row: {}", line);

            (
                fields[0].to_string(),
                FeatureRow {
                    intervals: fields[1].parse().unwrap(),
                    identity: fields[3].to_string(),
                    bases_per_interval: fields[6].to_string(),
                    matches_per_interval: fields[7].to_string(),
                },
            )
        })
        .collect()
}

fn parse_kmer_summary(summary: &str) -> HashMap<String, KmerRow> {
    summary
        .lines()
        .skip(1)
        .map(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            assert_eq!(fields.len(), 14, "unexpected row: {}", line);

            (
                fields[1].to_string(),
                KmerRow {
                    kmer_len: fields[0].parse().unwrap(),
                    intervals: fields[2].parse().unwrap(),
                    identical_intervals: fields[3].to_string(),
                    identity: fields[4].to_string(),
                    matches_per_interval: fields[8].to_string(),
                    mismatches_per_interval: fields[9].to_string(),
                    non_hp_ins_per_interval: fields[10].to_string(),
                    non_hp_del_per_interval: fields[11].to_string(),
                },
            )
        })
        .collect()
}

fn parse_kmer_position_summary(summary: &str) -> HashMap<(String, usize), KmerPositionRow> {
    summary
        .lines()
        .skip(1)
        .map(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            assert_eq!(fields.len(), 15, "unexpected row: {}", line);

            (
                (fields[1].to_string(), fields[2].parse().unwrap()),
                KmerPositionRow {
                    kmer_len: fields[0].parse().unwrap(),
                    kmer_base: fields[3].to_string(),
                    intervals: fields[4].parse().unwrap(),
                    identity: fields[5].to_string(),
                    matches: fields[7].parse().unwrap(),
                    mismatches: fields[8].parse().unwrap(),
                    non_hp_ins: fields[9].parse().unwrap(),
                    non_hp_del: fields[10].parse().unwrap(),
                },
            )
        })
        .collect()
}
