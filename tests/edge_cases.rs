use std::io::Write;
use std::path::PathBuf;

use fastqc_rs::config::FastQCConfig;
use fastqc_rs::sequence::Sequence;
use fastqc_rs::FastQCRunner;

// ── Empty / minimal input ────────────────────────────────────────────

#[test]
fn test_zero_sequences() {
    // No sequences at all — should not panic
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::empty()).unwrap();
    assert!(report.data_report.contains(">>Basic Statistics\tpass"));
    assert!(report.data_report.contains("Total Sequences\t0"));
}

#[test]
fn test_single_n_sequence() {
    // All-N sequence — should not panic, N content should be 100%
    let seq = Sequence::new(
        "test.fastq".into(),
        "NNNNNNNNNN".into(),
        "IIIIIIIIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::once(seq)).unwrap();
    assert!(report.data_report.contains("Per base N content"));
    // %GC should be 0 (no G or C)
    assert!(report.data_report.contains("%GC\t0"));
}

#[test]
fn test_min_length_filters_all() {
    // min_length > sequence length — should produce empty analysis, no panic
    let seq = Sequence::new(
        "test.fastq".into(),
        "ACGT".into(),
        "IIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        min_length: 100,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::once(seq)).unwrap();
    assert!(report.data_report.contains("Total Sequences\t0"));
}

#[test]
fn test_min_length_keeps_long_enough() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "AC".into(), "II".into(), "@short".into()),
        Sequence::new(
            "t.fq".into(),
            "ACGTACGTACGT".into(),
            "IIIIIIIIIIII".into(),
            "@long".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        min_length: 5,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(seqs.into_iter()).unwrap();
    // Only the 12bp read should be counted
    assert!(report.data_report.contains("Total Sequences\t1"));
    assert!(report.data_report.contains("Sequence length\t12"));
}

// ── FASTQ parsing edge cases ─────────────────────────────────────────

#[test]
fn test_gzip_fastq() {
    // Create a gzipped FASTQ in a temp file
    let dir = tempdir();
    let gz_path = dir.join("test.fastq.gz");
    {
        let file = std::fs::File::create(&gz_path).unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(b"@read1\nACGT\n+\nIIII\n").unwrap();
        gz.finish().unwrap();
    }

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_file(&gz_path).unwrap();
    assert!(report.data_report.contains("Total Sequences\t1"));
    assert!(report.data_report.contains("Sequence length\t4"));
}

#[test]
fn test_malformed_fastq_no_at() {
    let dir = tempdir();
    let path = dir.join("bad.fastq");
    std::fs::write(&path, "NOT_AN_AT_LINE\nACGT\n+\nIIII\n").unwrap();

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let result = runner.run_file(&path);
    assert!(result.is_err());
}

#[test]
fn test_malformed_fastq_seq_qual_mismatch() {
    let dir = tempdir();
    let path = dir.join("mismatch.fastq");
    std::fs::write(&path, "@read1\nACGTACGT\n+\nIII\n").unwrap();

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let result = runner.run_file(&path);
    assert!(result.is_err());
}

#[test]
fn test_truncated_fastq() {
    let dir = tempdir();
    let path = dir.join("truncated.fastq");
    std::fs::write(&path, "@read1\nACGT\n").unwrap(); // Missing + and quality

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let result = runner.run_file(&path);
    assert!(result.is_err());
}

// ── Module-specific tests ────────────────────────────────────────────

#[test]
fn test_gc_content_correct() {
    // 50% GC sequence
    let seq = Sequence::new(
        "t.fq".into(),
        "AACCGGTT".into(),
        "IIIIIIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::once(seq)).unwrap();
    assert!(report.data_report.contains("%GC\t50"));
}

#[test]
fn test_quality_encoding_sanger() {
    // Quality chars in Sanger range (ASCII 33-73)
    let seq = Sequence::new(
        "t.fq".into(),
        "ACGT".into(),
        "!!!!".into(), // ASCII 33 = lowest Sanger
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::once(seq)).unwrap();
    assert!(report.data_report.contains("Sanger / Illumina 1.9"));
}

#[test]
fn test_casava_filtering() {
    // Sequence with :Y: in header should be marked as filtered
    let seqs = vec![
        Sequence {
            id: "@read1:N:0:ACGT".into(),
            sequence: "ACGT".into(),
            quality: "IIII".into(),
            file_name: "t.fq".into(),
            colorspace: None,
            is_filtered: false,
        },
        Sequence {
            id: "@read2:Y:0:ACGT".into(),
            sequence: "ACGT".into(),
            quality: "IIII".into(),
            file_name: "t.fq".into(),
            colorspace: None,
            is_filtered: true, // Would be set by FASTQ parser in casava mode
        },
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(seqs.into_iter()).unwrap();
    // Both counted by BasicStats (ignoreFilteredSequences=false), but filtered one counted separately
    assert!(report.data_report.contains("Sequences flagged as poor quality\t1"));
    assert!(report.data_report.contains("Total Sequences\t1"));
}

#[test]
fn test_adapter_detection() {
    // PolyA sequence should be detected as adapter
    let seq = Sequence::new(
        "t.fq".into(),
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".into(), // 30bp polyA
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::once(seq)).unwrap();
    // Adapter Content section should show PolyA at 100%
    assert!(report.data_report.contains("Adapter Content"));
    assert!(report.data_report.contains("100.0"));
}

#[test]
fn test_variable_length_sequences() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "ACGT".into(), "IIII".into(), "@r1".into()),
        Sequence::new(
            "t.fq".into(),
            "ACGTACGTACGT".into(),
            "IIIIIIIIIIII".into(),
            "@r2".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(seqs.into_iter()).unwrap();
    // Should show range
    assert!(report.data_report.contains("Sequence length\t4-12"));
    // Should warn about non-uniform length
    assert!(report.data_report.contains("Sequence Length Distribution\twarn"));
}

#[test]
fn test_duplicate_sequences() {
    // 10 identical sequences should show up as duplicated
    let seqs: Vec<Sequence> = (0..10)
        .map(|i| {
            Sequence::new(
                "t.fq".into(),
                "ACGTACGTACGTACGT".into(),
                "IIIIIIIIIIIIIIII".into(),
                format!("@read{}", i),
            )
        })
        .collect();
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(seqs.into_iter()).unwrap();
    // Should show overrepresented
    assert!(report.data_report.contains("Overrepresented sequences"));
    assert!(report.data_report.contains("ACGTACGTACGTACGT"));
}

// ── format_java_double tests ─────────────────────────────────────────

#[test]
fn test_format_java_double_basics() {
    use fastqc_rs::modules::per_sequence_quality::format_java_double;
    assert_eq!(format_java_double(0.0), "0.0");
    assert_eq!(format_java_double(1.0), "1.0");
    assert_eq!(format_java_double(100.0), "100.0");
    assert_eq!(format_java_double(0.5), "0.5");
    assert_eq!(format_java_double(3.14), "3.14");
}

#[test]
fn test_format_java_double_scientific() {
    use fastqc_rs::modules::per_sequence_quality::format_java_double;
    // Small values should use scientific notation like Java
    let result = format_java_double(0.0005);
    assert!(result.contains("E"), "Expected scientific notation for 0.0005, got: {}", result);
    assert!(result.contains("5"), "Expected '5' in result for 0.0005, got: {}", result);
}

// ── Contaminant matching tests ───────────────────────────────────────

#[test]
fn test_contaminant_exact_match_short() {
    use fastqc_rs::contaminant::Contaminant;
    let c = Contaminant::new("TestAdapter", "AGATCGGAAGAG");
    // 12bp query matching the adapter exactly
    let hit = c.find_match("AGATCGGAAGAG");
    assert!(hit.is_some());
    let hit = hit.unwrap();
    assert_eq!(hit.length, 12);
    assert_eq!(hit.percent_id, 100);
}

#[test]
fn test_contaminant_no_match() {
    use fastqc_rs::contaminant::Contaminant;
    let c = Contaminant::new("TestAdapter", "AGATCGGAAGAG");
    // Completely different sequence
    let hit = c.find_match("TTTTTTTT");
    assert!(hit.is_none());
}

#[test]
fn test_contaminant_reverse_complement() {
    use fastqc_rs::contaminant::Contaminant;
    let c = Contaminant::new("TestAdapter", "AGATCGGAAGAG");
    // Reverse complement of AGATCGGAAGAG is CTCTTCCGATCT
    let hit = c.find_match("CTCTTCCGATCT");
    assert!(hit.is_some());
    let hit = hit.unwrap();
    assert_eq!(hit.percent_id, 100);
}

#[test]
fn test_contaminant_finder_known_adapter() {
    use fastqc_rs::contaminant::ContaminantFinder;
    let finder = ContaminantFinder::new(None); // Uses default contaminant list
    // Illumina Universal Adapter sequence
    let hit = finder.find_contaminant_hit("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    assert!(hit.is_some());
}

// ── HTML report tests ────────────────────────────────────────────────

#[test]
fn test_html_report_structure() {
    let seq = Sequence::new(
        "test.fastq".into(),
        "ACGTACGT".into(),
        "IIIIIIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_sequences(std::iter::once(seq)).unwrap();
    assert!(report.html_report.starts_with("<!DOCTYPE html>"));
    assert!(report.html_report.contains("<title>"));
    assert!(report.html_report.contains("Summary"));
    assert!(report.html_report.contains("Basic Statistics"));
    assert!(report.html_report.contains("</html>"));
}

// ── Helpers ──────────────────────────────────────────────────────────

fn tempdir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fastqc_rs_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}
