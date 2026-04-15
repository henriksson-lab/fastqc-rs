use std::io::{Read, Write};
use std::num::NonZero;
use std::path::PathBuf;

use fastqc_rs::config::FastQCConfig;
use fastqc_rs::sequence::bam_file::{BamFileReader, SamFileReader};
use fastqc_rs::sequence::Sequence;
use fastqc_rs::{run_fastqc, FastQCRunner};

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
fn test_n_content_no_n_all_n_and_position_spike() {
    let no_n = Sequence::new("t.fq".into(), "ACGT".into(), "IIII".into(), "@no_n".into());
    let all_n = Sequence::new("t.fq".into(), "NNNN".into(), "IIII".into(), "@all_n".into());
    let spike = Sequence::new("t.fq".into(), "ACNT".into(), "IIII".into(), "@spike".into());
    let config = FastQCConfig {
        quiet: true,
        nogroup: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences([no_n, all_n, spike].into_iter())
        .unwrap();
    let n_content = module_section(&report.data_report, "Per base N content");

    assert!(n_content.starts_with(">>Per base N content\tfail"));
    assert!(n_content.contains("1\t33.33333333333333"));
    assert!(n_content.contains("2\t33.33333333333333"));
    assert!(n_content.contains("3\t66.66666666666666"));
    assert!(n_content.contains("4\t33.33333333333333"));
}

#[test]
fn test_min_length_does_not_filter_reads() {
    // FastQC uses min_length to size base grouping, not to discard reads.
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
    assert!(report.data_report.contains("Total Sequences\t1"));
    assert!(report.data_report.contains("Sequence length\t4"));
    assert!(report.data_report.contains(">>Per base sequence quality"));
}

#[test]
fn test_sequence_length_distribution_all_filtered_input() {
    let seq = Sequence {
        id: "@filtered:Y:0:ACGT".into(),
        sequence: "ACGT".into(),
        quality: "IIII".into(),
        file_name: "t.fq".into(),
        colorspace: None,
        is_filtered: true,
    };
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();
    let length_section = module_section(&report.data_report, "Sequence Length Distribution");

    assert!(report.data_report.contains("Total Sequences\t0"));
    assert!(length_section.starts_with(">>Sequence Length Distribution\tpass"));
    assert!(length_section.contains("#Length\tCount"));
}

#[test]
fn test_sequence_length_distribution_zero_length_errors() {
    let seq = Sequence::new("t.fq".into(), "".into(), "".into(), "@empty".into());
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();
    let length_section = module_section(&report.data_report, "Sequence Length Distribution");

    assert!(length_section.starts_with(">>Sequence Length Distribution\tfail"));
    assert!(length_section.contains("0\t1.0"));
}

#[test]
fn test_min_length_preserves_short_and_long_reads() {
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
    assert!(report.data_report.contains("Total Sequences\t2"));
    assert!(report.data_report.contains("Sequence length\t2-12"));
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
fn test_multi_member_gzip_fastq() {
    let dir = tempdir();
    let gz_path = dir.join("multi.fastq.gz");
    {
        let file = std::fs::File::create(&gz_path).unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(b"@read1\nACGT\n+\nIIII\n").unwrap();
        gz.finish().unwrap();
    }
    {
        let file = std::fs::OpenOptions::new()
            .append(true)
            .open(&gz_path)
            .unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(b"@read2\nTGCA\n+\nHHHH\n").unwrap();
        gz.finish().unwrap();
    }

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config).run_file(&gz_path).unwrap();
    assert!(report.data_report.contains("Total Sequences\t2"));
    assert!(report.data_report.contains("Sequence length\t4"));
}

#[test]
fn test_bzip2_fastq() {
    let dir = tempdir();
    let bz2_path = dir.join("test.fastq.bz2");
    {
        let file = std::fs::File::create(&bz2_path).unwrap();
        let mut bz = bzip2::write::BzEncoder::new(file, bzip2::Compression::default());
        bz.write_all(b"@read1\nACGT\n+\nIIII\n").unwrap();
        bz.finish().unwrap();
    }

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let report = runner.run_file(&bz2_path).unwrap();
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
fn test_wrapped_fastq_is_rejected() {
    let dir = tempdir();
    let path = dir.join("wrapped.fastq");
    std::fs::write(&path, "@read1\nACGT\nACGT\n+\nIIII\nIIII\n").unwrap();

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let result = FastQCRunner::new(config).run_file(&path);
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
fn test_malformed_fastq_second_record_errors() {
    let dir = tempdir();
    let path = dir.join("bad_second.fastq");
    std::fs::write(&path, "@read1\nACGT\n+\nIIII\n@read2\nACGTACGT\n+\nIII\n").unwrap();

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let result = runner.run_file(&path);
    assert!(result.is_err());
}

#[test]
fn test_malformed_fastq_second_record_missing_plus_errors() {
    let dir = tempdir();
    let path = dir.join("bad_second_plus.fastq");
    std::fs::write(
        &path,
        "@read1\nACGT\n+\nIIII\n@read2\nACGT\nnot_plus\nIIII\n",
    )
    .unwrap();

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let runner = FastQCRunner::new(config);
    let result = runner.run_file(&path);
    assert!(result.is_err());
}

#[test]
fn test_truncated_fastq_second_record_errors() {
    let dir = tempdir();
    let path = dir.join("truncated_second.fastq");
    std::fs::write(&path, "@read1\nACGT\n+\nIIII\n@read2\nACGT\n").unwrap();

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

#[test]
fn test_sam_input_autodetect() {
    let dir = tempdir();
    let path = dir.join("reads.sam");
    std::fs::write(
        &path,
        "@SQ\tSN:chr1\tLN:100\nread1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
    )
    .unwrap();

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config).run_file(&path).unwrap();
    assert!(report.data_report.contains("Total Sequences\t1"));
    assert!(report.data_report.contains("Sequence length\t8"));
}

#[test]
fn test_sam_forced_format_on_non_sam_extension() {
    let dir = tempdir();
    let path = dir.join("reads.txt");
    std::fs::write(
        &path,
        "@SQ\tSN:chr1\tLN:100\nread1\t0\tchr1\t1\t60\t5M\t*\t0\t0\tACGTT\tIIIII\n",
    )
    .unwrap();

    let config = FastQCConfig {
        quiet: true,
        sequence_format: Some("sam".to_string()),
        ..Default::default()
    };
    let report = FastQCRunner::new(config).run_file(&path).unwrap();
    assert!(report.data_report.contains("Total Sequences\t1"));
}

#[test]
fn test_bam_input_autodetect_and_forced_format() {
    let dir = tempdir();
    let bam_path = dir.join("reads.bam");
    write_bam_fixture(&bam_path);

    let report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&bam_path)
    .unwrap();

    assert!(report.data_report.contains("Total Sequences\t3"));
    assert!(report.data_report.contains("Sequence length\t4-8"));

    let forced_path = dir.join("reads.dat");
    std::fs::copy(&bam_path, &forced_path).unwrap();
    let forced_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        sequence_format: Some("bam".to_string()),
        ..Default::default()
    })
    .run_file(&forced_path)
    .unwrap();

    assert!(forced_report.data_report.contains("Total Sequences\t3"));

    let mapped_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        sequence_format: Some("bam_mapped".to_string()),
        ..Default::default()
    })
    .run_file(&forced_path)
    .unwrap();

    assert!(mapped_report.data_report.contains("Total Sequences\t2"));
    assert!(mapped_report.data_report.contains("Sequence length\t4"));
}

#[test]
fn test_bam_mapped_reader_skips_unmapped_soft_clips_and_reverses() {
    let dir = tempdir();
    let path = dir.join("mapped.bam");
    write_bam_fixture(&path);

    let records = BamFileReader::open(&path, true)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "@soft");
    assert_eq!(records[0].sequence, "ACGT");
    assert_eq!(records[0].quality, "IIII");
    assert_eq!(records[1].id, "@rev");
    assert_eq!(records[1].sequence, "GACT");
    assert_eq!(records[1].quality, "DCBA");
}

#[test]
fn test_ubam_input_autodetects_as_bam() {
    let dir = tempdir();
    let path = dir.join("reads.ubam");
    write_bam_fixture(&path);

    let report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&path)
    .unwrap();

    assert!(report.data_report.contains("Total Sequences\t3"));
}

#[test]
fn test_sam_mapped_reader_skips_unmapped_and_soft_clips() {
    let dir = tempdir();
    let path = dir.join("mapped.sam");
    std::fs::write(
        &path,
        "@SQ\tSN:chr1\tLN:100\nsoft\t0\tchr1\t1\t60\t2S4M2S\t*\t0\t0\tTTACGTAA\tIIIIIIII\nunmapped\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n",
    )
    .unwrap();

    let mut reader = SamFileReader::open(&path, true).unwrap();
    let first = reader.next().unwrap().unwrap();
    assert_eq!(first.id, "@soft");
    assert_eq!(first.sequence, "ACGT");
    assert_eq!(first.quality, "IIII");
    assert!(reader.next().is_none());
}

#[test]
fn test_sam_reader_reverse_complements_reverse_records() {
    let dir = tempdir();
    let path = dir.join("reverse.sam");
    std::fs::write(
        &path,
        "@SQ\tSN:chr1\tLN:100\nrev\t16\tchr1\t1\t60\t4M\t*\t0\t0\tAGTC\tABCD\n",
    )
    .unwrap();

    let mut reader = SamFileReader::open(&path, false).unwrap();
    let first = reader.next().unwrap().unwrap();
    assert_eq!(first.id, "@rev");
    assert_eq!(first.sequence, "GACT");
    assert_eq!(first.quality, "DCBA");
}

#[test]
fn test_sam_missing_quality_does_not_panic() {
    let dir = tempdir();
    let path = dir.join("missing_quality.sam");
    std::fs::write(
        &path,
        "@SQ\tSN:chr1\tLN:100\nnoqual\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t*\n",
    )
    .unwrap();

    let config = FastQCConfig {
        quiet: true,
        sequence_format: Some("sam".to_string()),
        ..Default::default()
    };
    let report = FastQCRunner::new(config).run_file(&path).unwrap();
    assert!(report.data_report.contains("Total Sequences\t1"));
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
fn test_basic_stats_stdin_prefix_is_stripped() {
    let seq = Sequence::new(
        "stdin:stream.fastq".into(),
        "ACGT".into(),
        "IIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();
    let basic_stats = module_section(&report.data_report, "Basic Statistics");

    assert!(basic_stats.contains("Filename\tstream.fastq"));
}

#[test]
fn test_basic_stats_reports_colorspace_file_type() {
    let seq = Sequence::new_with_colorspace(
        "colors.csfastq".into(),
        "ACGT".into(),
        "A012".into(),
        "IIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();
    let basic_stats = module_section(&report.data_report, "Basic Statistics");

    assert!(basic_stats.contains("File type\tColorspace converted to bases"));
}

#[test]
fn test_basic_stats_filtered_only_encoding_falls_back_to_highest_char() {
    let seq = Sequence {
        id: "@filtered:Y:0:ACGT".into(),
        sequence: "ACGT".into(),
        quality: "!!!!".into(),
        file_name: "filtered.fq".into(),
        colorspace: None,
        is_filtered: true,
    };
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();
    let basic_stats = module_section(&report.data_report, "Basic Statistics");

    assert!(basic_stats.contains("Total Sequences\t0"));
    assert!(basic_stats.contains("Sequences flagged as poor quality\t1"));
    assert!(basic_stats.contains("Encoding\tIllumina 1.5"));
}

#[test]
fn test_per_base_quality_default_linear_grouping() {
    let seqs = repeated_quality_sequences(120, 100, '?');
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let quality_section = module_section(&report.data_report, "Per base sequence quality");

    assert!(quality_section.contains("\n1\t30.0\t30.0"));
    assert!(quality_section.contains("\n9\t30.0\t30.0"));
    assert!(quality_section.contains("\n10-11\t30.0\t30.0"));
}

#[test]
fn test_per_base_quality_nogroup_mode() {
    let seqs = repeated_quality_sequences(120, 100, '?');
    let config = FastQCConfig {
        quiet: true,
        nogroup: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let quality_section = module_section(&report.data_report, "Per base sequence quality");

    assert!(quality_section.contains("\n10\t30.0\t30.0"));
    assert!(quality_section.contains("\n100\t30.0\t30.0"));
    assert!(!quality_section.contains("\n10-11\t"));
}

#[test]
fn test_per_base_quality_expgroup_mode() {
    let seqs = repeated_quality_sequences(120, 100, '?');
    let config = FastQCConfig {
        quiet: true,
        expgroup: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let quality_section = module_section(&report.data_report, "Per base sequence quality");

    assert!(quality_section.contains("\n1\t30.0\t30.0"));
    assert!(quality_section.contains("\n9\t30.0\t30.0"));
    assert!(quality_section.contains("\n10-14\t30.0\t30.0"));
    assert!(quality_section.contains("\n95-99\t30.0\t30.0"));
}

#[test]
fn test_per_base_quality_warn_and_error_boundaries() {
    let warning_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(repeated_quality_sequences(120, 4, '5').into_iter())
    .unwrap();
    let warning_section = module_section(&warning_report.data_report, "Per base sequence quality");
    assert!(warning_section.starts_with(">>Per base sequence quality\twarn"));

    let failing_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(repeated_quality_sequences(120, 4, '4').into_iter())
    .unwrap();
    let failing_section = module_section(&failing_report.data_report, "Per base sequence quality");
    assert!(failing_section.starts_with(">>Per base sequence quality\tfail"));
}

#[test]
fn test_per_tile_quality_parses_illumina_tile_ids() {
    let seqs = vec![
        Sequence::new(
            "t.fq".into(),
            "ACGT".into(),
            "????".into(),
            "@INST:1:FCID:1:1101:100:100".into(),
        ),
        Sequence::new(
            "t.fq".into(),
            "ACGT".into(),
            "++++".into(),
            "@INST:1:FCID:1:1102:100:100".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        nogroup: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let tile_section = module_section(&report.data_report, "Per tile sequence quality");

    assert!(tile_section.starts_with(">>Per tile sequence quality\twarn"));
    assert!(tile_section.contains("#Tile\tBase\tMean"));
    assert!(tile_section.contains("1101\t1\t10.0"));
    assert!(tile_section.contains("1102\t1\t-10.0"));
}

#[test]
fn test_per_tile_quality_error_fixture() {
    let seqs = vec![
        Sequence::new(
            "t.fq".into(),
            "ACGT".into(),
            "????".into(),
            "@INST:1:FCID:1:1101:100:100".into(),
        ),
        Sequence::new(
            "t.fq".into(),
            "ACGT".into(),
            "!!!!".into(),
            "@INST:1:FCID:1:1102:100:100".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        nogroup: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let tile_section = module_section(&report.data_report, "Per tile sequence quality");

    assert!(tile_section.starts_with(">>Per tile sequence quality\tfail"));
    assert!(tile_section.contains("1101\t1\t15.0"));
    assert!(tile_section.contains("1102\t1\t-15.0"));
}

#[test]
fn test_per_tile_quality_non_tiled_headers_are_omitted() {
    let seq = Sequence::new(
        "t.fq".into(),
        "ACGT".into(),
        "????".into(),
        "@read_without_tile".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();

    assert!(!report.data_report.contains(">>Per tile sequence quality"));
}

#[test]
fn test_per_sequence_quality_low_and_mixed_distribution() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "ACGT".into(), "++++".into(), "@low1".into()),
        Sequence::new("t.fq".into(), "ACGT".into(), "++++".into(), "@low2".into()),
        Sequence::new("t.fq".into(), "ACGT".into(), "????".into(), "@high".into()),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let quality_section = module_section(&report.data_report, "Per sequence quality scores");

    assert!(quality_section.starts_with(">>Per sequence quality scores\tfail"));
    assert!(quality_section.contains("#Quality\tCount"));
    assert!(quality_section.contains("10\t2.0"));
    assert!(quality_section.contains("30\t1.0"));
}

#[test]
fn test_per_sequence_quality_histogram_includes_empty_bins() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "ACGT".into(), "!!!!".into(), "@q0".into()),
        Sequence::new("t.fq".into(), "ACGT".into(), "####".into(), "@q2".into()),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let quality_section = module_section(&report.data_report, "Per sequence quality scores");

    assert!(quality_section.contains("#Quality\tCount"));
    assert!(quality_section.contains("0\t1.0"));
    assert!(quality_section.contains("1\t0.0"));
    assert!(quality_section.contains("2\t1.0"));
}

#[test]
fn test_per_sequence_gc_content_all_n_all_gc_and_mixed_lengths() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "NNNN".into(), "????".into(), "@all_n".into()),
        Sequence::new(
            "t.fq".into(),
            "GGGG".into(),
            "????".into(),
            "@all_gc".into(),
        ),
        Sequence::new(
            "t.fq".into(),
            "ATGCATGC".into(),
            "????????".into(),
            "@mixed".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let gc_section = module_section(&report.data_report, "Per sequence GC content");

    assert!(gc_section.contains("#GC Content\tCount"));
    assert!(gc_section.contains("0\t1.0"));
    assert!(gc_section.contains("50\t1.0"));
    assert!(gc_section.contains("100\t1.0"));
}

#[test]
fn test_per_base_sequence_content_ignores_n_bases() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "AN".into(), "II".into(), "@read1".into()),
        Sequence::new("t.fq".into(), "CN".into(), "II".into(), "@read2".into()),
    ];
    let config = FastQCConfig {
        quiet: true,
        nogroup: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let sequence_content = module_section(&report.data_report, "Per base sequence content");

    assert!(sequence_content.contains("#Base\tG\tA\tT\tC"));
    assert!(sequence_content.contains("1\t0.0\t50.0\t0.0\t50.0"));
    assert!(sequence_content.contains("2\t0.0\t0.0\t0.0\t0.0"));
}

#[test]
fn test_per_base_sequence_content_grouped_and_ungrouped_positions() {
    let seq = Sequence::new(
        "t.fq".into(),
        "A".repeat(100),
        "?".repeat(100),
        "@read1".into(),
    );
    let default_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(std::iter::once(seq.clone()))
    .unwrap();
    let default_section = module_section(&default_report.data_report, "Per base sequence content");

    assert!(default_section.starts_with(">>Per base sequence content\tfail"));
    assert!(default_section.contains("\n1\t0.0\t100.0\t0.0\t0.0"));
    assert!(default_section.contains("\n10-11\t0.0\t100.0\t0.0\t0.0"));
    assert!(default_section.contains("\n100\t0.0\t100.0\t0.0\t0.0"));

    let ungrouped_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        nogroup: true,
        ..Default::default()
    })
    .run_sequences(std::iter::once(seq))
    .unwrap();
    let ungrouped_section =
        module_section(&ungrouped_report.data_report, "Per base sequence content");

    assert!(ungrouped_section.contains("\n1\t0.0\t100.0\t0.0\t0.0"));
    assert!(ungrouped_section.contains("\n100\t0.0\t100.0\t0.0\t0.0"));
    assert!(!ungrouped_section.contains("\n10-11\t"));
}

#[test]
fn test_per_base_sequence_content_warn_and_error_boundaries() {
    let pass_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(base_mix_sequences(10, 10).into_iter())
    .unwrap();
    let pass_section = module_section(&pass_report.data_report, "Per base sequence content");
    assert!(pass_section.starts_with(">>Per base sequence content\tpass"));

    let warning_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(base_mix_sequences(12, 8).into_iter())
    .unwrap();
    let warning_section = module_section(&warning_report.data_report, "Per base sequence content");
    assert!(warning_section.starts_with(">>Per base sequence content\twarn"));

    let failing_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(base_mix_sequences(13, 7).into_iter())
    .unwrap();
    let failing_section = module_section(&failing_report.data_report, "Per base sequence content");
    assert!(failing_section.starts_with(">>Per base sequence content\tfail"));
}

#[test]
fn test_adapter_svg_and_png_chart_modes_keep_same_data_report() {
    let dir = tempdir();
    let adapter_path = dir.join("custom_adapters.txt");
    std::fs::write(&adapter_path, "Custom Adapter\tAAA\n").unwrap();
    let seqs = repeated_sequences_with_base("t.fq", 120, 100, 'A', '?');

    let png_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        adapter_file: Some(adapter_path.clone()),
        ..Default::default()
    })
    .run_sequences(seqs.clone().into_iter())
    .unwrap();
    let svg_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        adapter_file: Some(adapter_path),
        svg_output: true,
        ..Default::default()
    })
    .run_sequences(seqs.into_iter())
    .unwrap();

    assert_eq!(png_report.data_report, svg_report.data_report);
    assert!(png_report.chart_images.iter().any(|image| {
        image.filename == "adapter_content.png" && image.mime_type == "image/png"
    }));
    assert!(svg_report.chart_images.iter().any(|image| {
        image.filename == "adapter_content.svg" && image.mime_type == "image/svg+xml"
    }));
}

#[test]
fn test_n_content_warn_and_error_boundaries() {
    let pass_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(n_spike_sequences(1, 20).into_iter())
    .unwrap();
    let pass_section = module_section(&pass_report.data_report, "Per base N content");
    assert!(pass_section.starts_with(">>Per base N content\tpass"));

    let warning_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(n_spike_sequences(4, 20).into_iter())
    .unwrap();
    let warning_section = module_section(&warning_report.data_report, "Per base N content");
    assert!(warning_section.starts_with(">>Per base N content\twarn"));

    let failing_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_sequences(n_spike_sequences(5, 20).into_iter())
    .unwrap();
    let failing_section = module_section(&failing_report.data_report, "Per base N content");
    assert!(failing_section.starts_with(">>Per base N content\tfail"));
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
    assert!(report
        .data_report
        .contains("Sequences flagged as poor quality\t1"));
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
fn test_custom_adapter_file_parsing() {
    let dir = tempdir();
    let adapter_path = dir.join("custom_adapters.txt");
    std::fs::write(&adapter_path, "# custom list\nCustom Adapter\tTTTTCCCC\n").unwrap();

    let seq = Sequence::new(
        "t.fq".into(),
        "AAAATTTTCCCC".into(),
        "IIIIIIIIIIII".into(),
        "@read1".into(),
    );
    let config = FastQCConfig {
        quiet: true,
        adapter_file: Some(adapter_path),
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(std::iter::once(seq))
        .unwrap();
    let adapter_section = module_section(&report.data_report, "Adapter Content");

    assert!(adapter_section.contains("#Position\tCustom Adapter"));
    assert!(adapter_section.contains("\t100.0"));
}

#[test]
fn test_adapter_cumulative_detection_by_position() {
    let dir = tempdir();
    let adapter_path = dir.join("custom_adapters.txt");
    std::fs::write(&adapter_path, "Custom Adapter\tAAA\n").unwrap();
    let seqs: Vec<Sequence> = (0..10)
        .map(|i| {
            Sequence::new(
                "t.fq".into(),
                "CCCCAAAACC".into(),
                "??????????".into(),
                format!("@read{}", i),
            )
        })
        .collect();
    let config = FastQCConfig {
        quiet: true,
        nogroup: true,
        adapter_file: Some(adapter_path),
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let adapter_section = module_section(&report.data_report, "Adapter Content");

    assert!(adapter_section.starts_with(">>Adapter Content\tfail"));
    assert!(adapter_section.contains("1\t0.0"));
    assert!(adapter_section.contains("4\t0.0"));
    assert!(adapter_section.contains("5\t100.0"));
    assert!(adapter_section.contains("8\t100.0"));
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
    assert!(report
        .data_report
        .contains("Sequence Length Distribution\twarn"));
}

#[test]
fn test_fixed_length_sequence_distribution_passes() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "ACGT".into(), "IIII".into(), "@r1".into()),
        Sequence::new("t.fq".into(), "TGCA".into(), "IIII".into(), "@r2".into()),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let length_section = module_section(&report.data_report, "Sequence Length Distribution");

    assert!(length_section.starts_with(">>Sequence Length Distribution\tpass"));
    assert!(length_section.contains("4\t2.0"));
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

#[test]
fn test_duplication_default_truncates_sequences_to_50_bases() {
    let shared_prefix = "A".repeat(50);
    let seqs = vec![
        Sequence::new(
            "t.fq".into(),
            format!("{}{}", shared_prefix, "C".repeat(10)),
            "I".repeat(60),
            "@read1".into(),
        ),
        Sequence::new(
            "t.fq".into(),
            format!("{}{}", shared_prefix, "G".repeat(10)),
            "I".repeat(60),
            "@read2".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let overrepresented_section = module_section(&report.data_report, "Overrepresented sequences");
    let duplication_section = module_section(&report.data_report, "Sequence Duplication Levels");

    assert!(overrepresented_section.contains(&format!("{}\t2\t100.0", shared_prefix)));
    assert!(!overrepresented_section.contains(&format!("{}{}", shared_prefix, "C")));
    assert!(duplication_section.contains("#Total Deduplicated Percentage\t50.0"));
    assert!(duplication_section.contains("2\t100.0"));
}

#[test]
fn test_dup_length_zero_uses_fastqc_50_base_fallback() {
    let shared_prefix = "T".repeat(50);
    let seqs = vec![
        Sequence::new(
            "t.fq".into(),
            format!("{}{}", shared_prefix, "A".repeat(10)),
            "I".repeat(60),
            "@read1".into(),
        ),
        Sequence::new(
            "t.fq".into(),
            format!("{}{}", shared_prefix, "C".repeat(10)),
            "I".repeat(60),
            "@read2".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        dup_length: 0,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let overrepresented_section = module_section(&report.data_report, "Overrepresented sequences");

    assert!(overrepresented_section.contains(&format!("{}\t2\t100.0", shared_prefix)));
}

#[test]
fn test_custom_dup_length_controls_truncated_sequence_key() {
    let seqs = vec![
        Sequence::new(
            "t.fq".into(),
            "ACGTACGTACGTAAAA".into(),
            "I".repeat(16),
            "@read1".into(),
        ),
        Sequence::new(
            "t.fq".into(),
            "ACGTACGTACGTCCCC".into(),
            "I".repeat(16),
            "@read2".into(),
        ),
    ];
    let config = FastQCConfig {
        quiet: true,
        dup_length: 12,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let overrepresented_section = module_section(&report.data_report, "Overrepresented sequences");
    let duplication_section = module_section(&report.data_report, "Sequence Duplication Levels");

    assert!(overrepresented_section.contains("ACGTACGTACGT\t2\t100.0"));
    assert!(!overrepresented_section.contains("ACGTACGTACGTAAAA"));
    assert!(duplication_section.contains("#Total Deduplicated Percentage\t50.0"));
}

#[test]
fn test_duplication_non_integer_dedup_percentage_formatting() {
    let seqs = vec![
        Sequence::new("t.fq".into(), "AAAA".into(), "????".into(), "@read1".into()),
        Sequence::new("t.fq".into(), "AAAA".into(), "????".into(), "@read2".into()),
        Sequence::new("t.fq".into(), "CCCC".into(), "????".into(), "@read3".into()),
    ];
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let duplication_section = module_section(&report.data_report, "Sequence Duplication Levels");

    assert!(duplication_section.contains("#Total Deduplicated Percentage\t66.66666666666666"));
    assert!(duplication_section.contains("1\t33.33333333333333"));
    assert!(duplication_section.contains("2\t66.66666666666666"));
}

#[test]
fn test_no_overrepresented_sequences_has_empty_section() {
    let seqs: Vec<Sequence> = (0..1024)
        .map(|i| {
            let sequence = dna_from_index(i, 10);
            Sequence::new(
                "t.fq".into(),
                sequence,
                "IIIIIIIIII".into(),
                format!("@read{}", i),
            )
        })
        .collect();
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let overrepresented_section = module_section(&report.data_report, "Overrepresented sequences");

    assert!(overrepresented_section.starts_with(">>Overrepresented sequences\tpass"));
    assert!(!overrepresented_section.contains("#Sequence\tCount\tPercentage\tPossible Source"));
}

// ── format_java_double tests ─────────────────────────────────────────

#[test]
fn test_format_java_double_basics() {
    use fastqc_rs::modules::per_sequence_quality::format_java_double;
    assert_eq!(format_java_double(0.0), "0.0");
    assert_eq!(format_java_double(1.0), "1.0");
    assert_eq!(format_java_double(100.0), "100.0");
    assert_eq!(format_java_double(0.5), "0.5");
    assert_eq!(format_java_double(3.125), "3.125");
}

#[test]
fn test_format_java_double_scientific() {
    use fastqc_rs::modules::per_sequence_quality::format_java_double;
    // Small values should use scientific notation like Java
    let result = format_java_double(0.0005);
    assert!(
        result.contains("E"),
        "Expected scientific notation for 0.0005, got: {}",
        result
    );
    assert!(
        result.contains("5"),
        "Expected '5' in result for 0.0005, got: {}",
        result
    );
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
    // Illumina Universal Adapter sequence
    let finder = ContaminantFinder::new(None); // Uses default contaminant list
    let hit = finder.find_contaminant_hit("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    assert!(hit.is_some());
}

#[test]
fn test_custom_contaminant_file_parsing_in_overrepresented_sequences() {
    let dir = tempdir();
    let contaminant_path = dir.join("custom_contaminants.txt");
    std::fs::write(
        &contaminant_path,
        "# custom list\nCustom Source\tAACCGGTTAACC\n",
    )
    .unwrap();

    let seqs: Vec<Sequence> = (0..20)
        .map(|i| {
            Sequence::new(
                "t.fq".into(),
                "AACCGGTTAACC".into(),
                "IIIIIIIIIIII".into(),
                format!("@read{}", i),
            )
        })
        .collect();
    let config = FastQCConfig {
        quiet: true,
        contaminant_file: Some(contaminant_path),
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let overrepresented_section = module_section(&report.data_report, "Overrepresented sequences");

    assert!(overrepresented_section.contains("AACCGGTTAACC\t20\t100.0"));
    assert!(overrepresented_section.contains("Custom Source"));
}

// ── K-mer tests ──────────────────────────────────────────────────────

#[test]
fn test_kmer_sequences_shorter_than_kmer_length() {
    let dir = tempdir();
    let limits_path = write_kmer_enabled_limits(&dir);
    let seqs: Vec<Sequence> = (0..100)
        .map(|i| {
            Sequence::new(
                "t.fq".into(),
                "ACG".into(),
                "III".into(),
                format!("@read{}", i),
            )
        })
        .collect();
    let config = FastQCConfig {
        quiet: true,
        kmer_size: 7,
        limits_file: Some(limits_path),
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let kmer_section = module_section(&report.data_report, "Kmer Content");

    assert!(kmer_section.starts_with(">>Kmer Content\tpass"));
    assert!(!kmer_section.contains("#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position"));
    assert!(!kmer_section.contains("ACG\t"));
}

#[test]
fn test_kmer_n_containing_kmers_are_not_reported() {
    let dir = tempdir();
    let limits_path = write_kmer_enabled_limits(&dir);
    let seqs: Vec<Sequence> = (0..100)
        .map(|i| {
            Sequence::new(
                "t.fq".into(),
                "NNNNNNNN".into(),
                "IIIIIIII".into(),
                format!("@read{}", i),
            )
        })
        .collect();
    let config = FastQCConfig {
        quiet: true,
        kmer_size: 3,
        limits_file: Some(limits_path),
        ..Default::default()
    };
    let report = FastQCRunner::new(config)
        .run_sequences(seqs.into_iter())
        .unwrap();
    let kmer_section = module_section(&report.data_report, "Kmer Content");

    assert!(kmer_section.starts_with(">>Kmer Content\tpass"));
    assert!(!kmer_section.contains("#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position"));
    assert!(!kmer_section.contains("NNN\t"));
}

// ── HTML report tests ────────────────────────────────────────────────

#[test]
fn test_html_report_structure() {
    let dir = tempdir();
    let input_path = dir.join("test.fastq");
    std::fs::write(&input_path, "@read1\nACGTACGT\n+\nIIIIIIII\n").unwrap();

    let report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();

    assert!(report.html_report.starts_with("<!DOCTYPE html>"));
    assert!(report
        .html_report
        .contains("<title>test.fastq FastQC Report</title>"));
    assert!(report.html_report.contains("Summary"));
    assert!(report.html_report.contains("Basic Statistics"));
    assert!(report
        .html_report
        .contains("<div id=\"header_filename\">\n    test.fastq"));
    assert!(report
        .html_report
        .contains("<a href=\"#M0\">Basic Statistics</a>"));
    assert!(report.html_report.contains("<h2 id=\"M0\"><img"));
    assert!(report.html_report.contains("alt=\"[PASS]\""));
    assert!(report.html_report.contains("<th>Measure</th>"));
    assert!(report.html_report.contains("<td>Filename</td>"));
    assert!(report.html_report.contains("<td>test.fastq</td>"));
    assert!(report.html_report.contains("</html>"));
}

#[test]
fn test_html_report_references_archive_images_by_default() {
    let dir = tempdir();
    let input_path = dir.join("html_images.fastq");
    write_repeated_fastq(&input_path, 120);

    let png_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();
    assert!(png_report
        .html_report
        .contains("src=\"Icons/fastqc_icon.png\""));
    assert!(png_report.html_report.contains("src=\"Icons/tick.png\""));
    assert!(png_report.html_report.contains("src=\"Images/"));
    assert!(png_report.html_report.contains(".png"));
    assert!(!png_report.html_report.contains("data:image/png;base64,"));

    let svg_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        svg_output: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();
    assert!(svg_report.html_report.contains("src=\"Images/"));
    assert!(svg_report.html_report.contains(".svg"));
    assert!(!svg_report
        .html_report
        .contains("data:image/svg+xml;base64,"));
}

#[test]
fn test_html_report_embed_images_option_keeps_data_uris() {
    let dir = tempdir();
    let input_path = dir.join("embedded_html_images.fastq");
    write_repeated_fastq(&input_path, 120);

    let report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        svg_output: true,
        embed_images: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();

    assert!(report.html_report.contains("data:image/png;base64,"));
    assert!(report.html_report.contains("data:image/svg+xml;base64,"));
    assert!(!report.html_report.contains("src=\"Icons/"));
    assert!(!report.html_report.contains("src=\"Images/"));
}

#[test]
fn test_html_report_escapes_file_and_table_text() {
    let dir = tempdir();
    let input_path = dir.join("weird<&\"'.fastq");
    write_repeated_fastq(&input_path, 1);

    let report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();

    assert!(report
        .html_report
        .contains("weird&lt;&amp;&quot;&#39;.fastq"));
    assert!(!report.html_report.contains("weird<&\"'.fastq"));
}

#[test]
fn test_fastqc_archive_layout() {
    let dir = tempdir();
    let input_path = dir.join("archive.fastq");
    write_repeated_fastq(&input_path, 120);

    let config = FastQCConfig {
        quiet: true,
        output_dir: Some(dir.clone()),
        ..Default::default()
    };

    run_fastqc(&[input_path.to_string_lossy().to_string()], &config).unwrap();

    assert!(dir.join("archive_fastqc.html").exists());
    assert!(dir.join("archive_fastqc.zip").exists());
    assert!(!dir.join("archive_fastqc_data.txt").exists());

    let names = zip_names(&dir.join("archive_fastqc.zip"));

    for expected in [
        "archive_fastqc/",
        "archive_fastqc/Icons/",
        "archive_fastqc/Images/",
        "archive_fastqc/Icons/fastqc_icon.png",
        "archive_fastqc/Icons/warning.png",
        "archive_fastqc/Icons/error.png",
        "archive_fastqc/Icons/tick.png",
        "archive_fastqc/summary.txt",
        "archive_fastqc/fastqc_data.txt",
        "archive_fastqc/fastqc_report.html",
    ] {
        assert!(
            names.iter().any(|name| name == expected),
            "missing {}",
            expected
        );
    }

    assert!(
        names
            .iter()
            .any(|name| name.starts_with("archive_fastqc/Images/")),
        "archive should include rendered chart images"
    );

    assert!(
        !names.iter().any(|name| name.ends_with("/fastqc.fo")),
        "fastqc.fo generation is out of scope"
    );
}

#[test]
fn test_public_archive_writer_helper() {
    let dir = tempdir();
    let input_path = dir.join("public_archive.fastq");
    write_repeated_fastq(&input_path, 120);

    let report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();

    let zip_path = dir.join("manual_fastqc.zip");
    fastqc_rs::report::write_fastqc_archive(
        &zip_path,
        "manual",
        &report.data_report,
        &report.html_report,
        &report.summary_report,
        &report.chart_images,
    )
    .unwrap();

    let names = zip_names(&zip_path);
    assert!(names.iter().any(|name| name == "manual_fastqc/"));
    assert!(names
        .iter()
        .any(|name| name == "manual_fastqc/fastqc_data.txt"));
    assert!(names.iter().any(|name| name == "manual_fastqc/summary.txt"));
    assert!(names
        .iter()
        .any(|name| name == "manual_fastqc/fastqc_report.html"));
    assert!(names
        .iter()
        .any(|name| name == "manual_fastqc/Icons/tick.png"));
    assert!(names
        .iter()
        .any(|name| name.starts_with("manual_fastqc/Images/")));
}

#[test]
fn test_casava_grouping_combines_files() {
    let dir = tempdir();
    let first = dir.join("sample_001.fastq");
    let second = dir.join("sample_002.fastq");
    std::fs::write(
        &first,
        "@read1:N:0:ACGT\nACGTACGT\n+\nIIIIIIII\n@read2:Y:0:ACGT\nACGTACGT\n+\nIIIIIIII\n",
    )
    .unwrap();
    std::fs::write(
        &second,
        "@read3:N:0:ACGT\nACGTACGT\n+\nIIIIIIII\n@read4:N:0:ACGT\nACGTACGT\n+\nIIIIIIII\n",
    )
    .unwrap();

    let config = FastQCConfig {
        quiet: true,
        output_dir: Some(dir.clone()),
        casava: true,
        ..Default::default()
    };

    run_fastqc(
        &[
            first.to_string_lossy().to_string(),
            second.to_string_lossy().to_string(),
        ],
        &config,
    )
    .unwrap();

    assert!(dir.join("sample_fastqc.html").exists());
    assert!(dir.join("sample_fastqc.zip").exists());
    assert!(!dir.join("sample_001_fastqc.zip").exists());
    assert!(!dir.join("sample_002_fastqc.zip").exists());

    let data = zip_entry_string(
        &dir.join("sample_fastqc.zip"),
        "sample_fastqc/fastqc_data.txt",
    );
    assert!(data.contains("Filename\tsample_001.fastq"));
    assert!(data.contains("Total Sequences\t3"));
    assert!(data.contains("Sequences flagged as poor quality\t1"));
}

#[test]
fn test_svg_output_embeds_and_archives_svg_images() {
    let dir = tempdir();
    let input_path = dir.join("svg.fastq");
    write_repeated_fastq(&input_path, 120);

    let png_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();
    let svg_report = FastQCRunner::new(FastQCConfig {
        quiet: true,
        svg_output: true,
        ..Default::default()
    })
    .run_file(&input_path)
    .unwrap();

    assert_eq!(png_report.data_report, svg_report.data_report);
    assert!(svg_report.html_report.contains("src=\"Images/"));
    assert!(svg_report.html_report.contains(".svg"));

    let config = FastQCConfig {
        quiet: true,
        output_dir: Some(dir.clone()),
        svg_output: true,
        ..Default::default()
    };
    run_fastqc(&[input_path.to_string_lossy().to_string()], &config).unwrap();

    let names = zip_names(&dir.join("svg_fastqc.zip"));
    assert!(
        names
            .iter()
            .any(|name| name.starts_with("svg_fastqc/Images/") && name.ends_with(".svg")),
        "archive should include SVG chart images"
    );
    assert!(
        !names
            .iter()
            .any(|name| name.starts_with("svg_fastqc/Images/") && name.ends_with(".png")),
        "SVG mode should not archive PNG chart images"
    );
}

#[test]
fn test_extract_and_delete_outputs() {
    let dir = tempdir();
    let input_path = dir.join("extract.fastq");
    write_repeated_fastq(&input_path, 120);

    let config = FastQCConfig {
        quiet: true,
        output_dir: Some(dir.clone()),
        do_unzip: true,
        ..Default::default()
    };
    run_fastqc(&[input_path.to_string_lossy().to_string()], &config).unwrap();

    let extract_dir = dir.join("extract_fastqc");
    assert!(dir.join("extract_fastqc.zip").exists());
    assert!(extract_dir.join("fastqc_data.txt").exists());
    assert!(extract_dir.join("fastqc_report.html").exists());
    assert!(extract_dir.join("summary.txt").exists());
    assert!(extract_dir.join("Icons").join("tick.png").exists());
    assert!(
        std::fs::read_dir(extract_dir.join("Images"))
            .unwrap()
            .next()
            .is_some(),
        "extracted Images directory should contain chart images"
    );

    let delete_dir = tempdir();
    let delete_input = delete_dir.join("delete.fastq");
    write_repeated_fastq(&delete_input, 120);
    let delete_config = FastQCConfig {
        quiet: true,
        output_dir: Some(delete_dir.clone()),
        do_unzip: true,
        delete_after_unzip: true,
        ..Default::default()
    };
    run_fastqc(
        &[delete_input.to_string_lossy().to_string()],
        &delete_config,
    )
    .unwrap();

    assert!(!delete_dir.join("delete_fastqc.zip").exists());
    assert!(delete_dir
        .join("delete_fastqc")
        .join("fastqc_data.txt")
        .exists());
}

#[test]
fn test_fast5_input_reports_unsupported() {
    let dir = tempdir();
    let input_path = dir.join("run_sample_001.fast5");
    std::fs::write(&input_path, b"not a real hdf5 file").unwrap();

    let config = FastQCConfig {
        quiet: true,
        nano: true,
        ..Default::default()
    };

    let err = match FastQCRunner::new(config).run_file(&input_path) {
        Ok(_) => panic!("Fast5 should fail with an explicit unsupported error"),
        Err(err) => err,
    };
    assert!(err
        .to_string()
        .contains("Fast5/Nanopore input is not implemented yet"));
}

// ── Helpers ──────────────────────────────────────────────────────────

fn write_repeated_fastq(path: &std::path::Path, count: usize) {
    let mut fastq = String::new();
    for i in 0..count {
        fastq.push_str(&format!("@read{}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n", i));
    }
    std::fs::write(path, fastq).unwrap();
}

fn write_bam_fixture(path: &std::path::Path) {
    use noodles::sam::{
        self,
        alignment::{
            record::{
                cigar::{op::Kind, Op},
                Flags,
            },
            record_buf::Cigar,
        },
        header::record::value::{map::ReferenceSequence, Map},
    };

    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(100).unwrap()),
        )
        .build();

    let file = std::fs::File::create(path).unwrap();
    let mut writer = noodles::bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    let soft_cigar: Cigar = [
        Op::new(Kind::SoftClip, 2),
        Op::new(Kind::Match, 4),
        Op::new(Kind::SoftClip, 2),
    ]
    .into_iter()
    .collect();
    let match_cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();

    for record in [
        bam_record("soft", Flags::empty(), soft_cigar, "TTACGTAA", "IIIIIIII"),
        bam_record(
            "unmapped",
            Flags::UNMAPPED,
            Cigar::default(),
            "GGGG",
            "IIII",
        ),
        bam_record(
            "rev",
            Flags::REVERSE_COMPLEMENTED,
            match_cigar,
            "AGTC",
            "ABCD",
        ),
    ] {
        noodles::sam::alignment::io::Write::write_alignment_record(&mut writer, &header, &record)
            .unwrap();
    }

    noodles::sam::alignment::io::Write::finish(&mut writer, &header).unwrap();
}

fn bam_record(
    name: &str,
    flags: noodles::sam::alignment::record::Flags,
    cigar: noodles::sam::alignment::record_buf::Cigar,
    sequence: &str,
    quality: &str,
) -> noodles::sam::alignment::RecordBuf {
    use noodles::core::Position;
    use noodles::sam::alignment::{
        record::MappingQuality,
        record_buf::{QualityScores, Sequence},
    };

    let quality_scores = QualityScores::from(quality.bytes().map(|b| b - b'!').collect::<Vec<_>>());
    let mut builder = noodles::sam::alignment::RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_sequence(Sequence::from(sequence.as_bytes()))
        .set_quality_scores(quality_scores);

    if !flags.is_unmapped() {
        builder = builder
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar);
    }

    builder.build()
}

fn repeated_quality_sequences(count: usize, len: usize, quality: char) -> Vec<Sequence> {
    let sequence = "A".repeat(len);
    let quality = quality.to_string().repeat(len);
    (0..count)
        .map(|i| {
            Sequence::new(
                "t.fq".into(),
                sequence.clone(),
                quality.clone(),
                format!("@read{}", i),
            )
        })
        .collect()
}

fn repeated_sequences_with_base(
    file_name: &str,
    count: usize,
    len: usize,
    base: char,
    quality: char,
) -> Vec<Sequence> {
    let sequence = base.to_string().repeat(len);
    let quality = quality.to_string().repeat(len);
    (0..count)
        .map(|i| {
            Sequence::new(
                file_name.to_string(),
                sequence.clone(),
                quality.clone(),
                format!("@read{}", i),
            )
        })
        .collect()
}

fn base_mix_sequences(a_count: usize, t_count: usize) -> Vec<Sequence> {
    let mut seqs = Vec::with_capacity(a_count + t_count);
    for i in 0..a_count {
        seqs.push(Sequence::new(
            "t.fq".into(),
            "A".into(),
            "?".into(),
            format!("@a{}", i),
        ));
    }
    for i in 0..t_count {
        seqs.push(Sequence::new(
            "t.fq".into(),
            "T".into(),
            "?".into(),
            format!("@t{}", i),
        ));
    }
    seqs
}

fn n_spike_sequences(n_count: usize, total_count: usize) -> Vec<Sequence> {
    (0..total_count)
        .map(|i| {
            let base = if i < n_count { "N" } else { "A" };
            Sequence::new(
                "t.fq".into(),
                base.into(),
                "?".into(),
                format!("@read{}", i),
            )
        })
        .collect()
}

fn zip_names(path: &std::path::Path) -> Vec<String> {
    let zip_file = std::fs::File::open(path).unwrap();
    let zip = zip::ZipArchive::new(zip_file).unwrap();
    zip.file_names().map(|name| name.to_string()).collect()
}

fn zip_entry_string(path: &std::path::Path, entry_name: &str) -> String {
    let zip_file = std::fs::File::open(path).unwrap();
    let mut zip = zip::ZipArchive::new(zip_file).unwrap();
    let mut entry = zip.by_name(entry_name).unwrap();
    let mut contents = String::new();
    entry.read_to_string(&mut contents).unwrap();
    contents
}

fn module_section<'a>(report: &'a str, module_name: &str) -> &'a str {
    let start_marker = format!(">>{}", module_name);
    let start = report
        .find(&start_marker)
        .unwrap_or_else(|| panic!("missing module section: {}", module_name));
    let rest = &report[start..];
    let end = rest
        .find(">>END_MODULE")
        .map(|i| i + ">>END_MODULE".len())
        .unwrap_or(rest.len());
    &rest[..end]
}

fn dna_from_index(mut index: usize, len: usize) -> String {
    let mut bases = vec!['A'; len];
    for pos in (0..len).rev() {
        bases[pos] = match index & 3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };
        index >>= 2;
    }
    bases.into_iter().collect()
}

fn write_kmer_enabled_limits(dir: &std::path::Path) -> PathBuf {
    let path = dir.join("limits.txt");
    std::fs::write(&path, "kmer ignore 0\n").unwrap();
    path
}

fn tempdir() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("fastqc_rs_test_{}_{}", std::process::id(), nanos));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}
