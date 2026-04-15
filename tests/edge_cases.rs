use std::io::{Read, Write};
use std::path::PathBuf;

use fastqc_rs::config::FastQCConfig;
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
    assert!(data.contains("Filename\tsample.fastq"));
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
    assert!(svg_report.html_report.contains("data:image/svg+xml"));

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

fn tempdir() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("fastqc_rs_test_{}_{}", std::process::id(), nanos));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}
