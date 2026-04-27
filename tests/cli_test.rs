use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastqc-compliant-rs")
}

fn tempdir() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("fastqc_rs_cli_{}_{}", std::process::id(), nanos));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_fastq(path: &Path, count: usize) {
    let mut fastq = String::new();
    for i in 0..count {
        fastq.push_str(&format!("@read{}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n", i));
    }
    std::fs::write(path, fastq).unwrap();
}

fn zip_names(path: &Path) -> Vec<String> {
    let zip_file = std::fs::File::open(path).unwrap();
    let zip = zip::ZipArchive::new(zip_file).unwrap();
    zip.file_names().map(|name| name.to_string()).collect()
}

fn zip_entry_string(path: &Path, entry_name: &str) -> String {
    let zip_file = std::fs::File::open(path).unwrap();
    let mut zip = zip::ZipArchive::new(zip_file).unwrap();
    let mut entry = zip.by_name(entry_name).unwrap();
    let mut contents = String::new();
    entry.read_to_string(&mut contents).unwrap();
    contents
}

#[test]
fn cli_default_writes_fastqc_outputs() {
    let dir = tempdir();
    let input = dir.join("sample.fastq");
    write_fastq(&input, 120);

    let output = Command::new(bin())
        .arg("--quiet")
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(dir.join("sample_fastqc.html").exists());
    assert!(dir.join("sample_fastqc.zip").exists());
    assert!(!dir.join("sample_fastqc_data.txt").exists());

    let data = zip_entry_string(
        &dir.join("sample_fastqc.zip"),
        "sample_fastqc/fastqc_data.txt",
    );
    assert!(data.contains("Total Sequences\t120"));
}

#[test]
fn cli_extract_delete_removes_zip_after_extraction() {
    let dir = tempdir();
    let input = dir.join("delete.fastq");
    write_fastq(&input, 120);

    let output = Command::new(bin())
        .arg("--quiet")
        .arg("--extract")
        .arg("--delete")
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(dir.join("delete_fastqc.html").exists());
    assert!(!dir.join("delete_fastqc.zip").exists());
    assert!(dir.join("delete_fastqc").join("fastqc_data.txt").exists());
    assert!(dir.join("delete_fastqc").join("summary.txt").exists());
}

#[test]
fn cli_svg_archives_svg_images() {
    let dir = tempdir();
    let input = dir.join("svg.fastq");
    write_fastq(&input, 120);

    let output = Command::new(bin())
        .arg("--quiet")
        .arg("--svg")
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let html = std::fs::read_to_string(dir.join("svg_fastqc.html")).unwrap();
    assert!(html.contains("src=\"Images/"));
    assert!(html.contains(".svg"));
    assert!(!html.contains("data:image/svg+xml"));

    let names = zip_names(&dir.join("svg_fastqc.zip"));
    assert!(names
        .iter()
        .any(|name| name.starts_with("svg_fastqc/Images/") && name.ends_with(".svg")));
}

#[test]
fn cli_embed_images_keeps_chart_data_uris() {
    let dir = tempdir();
    let input = dir.join("embed.fastq");
    write_fastq(&input, 120);

    let output = Command::new(bin())
        .arg("--quiet")
        .arg("--embed-images")
        .arg("--svg")
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let html = std::fs::read_to_string(dir.join("embed_fastqc.html")).unwrap();
    assert!(html.contains("data:image/png;base64,"));
    assert!(html.contains("data:image/svg+xml;base64,"));
    assert!(!html.contains("src=\"Images/"));
}

#[test]
fn cli_outdir_writes_outputs_to_existing_directory() {
    let dir = tempdir();
    let outdir = dir.join("results");
    std::fs::create_dir(&outdir).unwrap();
    let input = dir.join("outdir.fastq");
    write_fastq(&input, 2);

    let output = Command::new(bin())
        .arg("--quiet")
        .arg("--outdir")
        .arg(&outdir)
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(outdir.join("outdir_fastqc.html").exists());
    assert!(outdir.join("outdir_fastqc.zip").exists());
    assert!(!dir.join("outdir_fastqc.html").exists());
}

#[test]
fn cli_strips_fastqc_suffixes_for_output_names() {
    let cases = [
        ("reads.fastq", "reads_fastqc.html"),
        ("reads.fq", "reads_fastqc.html"),
        ("reads.fastq.gz", "reads_fastqc.html"),
        ("reads.fq.gz", "reads_fastqc.html"),
        ("reads.fastq.bz2", "reads_fastqc.html"),
        ("reads.fq.bz2", "reads_fastqc.html"),
        ("reads.txt", "reads_fastqc.html"),
        ("reads.csfastq", "reads_fastqc.html"),
        ("reads.sam", "reads_fastqc.html"),
        ("reads.bam", "reads_fastqc.html"),
        ("reads.ubam", "reads_fastqc.html"),
    ];

    for (input_name, expected_html) in cases {
        let dir = tempdir();
        let input = dir.join(input_name);
        if input_name.ends_with(".gz") {
            let file = std::fs::File::create(&input).unwrap();
            let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            gz.write_all(b"@read0\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
                .unwrap();
            gz.finish().unwrap();
        } else if input_name.ends_with(".bz2") {
            let file = std::fs::File::create(&input).unwrap();
            let mut bz = bzip2::write::BzEncoder::new(file, bzip2::Compression::default());
            bz.write_all(b"@read0\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
                .unwrap();
            bz.finish().unwrap();
        } else {
            write_fastq(&input, 1);
        }

        let output = Command::new(bin())
            .arg("--quiet")
            .arg("--format")
            .arg("fastq")
            .arg(&input)
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "case {} stderr: {}",
            input_name,
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            dir.join(expected_html).exists(),
            "case {} should write {}",
            input_name,
            expected_html
        );
    }
}

#[test]
fn cli_stdin_writes_stdin_outputs() {
    let dir = tempdir();
    let mut child = Command::new(bin())
        .current_dir(&dir)
        .arg("--quiet")
        .arg("stdin")
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();

    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(b"@read1\nACGTACGT\n+\nIIIIIIII\n")
        .unwrap();

    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(dir.join("stdin_fastqc.html").exists());
    assert!(dir.join("stdin_fastqc.zip").exists());

    let data = zip_entry_string(
        &dir.join("stdin_fastqc.zip"),
        "stdin_fastqc/fastqc_data.txt",
    );
    assert!(data.contains("Total Sequences\t1"));
}

#[test]
fn cli_fastqc_short_aliases_work() {
    let dir = tempdir();
    let outdir = dir.join("out");
    let tmpdir = dir.join("tmp");
    let adapters = dir.join("adapters.txt");
    let contaminants = dir.join("contaminants.txt");
    let limits = dir.join("limits.txt");
    let input = dir.join("aliases.fastq");

    std::fs::create_dir(&outdir).unwrap();
    std::fs::create_dir(&tmpdir).unwrap();
    std::fs::write(&adapters, "# none\n").unwrap();
    std::fs::write(&contaminants, "# none\n").unwrap();
    std::fs::write(&limits, include_str!("../src/resources/limits.txt")).unwrap();
    write_fastq(&input, 1);

    let output = Command::new(bin())
        .arg("-q")
        .arg("-o")
        .arg(&outdir)
        .arg("-t")
        .arg("1")
        .arg("-c")
        .arg(&contaminants)
        .arg("-a")
        .arg(&adapters)
        .arg("-l")
        .arg(&limits)
        .arg("-k")
        .arg("5")
        .arg("-f")
        .arg("fastq")
        .arg("-d")
        .arg(&tmpdir)
        .arg("-j")
        .arg("java")
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(outdir.join("aliases_fastqc.html").exists());
    assert!(outdir.join("aliases_fastqc.zip").exists());
}

#[test]
fn cli_accepts_kmer_range_boundaries() {
    for kmer_size in ["2", "10"] {
        let dir = tempdir();
        let input = dir.join(format!("kmers_{}.fastq", kmer_size));
        write_fastq(&input, 120);

        let output = Command::new(bin())
            .arg("--quiet")
            .arg("--kmers")
            .arg(kmer_size)
            .arg(&input)
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "kmer size {} stderr: {}",
            kmer_size,
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(dir.join(format!("kmers_{}_fastqc.zip", kmer_size)).exists());
    }
}

#[test]
fn cli_rejects_invalid_options() {
    let dir = tempdir();
    let input = dir.join("input.fastq");
    write_fastq(&input, 1);

    for (args, expected) in [
        (vec!["--kmers", "1"], "kmer size must be in the range 2-10"),
        (vec!["--kmers", "11"], "kmer size must be in the range 2-10"),
        (
            vec!["--threads", "0"],
            "number of threads must be a positive integer",
        ),
        (
            vec!["--nogroup", "--expgroup"],
            "you can't specify both --expgroup and --nogroup",
        ),
        (vec!["--format", "unknown"], "unrecognized format 'unknown'"),
    ] {
        let output = Command::new(bin())
            .args(args)
            .arg("--quiet")
            .arg(&input)
            .output()
            .unwrap();
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(expected),
            "expected {:?} in stderr: {}",
            expected,
            stderr
        );
    }
}

#[test]
fn cli_rejects_missing_paths() {
    let dir = tempdir();
    let input = dir.join("input.fastq");
    write_fastq(&input, 1);

    for (args, expected) in [
        (
            vec!["--outdir", "missing-results"],
            "specified output directory",
        ),
        (
            vec!["--contaminants", "missing-contaminants.txt"],
            "contaminant file",
        ),
        (vec!["--adapters", "missing-adapters.txt"], "adapter file"),
        (vec!["--limits", "missing-limits.txt"], "limits file"),
    ] {
        let output = Command::new(bin())
            .current_dir(&dir)
            .args(args)
            .arg("--quiet")
            .arg(&input)
            .output()
            .unwrap();
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(expected),
            "expected {:?} in stderr: {}",
            expected,
            stderr
        );
    }
}
