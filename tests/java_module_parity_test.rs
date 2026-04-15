use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use fastqc_rs::config::FastQCConfig;
use fastqc_rs::run_fastqc_on_file;

const AUDITED_MODULES: &[&str] = &[
    "Per sequence quality scores",
    "Per sequence GC content",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Kmer Content",
];

type FastqRecord = (String, String, String);
type Fixture = (&'static str, Vec<FastqRecord>);

fn tempdir() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "fastqc_rs_java_parity_{}_{}",
        std::process::id(),
        nanos
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn java_available() -> bool {
    Command::new("java")
        .arg("-version")
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
}

fn vendored_fastqc_available(base: &Path) -> bool {
    base.join("fastqc/bin").is_dir()
        && base.join("fastqc/jbzip2-0.9.jar").is_file()
        && base.join("fastqc/htsjdk.jar").is_file()
        && base.join("fastqc/cisd-jhdf5.jar").is_file()
}

fn vendored_fastqc_classpath(base: &Path) -> String {
    [
        "fastqc/bin",
        "fastqc/jbzip2-0.9.jar",
        "fastqc/sam-1.103.jar",
        "fastqc/cisd-jhdf5.jar",
        "fastqc/htsjdk.jar",
        "fastqc/lib/*",
    ]
    .into_iter()
    .map(|part| base.join(part).to_string_lossy().to_string())
    .collect::<Vec<_>>()
    .join(":")
}

fn write_fastq(path: &Path, records: &[FastqRecord]) {
    let mut fastq = String::new();
    for (name, sequence, quality) in records {
        fastq.push('@');
        fastq.push_str(name);
        fastq.push('\n');
        fastq.push_str(sequence);
        fastq.push_str("\n+\n");
        fastq.push_str(quality);
        fastq.push('\n');
    }
    std::fs::write(path, fastq).unwrap();
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

fn parity_fixtures() -> Vec<Fixture> {
    let prefix = "A".repeat(50);
    let contaminant = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_string();

    let mut contaminant_records = Vec::new();
    for i in 0..20 {
        contaminant_records.push((
            format!("c{}", i),
            contaminant.clone(),
            "I".repeat(contaminant.len()),
        ));
    }
    for i in 0..20 {
        contaminant_records.push((format!("u{}", i), dna_from_index(i, 36), "I".repeat(36)));
    }

    let mut kmer_records = Vec::new();
    for i in 0..100 {
        kmer_records.push((format!("k{}", i), "ACGT".repeat(25), "I".repeat(100)));
    }

    vec![
        (
            "per_seq_quality",
            vec![
                ("q0".into(), "ACGT".into(), "!!!!".into()),
                ("q2".into(), "ACGT".into(), "####".into()),
                ("q30".into(), "ACGT".into(), "????".into()),
            ],
        ),
        (
            "per_seq_gc",
            vec![
                ("all_n".into(), "NNNN".into(), "!!!!".into()),
                ("all_gc".into(), "GGGG".into(), "!!!!".into()),
                ("mixed".into(), "ATGCATGC".into(), "!!!!!!!!".into()),
            ],
        ),
        (
            "dup_overrep",
            vec![
                (
                    "d1".into(),
                    format!("{}{}", prefix, "C".repeat(10)),
                    "I".repeat(60),
                ),
                (
                    "d2".into(),
                    format!("{}{}", prefix, "G".repeat(10)),
                    "I".repeat(60),
                ),
                ("u1".into(), "C".repeat(60), "I".repeat(60)),
            ],
        ),
        ("contaminant", contaminant_records),
        ("kmer", kmer_records),
    ]
}

fn run_vendored_fastqc(base: &Path, input: &Path, outdir: &Path) -> String {
    let output = Command::new("java")
        .arg("-Xmx512m")
        .arg("-Dfastqc.quiet=true")
        .arg("-Dfastqc.unzip=false")
        .arg(format!("-Dfastqc.output_dir={}", outdir.display()))
        .arg("-Djava.awt.headless=true")
        .arg("-cp")
        .arg(vendored_fastqc_classpath(base))
        .arg("uk.ac.babraham.FastQC.FastQCApplication")
        .arg(input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "vendored FastQC failed for {}:\nstdout:\n{}\nstderr:\n{}",
        input.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let stem = input.file_stem().unwrap().to_string_lossy();
    let zip_path = outdir.join(format!("{}_fastqc.zip", stem));
    let entry_name = format!("{}_fastqc/fastqc_data.txt", stem);
    let zip_file = std::fs::File::open(&zip_path).unwrap();
    let mut zip = zip::ZipArchive::new(zip_file).unwrap();
    let mut entry = zip.by_name(&entry_name).unwrap();
    let mut contents = String::new();
    entry.read_to_string(&mut contents).unwrap();
    contents
}

fn sections(report: &str) -> HashMap<String, String> {
    let mut sections = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_lines = Vec::new();

    for line in report.lines() {
        if line.starts_with("##FastQC") {
            continue;
        }

        if line.starts_with(">>") && line != ">>END_MODULE" {
            if let Some(name) = current_name.take() {
                sections.insert(name, current_lines.join("\n"));
            }
            current_name = Some(line[2..].split('\t').next().unwrap().to_string());
            current_lines.clear();
            current_lines.push(line.to_string());
        } else if line == ">>END_MODULE" {
            if let Some(name) = current_name.take() {
                current_lines.push(line.to_string());
                sections.insert(name, current_lines.join("\n"));
                current_lines.clear();
            }
        } else if current_name.is_some() {
            current_lines.push(line.to_string());
        }
    }

    sections
}

fn canonical_section(module: &str, section: &str) -> String {
    if module != "Overrepresented sequences" {
        return section.to_string();
    }

    let mut lines: Vec<String> = section.lines().map(|line| line.to_string()).collect();
    if lines.len() > 2 {
        lines[2..].sort();
    }
    lines.join("\n")
}

#[test]
fn vendored_java_module_outputs_match_controlled_fixtures() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC parity audit; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    for (fixture_name, records) in parity_fixtures() {
        let input = dir.join(format!("{}.fastq", fixture_name));
        write_fastq(&input, &records);

        let java_out = dir.join(format!("java_{}", fixture_name));
        std::fs::create_dir(&java_out).unwrap();
        let java_report = run_vendored_fastqc(&base, &input, &java_out);

        let rust_report = run_fastqc_on_file(
            &input,
            &FastQCConfig {
                quiet: true,
                ..Default::default()
            },
        )
        .unwrap()
        .data_report;

        let java_sections = sections(&java_report);
        let rust_sections = sections(&rust_report);

        for module in AUDITED_MODULES {
            let Some(java_section) = java_sections.get(*module) else {
                assert!(
                    !rust_sections.contains_key(*module),
                    "fastqc-rs reported {} for {}, but vendored FastQC omitted it",
                    module,
                    fixture_name
                );
                continue;
            };
            let rust_section = rust_sections
                .get(*module)
                .unwrap_or_else(|| panic!("fastqc-rs missing {} for {}", module, fixture_name));
            assert_eq!(
                canonical_section(module, java_section),
                canonical_section(module, rust_section),
                "module {} differed for fixture {}",
                module,
                fixture_name
            );
        }
    }
}
