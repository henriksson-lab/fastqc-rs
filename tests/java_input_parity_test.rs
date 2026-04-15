use std::io::{Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::Command;

use fastqc_rs::config::FastQCConfig;
use fastqc_rs::{run_fastqc, run_fastqc_on_file};

type FastqRecord = (String, String, String);

fn tempdir() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "fastqc_rs_java_input_parity_{}_{}",
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

fn default_records() -> Vec<FastqRecord> {
    vec![
        ("read1".into(), "ACGTACGTACGT".into(), "IIIIIIIIIIII".into()),
        ("read2".into(), "GGGGCCCCAAAA".into(), "HHHHHHHHHHHH".into()),
        ("read3".into(), "TTTTAAAACCCC".into(), "FFFFFFFFFFFF".into()),
    ]
}

fn run_vendored_fastqc(base: &Path, inputs: &[PathBuf], outdir: &Path, properties: &[String]) {
    let mut command = Command::new("java");
    command
        .arg("-Xmx512m")
        .arg("-Dfastqc.quiet=true")
        .arg("-Dfastqc.unzip=false")
        .arg(format!("-Dfastqc.output_dir={}", outdir.display()))
        .arg("-Djava.awt.headless=true");
    for property in properties {
        command.arg(property);
    }
    let output = command
        .arg("-cp")
        .arg(vendored_fastqc_classpath(base))
        .arg("uk.ac.babraham.FastQC.FastQCApplication")
        .args(inputs)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "vendored FastQC failed for {:?}:\nstdout:\n{}\nstderr:\n{}",
        inputs,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn fastqc_data_from_outdir(outdir: &Path) -> String {
    let zips = std::fs::read_dir(outdir)
        .unwrap()
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.file_name()
                .unwrap()
                .to_string_lossy()
                .ends_with("_fastqc.zip")
        })
        .collect::<Vec<_>>();

    assert_eq!(
        zips.len(),
        1,
        "expected exactly one FastQC zip in {}, found {:?}",
        outdir.display(),
        zips
    );

    let zip_file = std::fs::File::open(&zips[0]).unwrap();
    let mut zip = zip::ZipArchive::new(zip_file).unwrap();
    let data_entry = (0..zip.len())
        .find(|&i| {
            let entry = zip.by_index(i).unwrap();
            entry.name().ends_with("/fastqc_data.txt")
        })
        .unwrap();

    let mut entry = zip.by_index(data_entry).unwrap();
    let mut contents = String::new();
    entry.read_to_string(&mut contents).unwrap();
    contents
}

fn normalize_report(report: &str) -> String {
    let mut normalized = Vec::new();
    let mut section = Vec::<String>::new();
    let mut section_name = String::new();

    for line in report.lines() {
        if line.starts_with("##FastQC") {
            continue;
        }

        if line.starts_with(">>") && line != ">>END_MODULE" {
            push_normalized_section(&mut normalized, &section_name, &mut section);
            section_name = line[2..].split('\t').next().unwrap().to_string();
            section.push(line.to_string());
        } else if line == ">>END_MODULE" {
            section.push(line.to_string());
            push_normalized_section(&mut normalized, &section_name, &mut section);
            section_name.clear();
        } else if !section.is_empty() {
            section.push(line.to_string());
        }
    }
    push_normalized_section(&mut normalized, &section_name, &mut section);

    normalized.join("\n")
}

fn push_normalized_section(
    normalized: &mut Vec<String>,
    section_name: &str,
    section: &mut Vec<String>,
) {
    if section.is_empty() {
        return;
    }

    if section_name == "Overrepresented sequences" && section.len() > 3 {
        let data_end = section.len() - 1;
        section[2..data_end].sort();
    }
    normalized.push(section.join("\n"));
    section.clear();
}

fn compare_reports(case_name: &str, java_report: &str, rust_report: &str) {
    assert_eq!(
        normalize_report(java_report),
        normalize_report(rust_report),
        "fastqc_data.txt differed for {}",
        case_name
    );
}

fn java_property(name: &str, value: impl std::fmt::Display) -> String {
    format!("-Dfastqc.{}={}", name, value)
}

fn run_java_single_case(
    base: &Path,
    case_dir: &Path,
    input: &Path,
    properties: &[String],
) -> String {
    let java_out = case_dir.join("java");
    std::fs::create_dir(&java_out).unwrap();
    run_vendored_fastqc(base, &[input.to_path_buf()], &java_out, properties);
    fastqc_data_from_outdir(&java_out)
}

#[test]
fn vendored_java_matches_compressed_fastq_inputs() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let cases = [
        ("gzip", dir.join("gzip.fastq.gz")),
        ("bzip2", dir.join("bzip2.fastq.bz2")),
    ];

    {
        let file = std::fs::File::create(&cases[0].1).unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        for (name, sequence, quality) in default_records() {
            writeln!(gz, "@{}\n{}\n+\n{}", name, sequence, quality).unwrap();
        }
        gz.finish().unwrap();
    }
    {
        let file = std::fs::File::create(&cases[1].1).unwrap();
        let mut bz = bzip2::write::BzEncoder::new(file, bzip2::Compression::default());
        for (name, sequence, quality) in default_records() {
            writeln!(bz, "@{}\n{}\n+\n{}", name, sequence, quality).unwrap();
        }
        bz.finish().unwrap();
    }

    for (case_name, input) in cases {
        let case_dir = dir.join(case_name);
        std::fs::create_dir(&case_dir).unwrap();
        let java_report = run_java_single_case(&base, &case_dir, &input, &[]);
        let rust_report = run_fastqc_on_file(
            &input,
            &FastQCConfig {
                quiet: true,
                ..Default::default()
            },
        )
        .unwrap()
        .data_report;
        compare_reports(case_name, &java_report, &rust_report);
    }
}

#[test]
fn vendored_java_matches_colorspace_fastq_input() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let input = dir.join("colorspace.csfastq");
    write_fastq(
        &input,
        &[
            ("cs1".into(), "T01230123".into(), "IIIIIIIII".into()),
            ("cs2".into(), "A12301230".into(), "HHHHHHHHH".into()),
            ("cs3".into(), "G23012301".into(), "FFFFFFFFF".into()),
        ],
    );

    let java_report = run_java_single_case(&base, &dir, &input, &[]);
    let rust_report = run_fastqc_on_file(
        &input,
        &FastQCConfig {
            quiet: true,
            ..Default::default()
        },
    )
    .unwrap()
    .data_report;
    compare_reports("colorspace", &java_report, &rust_report);
}

#[test]
fn vendored_java_matches_sam_and_mapped_sam_inputs() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let input = dir.join("reads.sam");
    std::fs::write(
        &input,
        concat!(
            "@SQ\tSN:chr1\tLN:100\n",
            "soft\t0\tchr1\t1\t60\t2S4M2S\t*\t0\t0\tTTACGTAA\tIIIIIIII\n",
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\tHHHH\n",
            "rev\t16\tchr1\t1\t60\t4M\t*\t0\t0\tAGTC\tABCD\n",
        ),
    )
    .unwrap();

    for (case_name, sequence_format) in [("sam", "sam"), ("sam_mapped", "sam_mapped")] {
        let case_dir = dir.join(case_name);
        std::fs::create_dir(&case_dir).unwrap();
        let properties = [java_property("sequence_format", sequence_format)];
        let java_report = run_java_single_case(&base, &case_dir, &input, &properties);
        let rust_report = run_fastqc_on_file(
            &input,
            &FastQCConfig {
                quiet: true,
                sequence_format: Some(sequence_format.to_string()),
                ..Default::default()
            },
        )
        .unwrap()
        .data_report;
        compare_reports(case_name, &java_report, &rust_report);
    }
}

#[test]
fn vendored_java_matches_bam_and_mapped_bam_inputs() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let input = dir.join("reads.bam");
    write_bam_fixture(&input);

    for (case_name, sequence_format) in [("bam", "bam"), ("bam_mapped", "bam_mapped")] {
        let case_dir = dir.join(case_name);
        std::fs::create_dir(&case_dir).unwrap();
        let properties = [java_property("sequence_format", sequence_format)];
        let java_report = run_java_single_case(&base, &case_dir, &input, &properties);
        let rust_report = run_fastqc_on_file(
            &input,
            &FastQCConfig {
                quiet: true,
                sequence_format: Some(sequence_format.to_string()),
                ..Default::default()
            },
        )
        .unwrap()
        .data_report;
        compare_reports(case_name, &java_report, &rust_report);
    }
}

#[test]
fn vendored_java_matches_casava_grouped_inputs() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let first = dir.join("sample_001.fastq");
    let second = dir.join("sample_002.fastq");
    write_fastq(
        &first,
        &[
            (
                "read1:N:0:ACGT".into(),
                "ACGTACGT".into(),
                "IIIIIIII".into(),
            ),
            (
                "read2:Y:0:ACGT".into(),
                "GGGGCCCC".into(),
                "IIIIIIII".into(),
            ),
        ],
    );
    write_fastq(
        &second,
        &[(
            "read3:N:0:ACGT".into(),
            "TTTTAAAA".into(),
            "HHHHHHHH".into(),
        )],
    );

    let java_out = dir.join("java");
    std::fs::create_dir(&java_out).unwrap();
    let properties = [java_property("casava", "true")];
    run_vendored_fastqc(
        &base,
        &[first.clone(), second.clone()],
        &java_out,
        &properties,
    );
    let java_report = fastqc_data_from_outdir(&java_out);

    let rust_out = dir.join("rust");
    std::fs::create_dir(&rust_out).unwrap();
    run_fastqc(
        &[
            first.to_string_lossy().to_string(),
            second.to_string_lossy().to_string(),
        ],
        &FastQCConfig {
            quiet: true,
            output_dir: Some(rust_out.clone()),
            casava: true,
            ..Default::default()
        },
    )
    .unwrap();
    let rust_report = fastqc_data_from_outdir(&rust_out);

    compare_reports("casava", &java_report, &rust_report);
}

#[test]
fn vendored_java_matches_min_length_and_dup_length_options() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let min_input = dir.join("min_length.fastq");
    write_fastq(
        &min_input,
        &[
            ("short".into(), "ACGT".into(), "IIII".into()),
            ("long".into(), "ACGTACGTACGT".into(), "IIIIIIIIIIII".into()),
        ],
    );
    let min_case = dir.join("min_length_case");
    std::fs::create_dir(&min_case).unwrap();
    let min_properties = [java_property("min_length", 5)];
    let java_report = run_java_single_case(&base, &min_case, &min_input, &min_properties);
    let rust_report = run_fastqc_on_file(
        &min_input,
        &FastQCConfig {
            quiet: true,
            min_length: 5,
            ..Default::default()
        },
    )
    .unwrap()
    .data_report;
    compare_reports("min_length", &java_report, &rust_report);

    let dup_input = dir.join("dup_length.fastq");
    write_fastq(
        &dup_input,
        &[
            (
                "same_prefix_1".into(),
                "ACGTACGTACGTAAAA".into(),
                "IIIIIIIIIIIIIIII".into(),
            ),
            (
                "same_prefix_2".into(),
                "ACGTACGTACGTCCCC".into(),
                "IIIIIIIIIIIIIIII".into(),
            ),
            (
                "different".into(),
                "TTTTGGGGAAAACCCC".into(),
                "IIIIIIIIIIIIIIII".into(),
            ),
        ],
    );
    let dup_case = dir.join("dup_length_case");
    std::fs::create_dir(&dup_case).unwrap();
    let dup_properties = [java_property("dup_length", 12)];
    let java_report = run_java_single_case(&base, &dup_case, &dup_input, &dup_properties);
    let rust_report = run_fastqc_on_file(
        &dup_input,
        &FastQCConfig {
            quiet: true,
            dup_length: 12,
            ..Default::default()
        },
    )
    .unwrap()
    .data_report;
    compare_reports("dup_length", &java_report, &rust_report);
}

#[test]
fn vendored_java_matches_custom_adapter_contaminant_and_limits_files() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !java_available() || !vendored_fastqc_available(&base) {
        eprintln!("skipping vendored Java FastQC input parity; java or FastQC classes unavailable");
        return;
    }

    let dir = tempdir();
    let input = dir.join("custom_config.fastq");
    write_fastq(
        &input,
        &[
            (
                "custom1".into(),
                "AACCGGTTAACC".into(),
                "IIIIIIIIIIII".into(),
            ),
            (
                "custom2".into(),
                "AACCGGTTAACC".into(),
                "IIIIIIIIIIII".into(),
            ),
            (
                "adapter".into(),
                "AAAACCCCGGGG".into(),
                "IIIIIIIIIIII".into(),
            ),
        ],
    );

    let adapter_file = dir.join("adapters.txt");
    let contaminant_file = dir.join("contaminants.txt");
    let limits_file = dir.join("limits.txt");
    std::fs::write(&adapter_file, "Custom Adapter\tAAAACCCC\n").unwrap();
    std::fs::write(&contaminant_file, "Custom Source\tAACCGGTTAACC\n").unwrap();
    std::fs::write(&limits_file, "adapter warn 0.0\nadapter error 100.0\n").unwrap();

    let properties = [
        java_property("adapter_file", adapter_file.display()),
        java_property("contaminant_file", contaminant_file.display()),
        java_property("limits_file", limits_file.display()),
    ];
    let java_report = run_java_single_case(&base, &dir, &input, &properties);
    let rust_report = run_fastqc_on_file(
        &input,
        &FastQCConfig {
            quiet: true,
            adapter_file: Some(adapter_file),
            contaminant_file: Some(contaminant_file),
            limits_file: Some(limits_file),
            ..Default::default()
        },
    )
    .unwrap()
    .data_report;

    compare_reports("custom_config", &java_report, &rust_report);
}

fn write_bam_fixture(path: &Path) {
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
            "HHHH",
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
