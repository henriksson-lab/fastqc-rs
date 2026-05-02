#![allow(clippy::needless_range_loop)]

pub mod analysis;
pub mod charts;
pub mod cli;
pub mod config;
pub mod contaminant;
pub mod encoding;
#[cfg(feature = "gui")]
pub mod gui;
pub mod modules;
pub mod report;
pub mod sequence;
pub mod statistics;
pub mod utils;

use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};

use config::FastQCConfig;
use modules::adapter_content::AdapterContent;
use modules::basic_stats::BasicStats;
use modules::duplication_level::DuplicationLevel;
use modules::kmer_content::KmerContent;
use modules::module_config::ModuleConfig;
use modules::n_content::NContent;
use modules::overrepresented_seqs::{OverRepresentedSeqs, SharedDuplicationData};
use modules::per_base_quality::PerBaseQuality;
use modules::per_base_sequence_content::PerBaseSequenceContent;
use modules::per_sequence_gc_content::PerSequenceGCContent;
use modules::per_sequence_quality::PerSequenceQualityScores;
use modules::per_tile_quality::PerTileQuality;
use modules::sequence_length_distribution::SequenceLengthDistribution;
use modules::QCModule;
use sequence::bam_file;
use sequence::fast5_file::Fast5FileReader;
use sequence::fastq_file::FastQFile;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Create the standard module list, matching Java ModuleFactory.getStandardModuleList().
pub fn create_modules(config: &FastQCConfig) -> Vec<Box<dyn QCModule>> {
    let limits = config.limits_file.as_deref();
    let mc = || ModuleConfig::new(limits);

    let shared_data = Arc::new(RwLock::new(SharedDuplicationData::default()));

    let overrep = OverRepresentedSeqs::new(mc(), config.clone(), Arc::clone(&shared_data));
    let duplication = DuplicationLevel::new(mc(), Arc::clone(&shared_data));

    vec![
        Box::new(BasicStats::new(&mc())),
        Box::new(PerBaseQuality::new(mc(), config.clone())),
        Box::new(PerTileQuality::new(mc(), config.clone())),
        Box::new(PerSequenceQualityScores::new(mc())),
        Box::new(PerBaseSequenceContent::new(mc(), config.clone())),
        Box::new(PerSequenceGCContent::new(mc())),
        Box::new(NContent::new(mc(), config.clone())),
        Box::new(SequenceLengthDistribution::new(mc(), config)),
        Box::new(duplication),
        Box::new(overrep),
        Box::new(AdapterContent::new(mc(), config.clone())),
        Box::new(KmerContent::new(mc(), config.clone())),
    ]
}

/// Run FastQC analysis on a single file and return the data report text.
/// Detect input format from config or file extension.
/// Matches Java SequenceFactory logic.
pub fn detect_format(path: &Path, config: &FastQCConfig) -> String {
    if let Some(ref fmt) = config.sequence_format {
        return match fmt.as_str() {
            "bam" | "sam" | "bam_mapped" | "sam_mapped" => fmt.clone(),
            _ => "fastq".to_string(),
        };
    }
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_lowercase())
        .unwrap_or_default();
    if name.ends_with(".bam") || name.ends_with(".ubam") {
        "bam".to_string()
    } else if name.ends_with(".sam") {
        "sam".to_string()
    } else if name.ends_with(".fast5") {
        "fast5".to_string()
    } else {
        "fastq".to_string()
    }
}

/// Per-module QC result with structured data for programmatic access.
#[derive(Debug, Clone)]
pub struct ModuleResult {
    /// Module display name (e.g. "Basic Statistics")
    pub name: String,
    /// Pass / Warn / Fail
    pub status: modules::QCStatus,
    /// Tab-separated text data for this module (same as in fastqc_data.txt)
    pub data_text: String,
    /// Chart data for rendering, if this module produces a chart
    pub chart_data: Option<charts::ChartData>,
}

/// Result of running FastQC on a single file or sequence set.
#[derive(Debug)]
pub struct FastQCReport {
    /// The text data report (fastqc_data.txt content)
    pub data_report: String,
    /// The HTML report (fastqc_report.html content)
    pub html_report: String,
    /// FastQC `summary.txt` content.
    pub summary_report: String,
    /// Rendered chart images for a FastQC-compatible archive.
    pub chart_images: Vec<ChartImage>,
    /// Per-module structured results for programmatic access
    pub modules: Vec<ModuleResult>,
}

/// Rendered chart image stored in the FastQC archive.
#[derive(Debug, Clone)]
pub struct ChartImage {
    /// Archive filename, usually under `Images/` when written to a FastQC zip.
    pub filename: String,
    /// MIME type for the rendered chart, e.g. `image/png` or `image/svg+xml`.
    pub mime_type: String,
    /// Encoded image bytes.
    pub bytes: Vec<u8>,
}

#[derive(Debug, Clone)]
struct AnalysisInput {
    paths: Vec<PathBuf>,
    report_name: String,
    output_stem: String,
}

fn strip_fastqc_suffixes(name: &str) -> String {
    let mut stem = name.to_string();
    for suffix in [
        ".gz", ".bz2", ".txt", ".fastq", ".fq", ".csfastq", ".sam", ".bam", ".ubam",
    ] {
        if stem.to_lowercase().ends_with(suffix) {
            let new_len = stem.len() - suffix.len();
            stem.truncate(new_len);
        }
    }
    stem
}

fn file_name_string(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

fn casava_basename(name: &str) -> Option<String> {
    if name.ends_with(".fastq.gz") && name.len() >= 13 {
        let marker = name.len() - 13;
        let digits_start = name.len() - 12;
        let digits_end = name.len() - 9;
        if &name[marker..marker + 1] == "_"
            && name[digits_start..digits_end]
                .chars()
                .all(|c| c.is_ascii_digit())
        {
            return Some(format!("{}.fastq.gz", &name[..marker]));
        }
    } else if name.ends_with(".fastq") && name.len() >= 10 {
        let marker = name.len() - 10;
        let digits_start = name.len() - 9;
        let digits_end = name.len() - 6;
        if &name[marker..marker + 1] == "_"
            && name[digits_start..digits_end]
                .chars()
                .all(|c| c.is_ascii_digit())
        {
            return Some(format!("{}.fastq", &name[..marker]));
        }
    }

    None
}

fn nanopore_basename(name: &str) -> Option<String> {
    let without_ext = name.strip_suffix(".fast5").unwrap_or(name);
    let parts: Vec<&str> = without_ext.split('_').collect();
    if parts.len() < 3 {
        None
    } else {
        Some(format!("{}_{}_{}", parts[0], parts[1], parts[2]))
    }
}

fn collect_analysis_inputs(
    files: &[String],
    config: &FastQCConfig,
) -> Result<Vec<AnalysisInput>, Box<dyn std::error::Error>> {
    if files.len() == 1 && files[0].starts_with("stdin") {
        return Ok(vec![AnalysisInput {
            paths: Vec::new(),
            report_name: "stdin".to_string(),
            output_stem: "stdin".to_string(),
        }]);
    }

    let mut paths = Vec::new();
    for file in files {
        let path = PathBuf::from(file);
        if !path.exists() {
            eprintln!("Skipping {} - file not found", file);
            continue;
        }

        if config.nano && path.is_dir() {
            for entry in std::fs::read_dir(&path)? {
                let entry = entry?;
                let entry_path = entry.path();
                if entry_path.is_file() && file_name_string(&entry_path).ends_with(".fast5") {
                    paths.push(entry_path);
                } else if entry_path.is_dir() {
                    for sub_entry in std::fs::read_dir(&entry_path)? {
                        let sub_path = sub_entry?.path();
                        if sub_path.is_file() && file_name_string(&sub_path).ends_with(".fast5") {
                            paths.push(sub_path);
                        }
                    }
                }
            }
        } else {
            paths.push(path);
        }
    }

    if config.casava {
        let mut groups = std::collections::BTreeMap::<String, Vec<PathBuf>>::new();
        for path in paths {
            let name = file_name_string(&path);
            let key = if let Some(base) = casava_basename(&name) {
                base
            } else {
                eprintln!("File '{}' didn't look like part of a CASAVA group", name);
                name
            };
            groups.entry(key).or_default().push(path);
        }

        Ok(groups
            .into_iter()
            .map(|(report_name, paths)| AnalysisInput {
                output_stem: strip_fastqc_suffixes(&report_name),
                report_name,
                paths,
            })
            .collect())
    } else if config.nano {
        let mut groups = std::collections::BTreeMap::<String, Vec<PathBuf>>::new();
        for path in paths {
            let name = file_name_string(&path);
            if name.contains("muxscan") {
                continue;
            }
            let key = nanopore_basename(&name).unwrap_or_else(|| {
                eprintln!("File '{}' didn't look like part of a nanopore group", name);
                name.clone()
            });
            groups.entry(key).or_default().push(path);
        }

        Ok(groups
            .into_iter()
            .map(|(report_name, paths)| AnalysisInput {
                output_stem: strip_fastqc_suffixes(&report_name),
                report_name,
                paths,
            })
            .collect())
    } else {
        Ok(paths
            .into_iter()
            .map(|path| {
                let report_name = file_name_string(&path);
                AnalysisInput {
                    output_stem: strip_fastqc_suffixes(&report_name),
                    report_name,
                    paths: vec![path],
                }
            })
            .collect())
    }
}

/// Run FastQC analysis on a single file and return both reports.
pub fn run_fastqc_on_file(
    path: &Path,
    config: &FastQCConfig,
) -> Result<FastQCReport, Box<dyn std::error::Error>> {
    let file_name = file_name_string(path);
    run_fastqc_on_paths(&[path.to_path_buf()], &file_name, config)
}

fn run_fastqc_on_paths(
    paths: &[PathBuf],
    report_name: &str,
    config: &FastQCConfig,
) -> Result<FastQCReport, Box<dyn std::error::Error>> {
    if !config.quiet {
        eprintln!("Started analysis of {}", report_name);
    }

    let mut modules = create_modules(config);
    let mut count = 0u64;
    let sequence_report_name = if paths.len() > 1 {
        paths
            .first()
            .map(|path| file_name_string(path))
            .unwrap_or_else(|| report_name.to_string())
    } else {
        report_name.to_string()
    };

    for path in paths {
        let format = detect_format(path, config);
        let only_mapped = format == "bam_mapped" || format == "sam_mapped";
        let sequence_report_name = sequence_report_name.clone();

        count += match format.as_str() {
            "fast5" => {
                let fast5 = Fast5FileReader::open(path)?;
                let sequences = fast5.map(move |r| {
                    r.map(|mut seq| {
                        seq.file_name = sequence_report_name.clone();
                        seq
                    })
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
                });
                analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
            }
            "bam" | "bam_mapped" => {
                let bam = bam_file::BamFileReader::open(path, only_mapped)?;
                let sequences = bam.map(move |r| {
                    r.map(|mut seq| {
                        seq.file_name = sequence_report_name.clone();
                        seq
                    })
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
                });
                analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
            }
            "sam" | "sam_mapped" => {
                let sam = bam_file::SamFileReader::open(path, only_mapped)?;
                let sequences = sam.map(move |r| {
                    r.map(|mut seq| {
                        seq.file_name = sequence_report_name.clone();
                        seq
                    })
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
                });
                analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
            }
            _ => {
                let fq = FastQFile::open(path, config.casava, config.nofilter)?;
                let sequences = fq.map(move |r| {
                    r.map(|mut seq| {
                        seq.file_name = sequence_report_name.clone();
                        seq
                    })
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
                });
                analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
            }
        };
    }

    if !config.quiet {
        eprintln!(
            "Analysis complete for {} ({} sequences)",
            report_name, count
        );
    }

    let data_report = report::generate_data_report(&mut modules, VERSION);
    let summary_report = report::generate_summary_report(&mut modules, report_name);
    let html_report = report::generate_html_report(
        &mut modules,
        report_name,
        VERSION,
        config.svg_output,
        config.embed_images,
    );
    let chart_images = collect_chart_images(&mut modules, config.svg_output);

    // Extract structured per-module results
    let module_results: Vec<ModuleResult> = modules
        .iter_mut()
        .filter(|m| !m.ignore_in_report())
        .map(|m| {
            let mut data_text = String::new();
            m.make_data_report(&mut data_text);
            ModuleResult {
                name: m.name().to_string(),
                status: m.status(),
                data_text,
                chart_data: m.chart_data(),
            }
        })
        .collect();

    Ok(FastQCReport {
        data_report,
        html_report,
        summary_report,
        chart_images,
        modules: module_results,
    })
}

fn collect_chart_images(modules: &mut [Box<dyn QCModule>], svg_output: bool) -> Vec<ChartImage> {
    modules
        .iter_mut()
        .filter_map(|module| {
            let chart_data = module.chart_data()?;
            let mut filename = chart_data.image_filename(module.name());
            let (mime_type, bytes) = if svg_output {
                filename = filename
                    .strip_suffix(".png")
                    .map(|stem| format!("{}.svg", stem))
                    .unwrap_or_else(|| format!("{}.svg", filename));
                (
                    "image/svg+xml".to_string(),
                    charts::render_chart_to_svg(&chart_data).ok()?,
                )
            } else {
                (
                    "image/png".to_string(),
                    charts::render_chart_to_png(&chart_data).ok()?,
                )
            };

            Some(ChartImage {
                filename,
                mime_type,
                bytes,
            })
        })
        .collect()
}

/// Process a single FastQC input group: analyze, write reports.
fn process_analysis_input(
    input: &AnalysisInput,
    config: &FastQCConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let report = run_fastqc_on_paths(&input.paths, &input.report_name, config)?;
    let stem = &input.output_stem;

    let output_dir = config.output_dir.as_deref().unwrap_or_else(|| {
        input
            .paths
            .first()
            .and_then(|path| path.parent())
            .unwrap_or(Path::new("."))
    });

    // Write top-level HTML report
    let html_path = output_dir.join(format!("{}_fastqc.html", stem));
    std::fs::write(&html_path, &report.html_report)?;

    // Create ZIP archive containing the FastQC report folder
    let zip_path = output_dir.join(format!("{}_fastqc.zip", stem));
    report::write_fastqc_archive(
        &zip_path,
        stem,
        &report.data_report,
        &report.html_report,
        &report.summary_report,
        &report.chart_images,
    )?;

    // Handle --extract: unzip the archive
    if config.do_unzip {
        let extract_dir = output_dir.join(format!("{}_fastqc", stem));
        std::fs::create_dir_all(&extract_dir)?;
        std::fs::create_dir_all(extract_dir.join("Icons"))?;
        std::fs::create_dir_all(extract_dir.join("Images"))?;
        std::fs::write(extract_dir.join("fastqc_data.txt"), &report.data_report)?;
        std::fs::write(extract_dir.join("fastqc_report.html"), &report.html_report)?;
        std::fs::write(extract_dir.join("summary.txt"), &report.summary_report)?;
        for (name, bytes) in report::fastqc_icon_files() {
            std::fs::write(extract_dir.join("Icons").join(name), bytes)?;
        }
        for image in &report.chart_images {
            std::fs::write(
                extract_dir.join("Images").join(&image.filename),
                &image.bytes,
            )?;
        }

        // Handle --delete: remove ZIP after extraction
        if config.delete_after_unzip {
            std::fs::remove_file(&zip_path)?;
        }
    }

    if !config.quiet {
        eprintln!("Report written to {}", html_path.display());
    }

    Ok(())
}

/// Run FastQC analysis on the given files with the provided configuration.
/// Uses rayon thread pool when --threads > 1 for parallel file processing.
pub fn run_fastqc(
    files: &[String],
    config: &FastQCConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    if !config.quiet {
        eprintln!("fastqc-compliant-rs v{}", VERSION);
    }

    // Handle stdin specially when it is the only input, matching FastQC.
    if files.len() == 1 && files[0].starts_with("stdin") {
        // Process stdin sequentially (can't parallelize stdin)
        let fq = FastQFile::from_stdin(config.casava, config.nofilter)?;
        let mut modules = create_modules(config);
        let sequences = fq.map(|r| {
            r.map(|mut seq| {
                seq.file_name = "stdin".to_string();
                seq
            })
            .map_err(|e| Box::new(e) as Box<dyn std::error::Error>)
        });
        analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?;
        let data_report = report::generate_data_report(&mut modules, VERSION);
        let summary_report = report::generate_summary_report(&mut modules, "stdin");
        let html_report = report::generate_html_report(
            &mut modules,
            "stdin",
            VERSION,
            config.svg_output,
            config.embed_images,
        );
        let chart_images = collect_chart_images(&mut modules, config.svg_output);

        let output_dir = config.output_dir.as_deref().unwrap_or(Path::new("."));
        let html_path = output_dir.join("stdin_fastqc.html");
        std::fs::write(&html_path, &html_report)?;
        let zip_path = output_dir.join("stdin_fastqc.zip");
        report::write_fastqc_archive(
            &zip_path,
            "stdin",
            &data_report,
            &html_report,
            &summary_report,
            &chart_images,
        )?;
        if config.do_unzip {
            let extract_dir = output_dir.join("stdin_fastqc");
            std::fs::create_dir_all(&extract_dir)?;
            std::fs::create_dir_all(extract_dir.join("Icons"))?;
            std::fs::create_dir_all(extract_dir.join("Images"))?;
            std::fs::write(extract_dir.join("fastqc_data.txt"), &data_report)?;
            std::fs::write(extract_dir.join("fastqc_report.html"), &html_report)?;
            std::fs::write(extract_dir.join("summary.txt"), &summary_report)?;
            for (name, bytes) in report::fastqc_icon_files() {
                std::fs::write(extract_dir.join("Icons").join(name), bytes)?;
            }
            for image in &chart_images {
                std::fs::write(
                    extract_dir.join("Images").join(&image.filename),
                    &image.bytes,
                )?;
            }
            if config.delete_after_unzip {
                std::fs::remove_file(&zip_path)?;
            }
        }
        if !config.quiet {
            eprintln!("Report written to {}", html_path.display());
        }
        return Ok(());
    }

    let inputs = collect_analysis_inputs(files, config)?;

    if config.threads > 1 {
        // Parallel processing with rayon
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|e| format!("Failed to create thread pool: {}", e))?;

        let errors: Vec<String> = pool.install(|| {
            use rayon::prelude::*;
            inputs
                .par_iter()
                .filter_map(|input| {
                    if let Err(e) = process_analysis_input(input, config) {
                        Some(format!("{}: {}", input.report_name, e))
                    } else {
                        None
                    }
                })
                .collect()
        });

        if !errors.is_empty() {
            for err in &errors {
                eprintln!("Error: {}", err);
            }
            return Err(format!("{} file(s) failed to process", errors.len()).into());
        }
    } else {
        // Sequential processing
        for input in &inputs {
            process_analysis_input(input, config)?;
        }
    }

    Ok(())
}

/// Public API for running QC on in-memory sequences.
pub struct FastQCRunner {
    config: FastQCConfig,
}

impl FastQCRunner {
    pub fn new(config: FastQCConfig) -> Self {
        Self { config }
    }

    /// Run QC on an iterator of sequences and return structured results.
    pub fn run_sequences(
        &self,
        sequences: impl Iterator<Item = sequence::Sequence>,
    ) -> Result<FastQCReport, Box<dyn std::error::Error>> {
        let mut modules = create_modules(&self.config);
        let sequences = sequences.map(Ok::<_, Box<dyn std::error::Error>>);
        analysis::run_analysis(
            sequences,
            &mut modules,
            self.config.quiet,
            self.config.min_length,
        )?;
        let data_report = report::generate_data_report(&mut modules, VERSION);
        let summary_report = report::generate_summary_report(&mut modules, "sequences");
        let html_report = report::generate_html_report(
            &mut modules,
            "sequences",
            VERSION,
            self.config.svg_output,
            self.config.embed_images,
        );
        let chart_images = collect_chart_images(&mut modules, self.config.svg_output);

        let module_results: Vec<ModuleResult> = modules
            .iter_mut()
            .filter(|m| !m.ignore_in_report())
            .map(|m| {
                let mut data_text = String::new();
                m.make_data_report(&mut data_text);
                ModuleResult {
                    name: m.name().to_string(),
                    status: m.status(),
                    data_text,
                    chart_data: m.chart_data(),
                }
            })
            .collect();

        Ok(FastQCReport {
            data_report,
            html_report,
            summary_report,
            chart_images,
            modules: module_results,
        })
    }

    /// Run QC on a file path and return both reports.
    pub fn run_file(&self, path: &Path) -> Result<FastQCReport, Box<dyn std::error::Error>> {
        run_fastqc_on_file(path, &self.config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tempdir() -> PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir =
            std::env::temp_dir().join(format!("fastqc_rs_lib_{}_{}", std::process::id(), nanos));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn test_strip_fastqc_suffixes() {
        for (input, expected) in [
            ("reads.fastq", "reads"),
            ("reads.fq", "reads"),
            ("reads.fastq.gz", "reads"),
            ("reads.fq.gz", "reads"),
            ("reads.fastq.bz2", "reads"),
            ("reads.fq.bz2", "reads"),
            ("reads.txt", "reads"),
            ("reads.csfastq", "reads"),
            ("reads.sam", "reads"),
            ("reads.bam", "reads"),
            ("reads.ubam", "reads"),
        ] {
            assert_eq!(strip_fastqc_suffixes(input), expected);
        }
    }

    #[test]
    fn test_casava_basename_parsing() {
        assert_eq!(
            casava_basename("sample_001.fastq.gz"),
            Some("sample.fastq.gz".to_string())
        );
        assert_eq!(
            casava_basename("sample_123.fastq"),
            Some("sample.fastq".to_string())
        );
        assert_eq!(casava_basename("sample.fastq.gz"), None);
        assert_eq!(casava_basename("sample_ABC.fastq"), None);
    }

    #[test]
    fn test_nanopore_basename_parsing() {
        assert_eq!(
            nanopore_basename("Computer_Sample_42_ch100_file7_strand.fast5"),
            Some("Computer_Sample_42".to_string())
        );
        assert_eq!(
            nanopore_basename("Computer_Sample_42.fast5"),
            Some("Computer_Sample_42".to_string())
        );
        assert_eq!(nanopore_basename("short.fast5"), None);
    }

    #[test]
    fn test_nanopore_directory_scan_depth_and_muxscan_filter() {
        let dir = tempdir();
        let nested = dir.join("nested");
        let deep = nested.join("deep");
        std::fs::create_dir(&nested).unwrap();
        std::fs::create_dir(&deep).unwrap();

        std::fs::write(dir.join("Run_Sample_001_ch1_file1_strand.fast5"), b"").unwrap();
        std::fs::write(nested.join("Run_Sample_001_ch1_file2_strand.fast5"), b"").unwrap();
        std::fs::write(nested.join("Run_Sample_001_muxscan.fast5"), b"").unwrap();
        std::fs::write(deep.join("Run_Sample_001_ch1_file3_strand.fast5"), b"").unwrap();

        let config = FastQCConfig {
            nano: true,
            ..Default::default()
        };
        let inputs =
            collect_analysis_inputs(&[dir.to_string_lossy().to_string()], &config).unwrap();

        assert_eq!(inputs.len(), 1);
        assert_eq!(inputs[0].report_name, "Run_Sample_001");
        assert_eq!(inputs[0].paths.len(), 2);
        assert!(inputs[0]
            .paths
            .iter()
            .all(|path| !file_name_string(path).contains("muxscan")));
        assert!(inputs[0]
            .paths
            .iter()
            .all(|path| path.parent() != Some(deep.as_path())));
    }

    #[test]
    fn test_casava_invalid_names_become_singletons() {
        let dir = tempdir();
        let grouped_a = dir.join("sample_001.fastq");
        let grouped_b = dir.join("sample_002.fastq");
        let singleton = dir.join("other.fastq");
        for path in [&grouped_a, &grouped_b, &singleton] {
            std::fs::write(path, b"").unwrap();
        }

        let config = FastQCConfig {
            casava: true,
            ..Default::default()
        };
        let inputs = collect_analysis_inputs(
            &[
                grouped_a.to_string_lossy().to_string(),
                grouped_b.to_string_lossy().to_string(),
                singleton.to_string_lossy().to_string(),
            ],
            &config,
        )
        .unwrap();

        let grouped = inputs
            .iter()
            .find(|input| input.report_name == "sample.fastq")
            .unwrap();
        let single = inputs
            .iter()
            .find(|input| input.report_name == "other.fastq")
            .unwrap();

        assert_eq!(grouped.paths.len(), 2);
        assert_eq!(single.paths.len(), 1);
    }
}
