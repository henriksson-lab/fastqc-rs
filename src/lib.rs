#![allow(clippy::needless_range_loop)]

pub mod analysis;
pub mod charts;
pub mod config;
pub mod contaminant;
#[cfg(feature = "gui")]
pub mod gui;
pub mod encoding;
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
pub struct FastQCReport {
    /// The text data report (fastqc_data.txt content)
    pub data_report: String,
    /// The HTML report with embedded charts (fastqc_report.html content)
    pub html_report: String,
    /// Per-module structured results for programmatic access
    pub modules: Vec<ModuleResult>,
}

/// Run FastQC analysis on a single file and return both reports.
pub fn run_fastqc_on_file(
    path: &Path,
    config: &FastQCConfig,
) -> Result<FastQCReport, Box<dyn std::error::Error>> {
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    if !config.quiet {
        eprintln!("Started analysis of {}", file_name);
    }

    let mut modules = create_modules(config);

    let format = detect_format(path, config);
    let only_mapped = format == "bam_mapped" || format == "sam_mapped";

    let count = match format.as_str() {
        "bam" | "bam_mapped" => {
            let bam = bam_file::BamFileReader::open(path, only_mapped)?;
            let sequences = bam.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
            analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
        }
        "sam" | "sam_mapped" => {
            let sam = bam_file::SamFileReader::open(path, only_mapped)?;
            let sequences = sam.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
            analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
        }
        _ => {
            let fq = FastQFile::open(path, config.casava, config.nofilter)?;
            let sequences = fq.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
            analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?
        }
    };

    if !config.quiet {
        eprintln!(
            "Analysis complete for {} ({} sequences)",
            file_name, count
        );
    }

    let data_report = report::generate_data_report(&mut modules, VERSION);
    let html_report = report::generate_html_report(&mut modules, &file_name, VERSION);

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
        modules: module_results,
    })
}

/// Create a ZIP archive containing both the data and HTML reports.
/// The ZIP structure mirrors Java FastQC: <stem>_fastqc/fastqc_data.txt and <stem>_fastqc/fastqc_report.html
fn create_zip_archive(
    zip_path: &Path,
    stem: &str,
    data_report: &str,
    html_report: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::io::Write;
    use zip::write::SimpleFileOptions;
    use zip::ZipWriter;

    let file = std::fs::File::create(zip_path)?;
    let mut zip = ZipWriter::new(file);
    let options = SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);

    let prefix = format!("{}_fastqc", stem);

    // Write fastqc_data.txt
    zip.start_file(format!("{}/fastqc_data.txt", prefix), options)?;
    zip.write_all(data_report.as_bytes())?;

    // Write fastqc_report.html
    zip.start_file(format!("{}/fastqc_report.html", prefix), options)?;
    zip.write_all(html_report.as_bytes())?;

    zip.finish()?;
    Ok(())
}

/// Process a single file: analyze, write reports.
fn process_single_file(
    path: &Path,
    config: &FastQCConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let report = run_fastqc_on_file(path, config)?;

    // Determine output path
    let stem = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "output".to_string());
    // Strip .fastq, .fq, etc. from stem
    let stem = stem
        .strip_suffix(".fastq")
        .or_else(|| stem.strip_suffix(".fq"))
        .unwrap_or(&stem)
        .to_string();

    let output_dir = config
        .output_dir
        .as_deref()
        .unwrap_or_else(|| path.parent().unwrap_or(Path::new(".")));

    // Write text data report
    let data_path = output_dir.join(format!("{}_fastqc_data.txt", stem));
    std::fs::write(&data_path, &report.data_report)?;

    // Write HTML report
    let html_path = output_dir.join(format!("{}_fastqc_report.html", stem));
    std::fs::write(&html_path, &report.html_report)?;

    // Create ZIP archive containing both reports
    let zip_path = output_dir.join(format!("{}_fastqc.zip", stem));
    create_zip_archive(&zip_path, &stem, &report.data_report, &report.html_report)?;

    // Handle --extract: unzip the archive
    if config.do_unzip {
        let extract_dir = output_dir.join(format!("{}_fastqc", stem));
        std::fs::create_dir_all(&extract_dir)?;
        std::fs::write(extract_dir.join("fastqc_data.txt"), &report.data_report)?;
        std::fs::write(
            extract_dir.join("fastqc_report.html"),
            &report.html_report,
        )?;

        // Handle --delete: remove ZIP after extraction
        if config.delete_after_unzip {
            std::fs::remove_file(&zip_path)?;
        }
    }

    if !config.quiet {
        eprintln!("Report written to {}", data_path.display());
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
        eprintln!("FastQC-rs v{}", VERSION);
    }

    // Handle stdin specially
    let stdin_files: Vec<&String> = files.iter().filter(|f| f.starts_with("stdin")).collect();
    if !stdin_files.is_empty() {
        // Process stdin sequentially (can't parallelize stdin)
        let fq = FastQFile::from_stdin(config.casava, config.nofilter)?;
        let mut modules = create_modules(config);
        let sequences = fq.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
        analysis::run_analysis(sequences, &mut modules, config.quiet, config.min_length)?;
        let data_report = report::generate_data_report(&mut modules, VERSION);
        let html_report = report::generate_html_report(&mut modules, "stdin", VERSION);

        let output_dir = config.output_dir.as_deref().unwrap_or(Path::new("."));
        let data_path = output_dir.join("stdin_fastqc_data.txt");
        std::fs::write(&data_path, &data_report)?;
        let html_path = output_dir.join("stdin_fastqc_report.html");
        std::fs::write(&html_path, &html_report)?;
        let zip_path = output_dir.join("stdin_fastqc.zip");
        create_zip_archive(&zip_path, "stdin", &data_report, &html_report)?;
        if !config.quiet {
            eprintln!("Report written to {}", data_path.display());
        }
        return Ok(());
    }

    // Filter to existing files
    let paths: Vec<PathBuf> = files
        .iter()
        .filter_map(|file| {
            let path = PathBuf::from(file);
            if path.exists() {
                Some(path)
            } else {
                eprintln!("Skipping {} - file not found", file);
                None
            }
        })
        .collect();

    if config.threads > 1 {
        // Parallel processing with rayon
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|e| format!("Failed to create thread pool: {}", e))?;

        let errors: Vec<String> = pool.install(|| {
            use rayon::prelude::*;
            paths
                .par_iter()
                .filter_map(|path| {
                    if let Err(e) = process_single_file(path, config) {
                        Some(format!("{}: {}", path.display(), e))
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
        for path in &paths {
            process_single_file(path, config)?;
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
        let sequences =
            sequences.map(Ok::<_, Box<dyn std::error::Error>>);
        analysis::run_analysis(sequences, &mut modules, self.config.quiet, self.config.min_length)?;
        let data_report = report::generate_data_report(&mut modules, VERSION);
        let html_report = report::generate_html_report(&mut modules, "sequences", VERSION);

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
            modules: module_results,
        })
    }

    /// Run QC on a file path and return both reports.
    pub fn run_file(&self, path: &Path) -> Result<FastQCReport, Box<dyn std::error::Error>> {
        run_fastqc_on_file(path, &self.config)
    }
}
