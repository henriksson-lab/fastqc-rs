use clap::Parser;
use std::ffi::OsString;
use std::path::PathBuf;

use crate::config::FastQCConfig;
use crate::run_fastqc;

/// FastQC - A quality control tool for high throughput sequence data
#[derive(Parser, Debug)]
#[command(name = "fastqc-compliant-rs", version, about)]
struct Cli {
    /// Input files (FASTQ, BAM, SAM, Fast5). Use "stdin" to read from standard input.
    #[arg(required = true)]
    files: Vec<String>,

    /// Create output files in specified directory
    #[arg(short = 'o', long = "outdir")]
    output_dir: Option<PathBuf>,

    /// Number of threads to use for processing files in parallel
    #[arg(short = 't', long = "threads", default_value = "1")]
    threads: usize,

    /// Custom contaminant list file
    #[arg(short = 'c', long = "contaminants")]
    contaminants: Option<PathBuf>,

    /// Custom adapter list file
    #[arg(short = 'a', long = "adapters")]
    adapters: Option<PathBuf>,

    /// Custom limits file for warn/error thresholds
    #[arg(short = 'l', long = "limits")]
    limits: Option<PathBuf>,

    /// Length of kmer to look for (2-10, default 7)
    #[arg(short = 'k', long = "kmers", default_value = "7")]
    kmers: usize,

    /// Suppress all progress messages
    #[arg(short = 'q', long = "quiet")]
    quiet: bool,

    /// Disable grouping of bases for reads >50bp
    #[arg(long = "nogroup")]
    nogroup: bool,

    /// Use exponential base grouping
    #[arg(long = "expgroup")]
    expgroup: bool,

    /// Force input file format (fastq, bam, sam, bam_mapped, sam_mapped)
    #[arg(short = 'f', long = "format")]
    format: Option<String>,

    /// Files come from raw CASAVA output (group and filter)
    #[arg(long = "casava")]
    casava: bool,

    /// Files come from Nanopore sequencing (scan directories)
    #[arg(long = "nano")]
    nano: bool,

    /// Don't remove poor quality reads when in CASAVA mode
    #[arg(long = "nofilter")]
    nofilter: bool,

    /// Unzip the output file after creating it
    #[arg(long = "extract")]
    extract: bool,

    /// Do not unzip the output file after creating it
    #[arg(long = "noextract", conflicts_with = "extract")]
    noextract: bool,

    /// Delete the zip file after extracting
    #[arg(long = "delete")]
    delete: bool,

    /// Base memory per file in megabytes. Accepted for FastQC CLI compatibility.
    #[arg(long = "memory")]
    memory: Option<usize>,

    /// Temporary directory. Accepted for FastQC CLI compatibility.
    #[arg(short = 'd', long = "dir")]
    temp_dir: Option<PathBuf>,

    /// Java executable. Accepted for FastQC CLI compatibility and ignored.
    #[arg(short = 'j', long = "java")]
    java: Option<PathBuf>,

    /// Minimum length of sequence to process
    #[arg(long = "min_length", default_value = "0")]
    min_length: usize,

    /// Length to truncate sequences for duplication detection
    #[arg(long = "dup_length", default_value = "50")]
    dup_length: usize,

    /// Generate SVG images instead of PNG
    #[arg(long = "svg")]
    svg: bool,

    /// Embed icons and charts as base64 data URIs in HTML
    #[arg(long = "embed-images")]
    embed_images: bool,
}

pub fn run_cli_from<I, T>(args: I) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let cli = Cli::parse_from(args);

    if cli.threads == 0 {
        eprintln!("Error: number of threads must be a positive integer");
        return 1;
    }

    if !(2..=10).contains(&cli.kmers) {
        eprintln!("Error: kmer size must be in the range 2-10");
        return 1;
    }

    if cli.nogroup && cli.expgroup {
        eprintln!("Error: you can't specify both --expgroup and --nogroup in the same run");
        return 1;
    }

    if let Some(ref fmt) = cli.format {
        let valid = ["fastq", "bam", "sam", "bam_mapped", "sam_mapped"];
        if !valid.contains(&fmt.as_str()) {
            eprintln!(
                "Error: unrecognized format '{}'. Valid formats: {}",
                fmt,
                valid.join(", ")
            );
            return 1;
        }
    }

    if let Some(ref output_dir) = cli.output_dir {
        if !output_dir.is_dir() {
            eprintln!(
                "Error: specified output directory '{}' does not exist",
                output_dir.display()
            );
            return 1;
        }
    }

    for (label, path) in [
        ("contaminant file", cli.contaminants.as_ref()),
        ("adapter file", cli.adapters.as_ref()),
        ("limits file", cli.limits.as_ref()),
    ] {
        if let Some(path) = path {
            if !path.is_file() {
                eprintln!("Error: {} '{}' does not exist", label, path.display());
                return 1;
            }
        }
    }

    if let Some(ref temp_dir) = cli.temp_dir {
        if !temp_dir.is_dir() {
            eprintln!(
                "Error: temp directory '{}' does not exist",
                temp_dir.display()
            );
            return 1;
        }
    }

    if let Some(memory) = cli.memory {
        if !(100..=10000).contains(&memory) {
            eprintln!(
                "Error: memory value {} MB was outside the allowed range (100 - 10000)",
                memory
            );
            return 1;
        }
    }

    let config = FastQCConfig {
        nogroup: cli.nogroup,
        expgroup: cli.expgroup,
        quiet: cli.quiet,
        kmer_size: cli.kmers,
        threads: cli.threads,
        output_dir: cli.output_dir,
        casava: cli.casava,
        nano: cli.nano,
        nofilter: cli.nofilter,
        do_unzip: cli.extract,
        delete_after_unzip: cli.delete,
        sequence_format: cli.format,
        contaminant_file: cli.contaminants,
        adapter_file: cli.adapters,
        limits_file: cli.limits,
        min_length: cli.min_length,
        dup_length: cli.dup_length,
        svg_output: cli.svg,
        embed_images: cli.embed_images,
    };

    match run_fastqc(&cli.files, &config) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error: {}", e);
            return 1;
        }
    }

    0
}
