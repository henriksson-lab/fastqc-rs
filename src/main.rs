use clap::Parser;
use std::path::PathBuf;
use std::process;

use fastqc_rs::config::FastQCConfig;
use fastqc_rs::run_fastqc;

/// FastQC - A quality control tool for high throughput sequence data
#[derive(Parser, Debug)]
#[command(name = "fastqc-rs", version, about)]
struct Cli {
    /// Input files (FASTQ, BAM, SAM). Use "stdin" to read from standard input.
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

    /// Delete the zip file after extracting
    #[arg(long = "delete")]
    delete: bool,

    /// Minimum length of sequence to process
    #[arg(long = "min_length", default_value = "0")]
    min_length: usize,

    /// Length to truncate sequences for duplication detection
    #[arg(long = "dup_length", default_value = "50")]
    dup_length: usize,

    /// Generate SVG images instead of PNG
    #[arg(long = "svg")]
    svg: bool,
}

fn main() {
    let cli = Cli::parse();

    if let Some(ref fmt) = cli.format {
        let valid = ["fastq", "bam", "sam", "bam_mapped", "sam_mapped"];
        if !valid.contains(&fmt.as_str()) {
            eprintln!(
                "Error: unrecognized format '{}'. Valid formats: {}",
                fmt,
                valid.join(", ")
            );
            process::exit(1);
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
    };

    match run_fastqc(&cli.files, &config) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    }
}
