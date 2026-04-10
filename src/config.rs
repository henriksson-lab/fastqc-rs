use std::path::PathBuf;

/// Configuration for FastQC analysis, mirroring Java FastQCConfig.
#[derive(Debug, Clone)]
pub struct FastQCConfig {
    /// Don't group bases for reads >50bp
    pub nogroup: bool,
    /// Use exponential base grouping
    pub expgroup: bool,
    /// Suppress progress messages
    pub quiet: bool,
    /// K-mer size for KmerContent module (default 7)
    pub kmer_size: usize,
    /// Number of threads for parallel file processing
    pub threads: usize,
    /// Output directory
    pub output_dir: Option<PathBuf>,
    /// CASAVA mode: group files and filter
    pub casava: bool,
    /// Nanopore mode: scan directories for Fast5
    pub nano: bool,
    /// Don't filter CASAVA reads
    pub nofilter: bool,
    /// Unzip output archives
    pub do_unzip: bool,
    /// Delete ZIP after unzipping
    pub delete_after_unzip: bool,
    /// Force input format (fastq, bam, sam, bam_mapped, sam_mapped)
    pub sequence_format: Option<String>,
    /// Custom contaminant file
    pub contaminant_file: Option<PathBuf>,
    /// Custom adapter file
    pub adapter_file: Option<PathBuf>,
    /// Custom limits file
    pub limits_file: Option<PathBuf>,
    /// Minimum sequence length to include
    pub min_length: usize,
    /// Length to truncate sequences for duplication detection
    pub dup_length: usize,
    /// Generate SVG images instead of PNG
    pub svg_output: bool,
}

impl Default for FastQCConfig {
    fn default() -> Self {
        Self {
            nogroup: false,
            expgroup: false,
            quiet: false,
            kmer_size: 7,
            threads: 1,
            output_dir: None,
            casava: false,
            nano: false,
            nofilter: false,
            do_unzip: false,
            delete_after_unzip: false,
            sequence_format: None,
            contaminant_file: None,
            adapter_file: None,
            limits_file: None,
            min_length: 0,
            dup_length: 50,
            svg_output: false,
        }
    }
}
