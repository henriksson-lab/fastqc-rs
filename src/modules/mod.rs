pub mod adapter_content;
pub mod basic_stats;
pub mod duplication_level;
pub mod gc_model;
pub mod kmer_content;
pub mod module_config;
pub mod n_content;
pub mod overrepresented_seqs;
pub mod per_base_quality;
pub mod per_base_sequence_content;
pub mod per_sequence_gc_content;
pub mod per_sequence_quality;
pub mod per_tile_quality;
pub mod sequence_length_distribution;

use crate::charts::ChartData;
use crate::sequence::Sequence;

/// Status of a QC module's result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QCStatus {
    Pass,
    Warn,
    Fail,
}

impl std::fmt::Display for QCStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QCStatus::Pass => write!(f, "pass"),
            QCStatus::Warn => write!(f, "warn"),
            QCStatus::Fail => write!(f, "fail"),
        }
    }
}

/// Trait for all QC analysis modules, mirroring Java QCModule interface.
pub trait QCModule {
    /// Module display name (e.g. "Basic Statistics")
    fn name(&self) -> &str;

    /// Module description
    fn description(&self) -> &str;

    /// Process a single sequence read
    fn process_sequence(&mut self, sequence: &Sequence);

    /// Reset accumulated data
    fn reset(&mut self);

    /// Whether this module raises an error (FAIL status)
    fn raises_error(&mut self) -> bool;

    /// Whether this module raises a warning (WARN status)
    fn raises_warning(&mut self) -> bool;

    /// Whether to skip filtered sequences for this module
    fn ignore_filtered_sequences(&self) -> bool;

    /// Whether to exclude this module from the report
    fn ignore_in_report(&self) -> bool;

    /// Generate the text data report section for this module
    fn make_data_report(&mut self, buf: &mut String);

    /// Return chart data for rendering. Returns None for table-only modules.
    fn chart_data(&mut self) -> Option<ChartData> {
        None
    }

    /// Get the current status
    fn status(&mut self) -> QCStatus {
        if self.raises_error() {
            QCStatus::Fail
        } else if self.raises_warning() {
            QCStatus::Warn
        } else {
            QCStatus::Pass
        }
    }
}
