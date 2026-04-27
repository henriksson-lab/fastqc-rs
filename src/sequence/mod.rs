pub mod bam_file;
pub mod fast5_file;
pub mod fastq_file;

/// A single sequence read with its quality scores and metadata.
#[derive(Debug, Clone)]
pub struct Sequence {
    /// Sequence identifier (from @ line in FASTQ)
    pub id: String,
    /// DNA bases (uppercase)
    pub sequence: String,
    /// Quality scores (ASCII encoded)
    pub quality: String,
    /// Source file name
    pub file_name: String,
    /// Original colorspace data if applicable
    pub colorspace: Option<String>,
    /// Whether this sequence was flagged as filtered (CASAVA Y flag)
    pub is_filtered: bool,
}

impl Sequence {
    /// Create a new Sequence with standard base calls.
    /// Sequence is converted to uppercase, matching Java behavior.
    pub fn new(file_name: String, sequence: String, quality: String, id: String) -> Self {
        Self {
            id,
            sequence: sequence.to_uppercase(),
            quality,
            file_name,
            colorspace: None,
            is_filtered: false,
        }
    }

    /// Create a new Sequence with colorspace data.
    pub fn new_with_colorspace(
        file_name: String,
        sequence: String,
        colorspace: String,
        quality: String,
        id: String,
    ) -> Self {
        Self {
            id,
            sequence,
            quality,
            file_name,
            colorspace: Some(colorspace),
            is_filtered: false,
        }
    }
}
