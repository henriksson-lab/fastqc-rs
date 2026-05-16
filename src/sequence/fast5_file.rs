use std::path::Path;

use hdf5_pure_rust::{Dataset, File};
use thiserror::Error;

use crate::sequence::Sequence;

const FASTQ_DATASET_PATHS: [&str; 4] = [
    "Analyses/Basecall_2D_000/BaseCalled_template/Fastq",
    "Analyses/Basecall_2D_000/BaseCalled_2D/Fastq",
    "Analyses/Basecall_1D_000/BaseCalled_template/Fastq",
    "Analyses/Basecall_1D_000/BaseCalled_1D/Fastq",
];

/// Error type for Fast5 (Oxford Nanopore HDF5) parsing.
#[derive(Debug, Error)]
pub enum Fast5Error {
    #[error("HDF5 error: {0}")]
    Hdf5(#[from] hdf5_pure_rust::Error),
    #[error("Fast5 format error: {0}")]
    Format(String),
}

/// Iterator that yields Sequence records from a Fast5 (HDF5) basecalled file.
pub struct Fast5FileReader {
    file: File,
    file_name: String,
    read_prefixes: Vec<String>,
    next_index: usize,
}

impl Fast5FileReader {
    /// Open a Fast5 file and enumerate its read groups.
    pub fn open(path: &Path) -> Result<Self, Fast5Error> {
        let file = File::open(path)?;
        let file_name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| path.display().to_string());
        let read_prefixes = read_prefixes(&file)?;

        Ok(Self {
            file,
            file_name,
            read_prefixes,
            next_index: 0,
        })
    }
}

impl Iterator for Fast5FileReader {
    type Item = Result<Sequence, Fast5Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_index >= self.read_prefixes.len() {
            return None;
        }

        let prefix = self.read_prefixes[self.next_index].clone();
        self.next_index += 1;

        Some(read_sequence(&self.file, &prefix, &self.file_name))
    }
}

/// List the `read_*` group prefixes at the root, or a single empty prefix for legacy single-read files.
fn read_prefixes(file: &File) -> Result<Vec<String>, Fast5Error> {
    let read_folders: Vec<String> = file
        .member_names()?
        .into_iter()
        .filter(|name| name.starts_with("read_"))
        .map(|name| format!("{}/", name))
        .collect();

    if read_folders.is_empty() {
        Ok(vec![String::new()])
    } else {
        Ok(read_folders)
    }
}

/// Locate an embedded FASTQ dataset for one read group and decode it into a Sequence.
fn read_sequence(file: &File, prefix: &str, file_name: &str) -> Result<Sequence, Fast5Error> {
    for suffix in FASTQ_DATASET_PATHS {
        let path = format!("{}{}", prefix, suffix);
        if let Ok(dataset) = file.dataset(&path) {
            let fastq = read_dataset_string(&dataset)?;
            return parse_embedded_fastq(&fastq, file_name);
        }
    }

    Err(Fast5Error::Format(format!(
        "No valid fastq paths found in Fast5 read '{}'",
        prefix
    )))
}

/// Read a HDF5 dataset as a UTF-8 string, falling back to raw bytes for non-string-typed datasets.
fn read_dataset_string(dataset: &Dataset) -> Result<String, Fast5Error> {
    if let Ok(value) = dataset.read_string() {
        return Ok(value);
    }

    let bytes = dataset.read_raw()?;
    String::from_utf8(bytes)
        .map_err(|err| Fast5Error::Format(format!("Fast5 Fastq dataset was not UTF-8: {}", err)))
}

/// Parse the 4-line FASTQ record embedded in a Fast5 dataset into a Sequence.
fn parse_embedded_fastq(fastq: &str, file_name: &str) -> Result<Sequence, Fast5Error> {
    let lines: Vec<&str> = fastq.trim_end_matches(['\r', '\n']).lines().collect();

    if lines.len() != 4 {
        return Err(Fast5Error::Format(format!(
            "Didn't get 4 sections from embedded FastQ: {}",
            fastq
        )));
    }
    if !lines[0].starts_with('@') {
        return Err(Fast5Error::Format(format!(
            "Embedded FastQ header did not start with '@': {}",
            lines[0]
        )));
    }
    if !lines[2].starts_with('+') {
        return Err(Fast5Error::Format(format!(
            "Embedded FastQ separator did not start with '+': {}",
            lines[2]
        )));
    }
    if lines[1].len() != lines[3].len() {
        return Err(Fast5Error::Format(format!(
            "Embedded FastQ sequence length ({}) does not match quality length ({})",
            lines[1].len(),
            lines[3].len()
        )));
    }

    Ok(Sequence::new(
        file_name.to_string(),
        lines[1].to_string(),
        lines[3].to_string(),
        lines[0].to_string(),
    ))
}
