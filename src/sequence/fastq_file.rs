use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use crate::sequence::Sequence;

/// Error type for FASTQ parsing.
#[derive(Debug, thiserror::Error)]
pub enum FastQError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Format error at line {line}: {message}")]
    Format { line: u64, message: String },
}

/// FASTQ file reader supporting plain, gzip, and bzip2 formats.
/// Mirrors Java FastQFile behavior exactly.
pub struct FastQFile {
    reader: BufReader<Box<dyn Read>>,
    file_name: String,
    line_number: u64,
    is_colorspace: bool,
    colorspace_checked: bool,
    casava_mode: bool,
    nofilter: bool,
    next_sequence: Option<Sequence>,
    pending_error: Option<FastQError>,
    finished: bool,
}

impl FastQFile {
    /// Open a FASTQ file, auto-detecting compression from extension.
    pub fn open(path: &Path, casava: bool, nofilter: bool) -> Result<Self, FastQError> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let file = File::open(path)?;
        let reader: Box<dyn Read> = if name.to_lowercase().ends_with(".gz") {
            Box::new(flate2::read::MultiGzDecoder::new(file))
        } else if name.to_lowercase().ends_with(".bz2") {
            Box::new(bzip2::read::BzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let mut fq = Self {
            reader: BufReader::new(reader),
            file_name: name,
            line_number: 0,
            is_colorspace: false,
            colorspace_checked: false,
            casava_mode: casava,
            nofilter,
            next_sequence: None,
            pending_error: None,
            finished: false,
        };

        fq.read_next()?;
        Ok(fq)
    }

    /// Open from stdin.
    pub fn from_stdin(casava: bool, nofilter: bool) -> Result<Self, FastQError> {
        let reader: Box<dyn Read> = Box::new(io::stdin());
        let mut fq = Self {
            reader: BufReader::new(reader),
            file_name: "stdin".to_string(),
            line_number: 0,
            is_colorspace: false,
            colorspace_checked: false,
            casava_mode: casava,
            nofilter,
            next_sequence: None,
            pending_error: None,
            finished: false,
        };
        fq.read_next()?;
        Ok(fq)
    }

    pub fn name(&self) -> &str {
        &self.file_name
    }

    pub fn is_colorspace(&self) -> bool {
        self.is_colorspace
    }

    fn read_line(&mut self) -> Result<Option<String>, io::Error> {
        let mut line = String::new();
        let bytes = self.reader.read_line(&mut line)?;
        self.line_number += 1;
        if bytes == 0 {
            Ok(None)
        } else {
            // Remove trailing newline
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            Ok(Some(line))
        }
    }

    /// Read the next FASTQ record. Mirrors Java FastQFile.readNext() exactly.
    fn read_next(&mut self) -> Result<(), FastQError> {
        // Read ID line, skipping blank lines
        let id = loop {
            match self.read_line()? {
                None => {
                    self.next_sequence = None;
                    self.finished = true;
                    return Ok(());
                }
                Some(line) if line.is_empty() => continue,
                Some(line) => break line,
            }
        };

        if !id.starts_with('@') {
            self.next_sequence = None;
            return Err(FastQError::Format {
                line: self.line_number,
                message: "ID line didn't start with '@'".to_string(),
            });
        }

        // Read sequence line
        let seq = self.read_line()?.ok_or_else(|| FastQError::Format {
            line: self.line_number,
            message:
                "Ran out of data in the middle of a fastq entry. Your file is probably truncated"
                    .to_string(),
        })?;

        // Read mid line (separator starting with +)
        let mid_line = self.read_line()?.ok_or_else(|| FastQError::Format {
            line: self.line_number,
            message:
                "Ran out of data in the middle of a fastq entry. Your file is probably truncated"
                    .to_string(),
        })?;

        if !mid_line.starts_with('+') {
            return Err(FastQError::Format {
                line: self.line_number,
                message: format!("Midline '{}' didn't start with '+'", mid_line),
            });
        }

        // Read quality line
        let quality = self.read_line()?.ok_or_else(|| FastQError::Format {
            line: self.line_number,
            message:
                "Ran out of data in the middle of a fastq entry. Your file is probably truncated"
                    .to_string(),
        })?;

        // Validate sequence and quality lengths match
        if seq.len() != quality.len() {
            return Err(FastQError::Format {
                line: self.line_number,
                message: format!(
                    "Sequence length ({}) doesn't match quality length ({})",
                    seq.len(),
                    quality.len()
                ),
            });
        }

        // Check for colorspace on first sequence only
        if !self.colorspace_checked {
            self.colorspace_checked = true;
            self.is_colorspace = check_colorspace(&seq);
        }

        let sequence = if self.is_colorspace {
            let upper_seq = seq.to_uppercase();
            let converted = convert_colorspace_to_bases(&upper_seq);
            Sequence::new_with_colorspace(self.file_name.clone(), converted, upper_seq, quality, id)
        } else {
            Sequence::new(self.file_name.clone(), seq, quality, id)
        };

        let mut sequence = sequence;

        // CASAVA filtering: mark sequences with :Y: as filtered
        if self.casava_mode && !self.nofilter && sequence.id.contains(":Y:") {
            sequence.is_filtered = true;
        }

        self.next_sequence = Some(sequence);
        Ok(())
    }
}

impl Iterator for FastQFile {
    type Item = Result<Sequence, FastQError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        if let Some(err) = self.pending_error.take() {
            self.finished = true;
            self.next_sequence = None;
            return Some(Err(err));
        }

        let seq = self.next_sequence.take()?;
        if let Err(err) = self.read_next() {
            self.pending_error = Some(err);
        }
        Some(Ok(seq))
    }
}

/// Check if a sequence is colorspace data.
/// Matches Java: regex `^[GATCNgatcn][\\.0123456]+$`
fn check_colorspace(seq: &str) -> bool {
    if seq.len() < 2 {
        return false;
    }
    let bytes = seq.as_bytes();
    let first = bytes[0].to_ascii_uppercase();
    if !matches!(first, b'G' | b'A' | b'T' | b'C' | b'N') {
        return false;
    }
    bytes[1..]
        .iter()
        .all(|&b| matches!(b, b'.' | b'0' | b'1' | b'2' | b'3' | b'4' | b'5' | b'6'))
}

/// Convert colorspace data to base calls.
/// Exact port of Java FastQFile.convertColorspaceToBases().
fn convert_colorspace_to_bases(s: &str) -> String {
    let cs: Vec<char> = s.to_uppercase().chars().collect();
    if cs.is_empty() {
        return String::new();
    }

    let mut bp = vec![' '; cs.len() - 1];

    for i in 1..cs.len() {
        let ref_base = if i == 1 { cs[i - 1] } else { bp[i - 2] };

        if !matches!(ref_base, 'G' | 'A' | 'T' | 'C') {
            panic!(
                "Colourspace sequence data should always start with a real DNA letter. Line '{}' started with {} at position {}",
                s, ref_base, i
            );
        }

        match cs[i] {
            '0' => {
                bp[i - 1] = ref_base; // Same as reference
            }
            '1' => {
                bp[i - 1] = match ref_base {
                    'G' => 'T',
                    'A' => 'C',
                    'T' => 'G',
                    'C' => 'A',
                    _ => unreachable!(),
                };
            }
            '2' => {
                bp[i - 1] = match ref_base {
                    'G' => 'A',
                    'A' => 'G',
                    'T' => 'C',
                    'C' => 'T',
                    _ => unreachable!(),
                };
            }
            '3' => {
                bp[i - 1] = match ref_base {
                    'G' => 'C',
                    'A' => 'T',
                    'T' => 'A',
                    'C' => 'G',
                    _ => unreachable!(),
                };
            }
            '.' | '4' | '5' | '6' => {
                // Fill remaining with N
                for j in i..cs.len() {
                    bp[j - 1] = 'N';
                }
                break;
            }
            other => {
                panic!("Unexpected cs char {}", other);
            }
        }
    }

    bp.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_colorspace_check() {
        assert!(check_colorspace("A0123"));
        assert!(check_colorspace("G01230123"));
        assert!(!check_colorspace("ACGTACGT"));
        assert!(!check_colorspace(""));
    }

    #[test]
    fn test_colorspace_conversion() {
        // A0 -> A (same), A01 -> AC, etc.
        assert_eq!(convert_colorspace_to_bases("A0"), "A");
        assert_eq!(convert_colorspace_to_bases("A01"), "AC");
    }

    #[test]
    fn test_read_minimal_fastq() {
        let path = PathBuf::from("fastqc/test/data/minimal.fastq");
        if !path.exists() {
            return; // Skip if not in right directory
        }
        let fq = FastQFile::open(&path, false, false).unwrap();
        let sequences: Vec<Sequence> = fq.filter_map(|r| r.ok()).collect();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].id, "@READ0001");
        assert_eq!(sequences[0].sequence, "AAAAAAAAAAAAAAAA");
        assert_eq!(sequences[0].quality, "IIIIIIIIIIIIIIII");
    }

    #[test]
    fn test_read_complex_fastq() {
        let path = PathBuf::from("fastqc/test/data/complex.fastq");
        if !path.exists() {
            return;
        }
        let fq = FastQFile::open(&path, false, false).unwrap();
        let sequences: Vec<Sequence> = fq.filter_map(|r| r.ok()).collect();
        assert_eq!(sequences.len(), 5);
    }
}
