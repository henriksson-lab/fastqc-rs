use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::QualityScores;

use crate::sequence::Sequence;

/// Error type for BAM/SAM parsing.
#[derive(Debug, thiserror::Error)]
pub enum BamError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("BAM format error: {0}")]
    Format(String),
}

/// Iterator that yields Sequence records from a BAM file.
pub struct BamFileReader {
    reader: noodles::bam::io::Reader<noodles::bgzf::Reader<BufReader<File>>>,
    record: noodles::bam::Record,
    file_name: String,
    only_mapped: bool,
    finished: bool,
}

impl BamFileReader {
    /// Open a BAM file and consume its header. If `only_mapped` is true, unmapped reads
    /// are skipped and soft-clipped portions are trimmed from the returned sequence.
    pub fn open(path: &Path, only_mapped: bool) -> Result<Self, BamError> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let file = File::open(path)?;
        let mut reader = noodles::bam::io::Reader::new(BufReader::new(file));
        let _header = reader
            .read_header()
            .map_err(|e| BamError::Format(e.to_string()))?;

        Ok(Self {
            reader,
            record: noodles::bam::Record::default(),
            file_name: name,
            only_mapped,
            finished: false,
        })
    }
}

impl Iterator for BamFileReader {
    type Item = Result<Sequence, BamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        loop {
            match self.reader.read_record(&mut self.record) {
                Ok(0) => {
                    self.finished = true;
                    return None;
                }
                Ok(_) => {}
                Err(e) => {
                    self.finished = true;
                    return Some(Err(BamError::Format(e.to_string())));
                }
            }

            let flags = self.record.flags();
            if self.only_mapped && flags.is_unmapped() {
                continue;
            }

            // Extract sequence
            let seq_data = self.record.sequence();
            let seq_len = seq_data.len();
            let mut sequence: String = (0..seq_len)
                .map(|i| seq_data.get(i).map(|b| b as char).unwrap_or('N'))
                .collect();

            // Extract quality scores
            let qual_data = self.record.quality_scores();
            let mut qualities: String = qual_data
                .iter()
                .filter_map(|r| r.ok())
                .map(|q| (q + 33) as char)
                .collect();

            // Pad quality string if shorter than sequence (can happen if quality parsing fails)
            if qualities.len() < sequence.len() {
                qualities.extend(std::iter::repeat_n('!', sequence.len() - qualities.len()));
            }

            // Handle soft clipping for mapped-only mode
            if self.only_mapped {
                let cigar = self.record.cigar();
                let ops: Vec<_> = cigar.iter().filter_map(|r| r.ok()).collect();
                if !ops.is_empty() {
                    let last = ops[ops.len() - 1];
                    if last.kind() == Kind::SoftClip {
                        let clip_len = last.len();
                        let slen = sequence.len();
                        if clip_len < slen {
                            sequence.truncate(slen - clip_len);
                            qualities.truncate(slen - clip_len);
                        }
                    }
                    let first = ops[0];
                    if first.kind() == Kind::SoftClip {
                        let clip_len = first.len();
                        if clip_len < sequence.len() {
                            sequence = sequence[clip_len..].to_string();
                            qualities = qualities[clip_len..].to_string();
                        }
                    }
                }
            }

            if flags.is_reverse_complemented() {
                sequence = reverse_complement(&sequence);
                qualities = reverse_string(&qualities);
            }

            let read_name = self
                .record
                .name()
                .map(|n| n.to_string())
                .unwrap_or_else(|| "unknown".to_string());

            return Some(Ok(Sequence::new(
                self.file_name.clone(),
                sequence,
                qualities,
                format!("@{}", read_name),
            )));
        }
    }
}

/// Iterator that yields Sequence records from a SAM file by parsing lines directly.
pub struct SamFileReader {
    reader: BufReader<File>,
    file_name: String,
    only_mapped: bool,
    finished: bool,
}

impl SamFileReader {
    /// Open a SAM file for line-by-line parsing. If `only_mapped` is true, unmapped reads
    /// are skipped and soft-clipped portions are trimmed from the returned sequence.
    pub fn open(path: &Path, only_mapped: bool) -> Result<Self, BamError> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let file = File::open(path)?;

        Ok(Self {
            reader: BufReader::new(file),
            file_name: name,
            only_mapped,
            finished: false,
        })
    }
}

impl Iterator for SamFileReader {
    type Item = Result<Sequence, BamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        loop {
            let mut buf = String::new();
            match self.reader.read_line(&mut buf) {
                Ok(0) => {
                    self.finished = true;
                    return None;
                }
                Ok(_) => {}
                Err(e) => {
                    self.finished = true;
                    return Some(Err(BamError::Io(e)));
                }
            }

            let line = buf.trim();
            if line.is_empty() || line.starts_with('@') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 11 {
                return Some(Err(BamError::Format(format!(
                    "SAM line has {} fields, expected at least 11",
                    fields.len()
                ))));
            }

            let read_name = fields[0];
            let flag: u16 = fields[1].parse().unwrap_or(0);
            let flags = Flags::from(flag);

            if self.only_mapped && flags.is_unmapped() {
                continue;
            }

            let cigar_str = fields[5];
            let mut sequence = fields[9].to_string();
            let mut qualities = fields[10].to_string();

            if self.only_mapped && cigar_str != "*" {
                let ops = parse_cigar(cigar_str);
                if !ops.is_empty() {
                    let last = ops[ops.len() - 1];
                    if last.0 == 'S' {
                        let clip_len = last.1;
                        let seq_len = sequence.len();
                        if clip_len < seq_len {
                            sequence.truncate(seq_len - clip_len);
                            qualities.truncate(seq_len - clip_len);
                        }
                    }
                    let first = ops[0];
                    if first.0 == 'S' {
                        let clip_len = first.1;
                        if clip_len < sequence.len() {
                            sequence = sequence[clip_len..].to_string();
                            qualities = qualities[clip_len..].to_string();
                        }
                    }
                }
            }

            if flags.is_reverse_complemented() {
                sequence = reverse_complement(&sequence);
                qualities = reverse_string(&qualities);
            }

            return Some(Ok(Sequence::new(
                self.file_name.clone(),
                sequence,
                qualities,
                format!("@{}", read_name),
            )));
        }
    }
}

/// Parse a CIGAR string into (op, length) pairs (e.g. "10S40M" -> [('S',10),('M',40)]).
fn parse_cigar(cigar: &str) -> Vec<(char, usize)> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            if let Ok(len) = num.parse::<usize>() {
                ops.push((ch, len));
            }
            num.clear();
        }
    }
    ops
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'G' => 'C',
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            other => other,
        })
        .collect()
}

fn reverse_string(s: &str) -> String {
    s.chars().rev().collect()
}
