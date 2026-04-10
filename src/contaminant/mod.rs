use std::path::Path;

const DEFAULT_CONTAMINANTS: &str = include_str!("../resources/contaminant_list.txt");

/// Direction of a contaminant match.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Reverse,
}

/// Result of matching a sequence against a contaminant.
#[derive(Debug, Clone)]
pub struct ContaminantHit {
    pub contaminant_name: String,
    pub direction: Direction,
    pub length: usize,
    pub percent_id: i32,
}

impl std::fmt::Display for ContaminantHit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} ({}% over {}bp)",
            self.contaminant_name, self.percent_id, self.length
        )
    }
}

/// A known contaminant sequence with precomputed reverse complement.
#[derive(Debug, Clone)]
pub struct Contaminant {
    name: String,
    forward: Vec<u8>,
    reverse: Vec<u8>,
}

impl Contaminant {
    pub fn new(name: &str, sequence: &str) -> Self {
        let forward: Vec<u8> = sequence.to_uppercase().bytes().collect();
        let mut reverse = vec![0u8; forward.len()];
        for (c, &base) in forward.iter().enumerate() {
            let rev_pos = forward.len() - 1 - c;
            reverse[rev_pos] = match base {
                b'G' => b'C',
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                _ => panic!(
                    "Contaminant contained the illegal character '{}'",
                    base as char
                ),
            };
        }
        Self {
            name: name.to_string(),
            forward,
            reverse,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Find the best match of this contaminant against a query sequence.
    /// Matches Java Contaminant.findMatch() exactly.
    pub fn find_match(&self, query: &str) -> Option<ContaminantHit> {
        let query = query.to_uppercase();

        // Special case for queries between 8-20bp: exact substring match
        if query.len() < 20 && query.len() >= 8 {
            let forward_str =
                std::str::from_utf8(&self.forward).unwrap();
            let reverse_str =
                std::str::from_utf8(&self.reverse).unwrap();

            if forward_str.contains(&query) {
                return Some(ContaminantHit {
                    contaminant_name: self.name.clone(),
                    direction: Direction::Forward,
                    length: query.len(),
                    percent_id: 100,
                });
            }
            if reverse_str.contains(&query) {
                return Some(ContaminantHit {
                    contaminant_name: self.name.clone(),
                    direction: Direction::Reverse,
                    length: query.len(),
                    percent_id: 100,
                });
            }
        }

        let q: Vec<u8> = query.bytes().collect();
        let mut best_hit: Option<ContaminantHit> = None;

        // Search forward strand
        let min_offset = -(self.forward.len() as i64 - 20);
        let max_offset = q.len() as i64 - 20;
        for offset in min_offset..max_offset {
            if let Some(hit) = self.find_match_inner(&self.forward, &q, offset, Direction::Forward)
            {
                if best_hit.is_none() || hit.length > best_hit.as_ref().unwrap().length {
                    best_hit = Some(hit);
                }
            }
        }

        // Search reverse strand
        for offset in min_offset..max_offset {
            if let Some(hit) = self.find_match_inner(&self.reverse, &q, offset, Direction::Reverse)
            {
                if best_hit.is_none() || hit.length > best_hit.as_ref().unwrap().length {
                    best_hit = Some(hit);
                }
            }
        }

        best_hit
    }

    /// Inner matching logic. Exact port of Java Contaminant.findMatch(char[], char[], int, int).
    fn find_match_inner(
        &self,
        ca: &[u8],
        cb: &[u8],
        offset: i64,
        direction: Direction,
    ) -> Option<ContaminantHit> {
        let mut best_hit: Option<ContaminantHit> = None;
        let mut mismatch_count: i64 = 0;
        let mut start: i64 = 0;
        let mut end: i64 = 0;

        for (i, &ca_byte) in ca.iter().enumerate() {
            let i = i as i64;
            let j = i + offset;
            if j < 0 {
                start = i + 1;
                continue;
            }
            if j >= cb.len() as i64 {
                break;
            }

            if ca_byte == cb[j as usize] {
                end = i;
            } else {
                mismatch_count += 1;
                if mismatch_count > 1 {
                    // End of this match span
                    let span_len = 1 + end - start;
                    if span_len > 20 {
                        let id = ((span_len - (mismatch_count - 1)) * 100) / span_len;
                        let id = id as i32;
                        let span_len = span_len as usize;
                        if best_hit.is_none()
                            || best_hit.as_ref().unwrap().length < span_len
                            || (best_hit.as_ref().unwrap().length == span_len
                                && best_hit.as_ref().unwrap().percent_id < id)
                        {
                            best_hit = Some(ContaminantHit {
                                contaminant_name: self.name.clone(),
                                direction,
                                length: span_len,
                                percent_id: id,
                            });
                        }
                    }
                    start = i + 1;
                    end = i + 1;
                    mismatch_count = 0;
                }
            }
        }

        // Check if we ended with a match
        let span_len = 1 + end - start;
        if span_len > 20 {
            let id = ((span_len - mismatch_count) * 100) / span_len;
            let id = id as i32;
            let span_len = span_len as usize;
            if best_hit.is_none()
                || best_hit.as_ref().unwrap().length < span_len
                || (best_hit.as_ref().unwrap().length == span_len
                    && best_hit.as_ref().unwrap().percent_id < id)
            {
                best_hit = Some(ContaminantHit {
                    contaminant_name: self.name.clone(),
                    direction,
                    length: span_len,
                    percent_id: id,
                });
            }
        }

        best_hit
    }
}

/// Finds the best contaminant match for a given sequence.
pub struct ContaminantFinder {
    contaminants: Vec<Contaminant>,
}

impl ContaminantFinder {
    pub fn new(contaminant_file: Option<&Path>) -> Self {
        let text = if let Some(path) = contaminant_file {
            std::fs::read_to_string(path).unwrap_or_else(|e| {
                eprintln!(
                    "Warning: could not read contaminant file {:?}: {}",
                    path, e
                );
                DEFAULT_CONTAMINANTS.to_string()
            })
        } else {
            DEFAULT_CONTAMINANTS.to_string()
        };

        let contaminants = Self::parse_contaminant_list(&text);
        Self { contaminants }
    }

    fn parse_contaminant_list(text: &str) -> Vec<Contaminant> {
        let mut contaminants = Vec::new();
        for line in text.lines() {
            let line = line.trim();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let sections: Vec<&str> = line.splitn(2, '\t').collect();
            if sections.len() != 2 {
                eprintln!(
                    "Expected 2 sections for contaminant line but got {} from {}",
                    sections.len(),
                    line
                );
                continue;
            }
            contaminants.push(Contaminant::new(sections[0].trim(), sections[1].trim()));
        }
        contaminants
    }

    /// Find the best contaminant hit for a sequence.
    /// Matches Java ContaminentFinder.findContaminantHit() exactly.
    pub fn find_contaminant_hit(&self, sequence: &str) -> Option<ContaminantHit> {
        let mut best_hit: Option<ContaminantHit> = None;

        for contaminant in &self.contaminants {
            if let Some(hit) = contaminant.find_match(sequence) {
                if best_hit.is_none()
                    || hit.length > best_hit.as_ref().unwrap().length
                    || (hit.length == best_hit.as_ref().unwrap().length
                        && hit.percent_id > best_hit.as_ref().unwrap().percent_id)
                {
                    best_hit = Some(hit);
                }
            }
        }

        best_hit
    }
}
