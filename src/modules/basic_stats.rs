use crate::encoding::phred::PhredEncoding;
use crate::modules::module_config::ModuleConfig;
use crate::modules::QCModule;
use crate::sequence::Sequence;

pub struct BasicStats {
    name: Option<String>,
    actual_count: u64,
    filtered_count: u64,
    min_length: usize,
    max_length: usize,
    total_bases: u64,
    g_count: u64,
    c_count: u64,
    a_count: u64,
    t_count: u64,
    n_count: u64,
    lowest_char: u8,
    file_type: Option<String>,
}

impl BasicStats {
    pub fn new(_config: &ModuleConfig) -> Self {
        Self {
            name: None,
            actual_count: 0,
            filtered_count: 0,
            min_length: 0,
            max_length: 0,
            total_bases: 0,
            g_count: 0,
            c_count: 0,
            a_count: 0,
            t_count: 0,
            n_count: 0,
            lowest_char: 126,
            file_type: None,
        }
    }

    pub fn set_file_name(&mut self, name: &str) {
        self.name = Some(name.replacen("stdin:", "", 1));
    }

    pub fn lowest_char(&self) -> char {
        self.lowest_char as char
    }

    /// Format a base count with appropriate units.
    /// Exact port of Java BasicStats.formatLength().
    pub fn format_length(original_length: u64) -> String {
        let mut length = original_length as f64;
        let unit;

        if length >= 1_000_000_000.0 {
            length /= 1_000_000_000.0;
            unit = " Gbp";
        } else if length >= 1_000_000.0 {
            length /= 1_000_000.0;
            unit = " Mbp";
        } else if length >= 1000.0 {
            length /= 1000.0;
            unit = " kbp";
        } else {
            unit = " bp";
        }

        let raw = format!("{}", length);
        let chars: Vec<char> = raw.chars().collect();

        // Find the dot position
        let mut last_index = 0;
        for (i, &ch) in chars.iter().enumerate() {
            last_index = i;
            if ch == '.' {
                break;
            }
        }

        // Keep one more char if it's non-zero
        if last_index + 1 < chars.len() && chars[last_index + 1] != '0' {
            last_index += 1;
        } else if last_index > 0 && chars[last_index] == '.' {
            last_index -= 1; // Remove trailing dot
        }

        let formatted: String = chars[..=last_index].iter().collect();
        format!("{}{}", formatted, unit)
    }
}

impl QCModule for BasicStats {
    fn name(&self) -> &str {
        "Basic Statistics"
    }

    fn description(&self) -> &str {
        "Calculates some basic statistics about the file"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        false
    }

    fn ignore_in_report(&self) -> bool {
        false
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        if self.name.is_none() {
            self.set_file_name(&sequence.file_name);
        }

        if sequence.is_filtered {
            self.filtered_count += 1;
            return;
        }

        self.actual_count += 1;
        let seq_len = sequence.sequence.len();
        self.total_bases += seq_len as u64;

        if self.file_type.is_none() {
            self.file_type = Some(if sequence.colorspace.is_some() {
                "Colorspace converted to bases".to_string()
            } else {
                "Conventional base calls".to_string()
            });
        }

        if self.actual_count == 1 {
            self.min_length = seq_len;
            self.max_length = seq_len;
        } else {
            if seq_len < self.min_length {
                self.min_length = seq_len;
            }
            if seq_len > self.max_length {
                self.max_length = seq_len;
            }
        }

        for ch in sequence.sequence.bytes() {
            match ch {
                b'G' => self.g_count += 1,
                b'A' => self.a_count += 1,
                b'T' => self.t_count += 1,
                b'C' => self.c_count += 1,
                b'N' => self.n_count += 1,
                _ => {}
            }
        }

        for &ch in sequence.quality.as_bytes() {
            if ch < self.lowest_char {
                self.lowest_char = ch;
            }
        }
    }

    fn reset(&mut self) {
        self.min_length = 0;
        self.max_length = 0;
        self.g_count = 0;
        self.c_count = 0;
        self.a_count = 0;
        self.t_count = 0;
        self.n_count = 0;
    }

    fn raises_error(&mut self) -> bool {
        false
    }

    fn raises_warning(&mut self) -> bool {
        false
    }

    fn make_data_report(&mut self, buf: &mut String) {
        buf.push_str("#Measure\tValue\n");

        buf.push_str("Filename\t");
        buf.push_str(self.name.as_deref().unwrap_or("unknown"));
        buf.push('\n');

        buf.push_str("File type\t");
        buf.push_str(
            self.file_type
                .as_deref()
                .unwrap_or("Conventional base calls"),
        );
        buf.push('\n');

        buf.push_str("Encoding\t");
        let encoding = PhredEncoding::get_fastq_encoding_offset(self.lowest_char as char)
            .unwrap_or(PhredEncoding {
                name: "Unknown".to_string(),
                offset: 33,
            });
        buf.push_str(&encoding.name);
        buf.push('\n');

        buf.push_str("Total Sequences\t");
        buf.push_str(&self.actual_count.to_string());
        buf.push('\n');

        buf.push_str("Total Bases\t");
        buf.push_str(&Self::format_length(self.total_bases));
        buf.push('\n');

        buf.push_str("Sequences flagged as poor quality\t");
        buf.push_str(&self.filtered_count.to_string());
        buf.push('\n');

        buf.push_str("Sequence length\t");
        if self.min_length == self.max_length {
            buf.push_str(&self.min_length.to_string());
        } else {
            buf.push_str(&format!("{}-{}", self.min_length, self.max_length));
        }
        buf.push('\n');

        buf.push_str("%GC\t");
        let total_atgc = self.a_count + self.t_count + self.g_count + self.c_count;
        if total_atgc > 0 {
            // Java integer division: ((gCount+cCount)*100)/(aCount+tCount+gCount+cCount)
            let gc_percent = ((self.g_count + self.c_count) * 100) / total_atgc;
            buf.push_str(&gc_percent.to_string());
        } else {
            buf.push('0');
        }
        buf.push('\n');
    }
}
