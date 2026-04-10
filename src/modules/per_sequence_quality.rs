use std::collections::HashMap;

use crate::charts::{ChartData, Series};
use crate::encoding::phred::PhredEncoding;
use crate::modules::module_config::ModuleConfig;
use crate::modules::QCModule;
use crate::sequence::Sequence;

pub struct PerSequenceQualityScores {
    average_score_counts: HashMap<i32, u64>,
    quality_distribution: Vec<f64>,
    x_categories: Vec<i32>,
    lowest_char: char,
    max_count: f64,
    most_frequent_score: i32,
    calculated: bool,
    config: ModuleConfig,
}

impl PerSequenceQualityScores {
    pub fn new(config: ModuleConfig) -> Self {
        Self {
            average_score_counts: HashMap::new(),
            quality_distribution: Vec::new(),
            x_categories: Vec::new(),
            lowest_char: 126 as char,
            max_count: 0.0,
            most_frequent_score: 0,
            calculated: false,
            config,
        }
    }

    fn calculate_distribution(&mut self) {
        if self.average_score_counts.is_empty() {
            self.calculated = true;
            return;
        }

        let encoding = PhredEncoding::get_fastq_encoding_offset(self.lowest_char)
            .unwrap_or(PhredEncoding {
                name: "Unknown".to_string(),
                offset: 33,
            });

        let mut raw_scores: Vec<i32> = self.average_score_counts.keys().copied().collect();
        raw_scores.sort();

        let min_score = raw_scores[0];
        let max_score = raw_scores[raw_scores.len() - 1];
        let len = (1 + max_score - min_score) as usize;

        self.quality_distribution = vec![0.0; len];
        self.x_categories = vec![0; len];

        for i in 0..len {
            let score = min_score + i as i32;
            self.x_categories[i] = score - encoding.offset;
            if let Some(&count) = self.average_score_counts.get(&score) {
                self.quality_distribution[i] = count as f64;
            }
        }

        self.max_count = 0.0;
        for i in 0..self.quality_distribution.len() {
            if self.quality_distribution[i] > self.max_count {
                self.max_count = self.quality_distribution[i];
                self.most_frequent_score = self.x_categories[i];
            }
        }

        self.calculated = true;
    }
}

impl QCModule for PerSequenceQualityScores {
    fn name(&self) -> &str {
        "Per sequence quality scores"
    }

    fn description(&self) -> &str {
        "Shows the distribution of average quality scores for whole sequences"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("quality_sequence", "ignore") > 0.0
            || self.average_score_counts.is_empty()
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;
        let quality = sequence.quality.as_bytes();
        let mut average_quality: i32 = 0;

        for &q in quality {
            if (q as char) < self.lowest_char {
                self.lowest_char = q as char;
            }
            average_quality += q as i32;
        }

        if !quality.is_empty() {
            // Java integer division
            average_quality /= quality.len() as i32;
            *self
                .average_score_counts
                .entry(average_quality)
                .or_insert(0) += 1;
        }
    }

    fn reset(&mut self) {
        self.average_score_counts.clear();
        self.lowest_char = 126 as char;
        self.max_count = 0.0;
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.calculate_distribution();
        }
        self.most_frequent_score <= self.config.get_param("quality_sequence", "error") as i32
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.calculate_distribution();
        }
        self.most_frequent_score <= self.config.get_param("quality_sequence", "warn") as i32
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.calculate_distribution();
        }

        buf.push_str("#Quality\tCount\n");
        for i in 0..self.x_categories.len() {
            buf.push_str(&self.x_categories[i].to_string());
            buf.push('\t');
            buf.push_str(&format_java_double(self.quality_distribution[i]));
            buf.push('\n');
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.calculate_distribution();
        }
        if self.quality_distribution.is_empty() {
            return None;
        }

        let x_labels: Vec<String> = self.x_categories.iter().map(|q| q.to_string()).collect();
        let max_y = self.max_count;

        Some(ChartData::LineGraph {
            series: vec![Series {
                name: "Average Quality per read".to_string(),
                data: self.quality_distribution.clone(),
            }],
            x_labels,
            x_axis_label: "Mean Sequence Quality (Phred Score)".to_string(),
            title: "Quality score distribution over all sequences".to_string(),
            min_y: 0.0,
            max_y,
        })
    }
}

/// Format a f64 to match Java's Double.toString() behavior.
/// Java always includes ".0" for integer values (e.g., "1.0", "100.0").
/// Java uses scientific notation for |v| < 10^-3 or |v| >= 10^7.
pub fn format_java_double(v: f64) -> String {
    if !v.is_finite() {
        return format!("{}", v);
    }
    let abs = v.abs();
    if abs == 0.0 {
        return "0.0".to_string();
    }
    // Java Double.toString() uses scientific notation for very small/large values
    if !(1e-3..1e7).contains(&abs) {
        // Java format: e.g., "5.0E-4", "1.5E7"
        return format_java_scientific(v);
    }
    if v.fract() == 0.0 {
        format!("{:.1}", v)
    } else {
        format!("{}", v)
    }
}

/// Format like Java's scientific notation: e.g., 0.0005 -> "5.0E-4"
fn format_java_scientific(v: f64) -> String {
    // Java's Double.toString uses format like "5.0E-4" or "1.5E7"
    let formatted = format!("{:E}", v);
    // Rust formats as "5E-4", Java as "5.0E-4" — ensure .0 for integers
    // Also Rust uses "e" lowercase by default with {:e}, or uppercase with {:E}
    // Java always uses uppercase E
    if let Some(pos) = formatted.find('E') {
        let mantissa = &formatted[..pos];
        let exponent = &formatted[pos + 1..];
        // Ensure mantissa has decimal point
        let mantissa = if mantissa.contains('.') {
            mantissa.to_string()
        } else {
            format!("{}.0", mantissa)
        };
        format!("{}E{}", mantissa, exponent)
    } else {
        formatted
    }
}
