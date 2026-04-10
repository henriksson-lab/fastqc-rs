use crate::charts::ChartData;
use crate::config::FastQCConfig;
use crate::encoding::phred::PhredEncoding;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::utils::base_group::BaseGroup;
use crate::utils::quality_count::QualityCount;

pub struct PerBaseQuality {
    quality_counts: Vec<QualityCount>,
    means: Vec<f64>,
    medians: Vec<f64>,
    lower_quartile: Vec<f64>,
    upper_quartile: Vec<f64>,
    lowest: Vec<f64>,
    highest: Vec<f64>,
    x_labels: Vec<String>,
    calculated: bool,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
}

impl PerBaseQuality {
    pub fn new(config: ModuleConfig, fqc_config: FastQCConfig) -> Self {
        Self {
            quality_counts: Vec::new(),
            means: Vec::new(),
            medians: Vec::new(),
            lower_quartile: Vec::new(),
            upper_quartile: Vec::new(),
            lowest: Vec::new(),
            highest: Vec::new(),
            x_labels: Vec::new(),
            calculated: false,
            config,
            fqc_config,
        }
    }

    fn calculate_offsets(&self) -> (char, char) {
        let mut min_char: char = '\0';
        let mut max_char: char = '\0';

        for (q, qc) in self.quality_counts.iter().enumerate() {
            if q == 0 {
                min_char = qc.get_min_char();
                max_char = qc.get_max_char();
            } else {
                if qc.get_min_char() < min_char {
                    min_char = qc.get_min_char();
                }
                if qc.get_max_char() > max_char {
                    max_char = qc.get_max_char();
                }
            }
        }

        (min_char, max_char)
    }

    fn get_percentages(&mut self) {
        let (min_char, _max_char) = self.calculate_offsets();
        let encoding = PhredEncoding::get_fastq_encoding_offset(min_char)
            .unwrap_or(PhredEncoding {
                name: "Unknown".to_string(),
                offset: 33,
            });
        let offset = encoding.offset;

        let groups = BaseGroup::make_base_groups(self.quality_counts.len(), &self.fqc_config);

        self.means = vec![0.0; groups.len()];
        self.medians = vec![0.0; groups.len()];
        self.lowest = vec![0.0; groups.len()];
        self.highest = vec![0.0; groups.len()];
        self.lower_quartile = vec![0.0; groups.len()];
        self.upper_quartile = vec![0.0; groups.len()];
        self.x_labels = Vec::with_capacity(groups.len());

        for (i, group) in groups.iter().enumerate() {
            self.x_labels.push(group.to_string());
            let min_base = group.lower_count;
            let max_base = group.upper_count;
            self.lowest[i] = self.get_percentile(min_base, max_base, offset, 10);
            self.highest[i] = self.get_percentile(min_base, max_base, offset, 90);
            self.means[i] = self.get_mean(min_base, max_base, offset);
            self.medians[i] = self.get_percentile(min_base, max_base, offset, 50);
            self.lower_quartile[i] = self.get_percentile(min_base, max_base, offset, 25);
            self.upper_quartile[i] = self.get_percentile(min_base, max_base, offset, 75);
        }

        self.calculated = true;
    }

    fn get_percentile(&self, min_bp: usize, max_bp: usize, offset: i32, percentile: i32) -> f64 {
        let mut count = 0;
        let mut total = 0.0;

        // Java: for (int i=minbp-1;i<maxbp;i++)
        for i in (min_bp - 1)..max_bp {
            if self.quality_counts[i].get_total_count() > 100 {
                count += 1;
                total += self.quality_counts[i].get_percentile(offset, percentile);
            }
        }

        if count > 0 {
            total / count as f64
        } else {
            f64::NAN
        }
    }

    fn get_mean(&self, min_bp: usize, max_bp: usize, offset: i32) -> f64 {
        let mut count = 0;
        let mut total = 0.0;

        for i in (min_bp - 1)..max_bp {
            if self.quality_counts[i].get_total_count() > 0 {
                count += 1;
                total += self.quality_counts[i].get_mean(offset);
            }
        }

        if count > 0 {
            total / count as f64
        } else {
            0.0
        }
    }
}

impl QCModule for PerBaseQuality {
    fn name(&self) -> &str {
        "Per base sequence quality"
    }

    fn description(&self) -> &str {
        "Shows the Quality scores of all bases at a given position in a sequencing run"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("quality_base", "ignore") > 0.0 || self.quality_counts.is_empty()
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;
        let qual = &sequence.quality;

        if self.quality_counts.len() < qual.len() {
            self.quality_counts.resize_with(qual.len(), QualityCount::new);
        }

        for (i, ch) in qual.chars().enumerate() {
            self.quality_counts[i].add_value(ch);
        }
    }

    fn reset(&mut self) {
        self.quality_counts.clear();
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }

        for i in 0..self.lower_quartile.len() {
            if self.lower_quartile[i].is_nan() {
                continue;
            }
            if self.lower_quartile[i] < self.config.get_param("quality_base_lower", "error")
                || self.medians[i] < self.config.get_param("quality_base_median", "error")
            {
                return true;
            }
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }

        for i in 0..self.lower_quartile.len() {
            if self.lower_quartile[i].is_nan() {
                continue;
            }
            if self.lower_quartile[i] < self.config.get_param("quality_base_lower", "warn")
                || self.medians[i] < self.config.get_param("quality_base_median", "warn")
            {
                return true;
            }
        }
        false
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.get_percentages();
        }

        buf.push_str("#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n");
        for i in 0..self.means.len() {
            buf.push_str(&self.x_labels[i]);
            buf.push('\t');
            buf.push_str(&format_java_double(self.means[i]));
            buf.push('\t');
            buf.push_str(&format_java_double(self.medians[i]));
            buf.push('\t');
            buf.push_str(&format_java_double(self.lower_quartile[i]));
            buf.push('\t');
            buf.push_str(&format_java_double(self.upper_quartile[i]));
            buf.push('\t');
            buf.push_str(&format_java_double(self.lowest[i]));
            buf.push('\t');
            buf.push_str(&format_java_double(self.highest[i]));
            buf.push('\n');
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.get_percentages();
        }
        if self.means.is_empty() {
            return None;
        }

        let (min_char, _) = self.calculate_offsets();
        let encoding = PhredEncoding::get_fastq_encoding_offset(min_char)
            .unwrap_or(PhredEncoding {
                name: "Unknown".to_string(),
                offset: 33,
            });
        // Sanger max quality ~41, Illumina max ~40
        let y_max = if encoding.offset == 33 { 41.0 } else { 40.0 };

        Some(ChartData::QualityBoxPlot {
            means: self.means.clone(),
            medians: self.medians.clone(),
            lower_quartile: self.lower_quartile.clone(),
            upper_quartile: self.upper_quartile.clone(),
            tenth_percentile: self.lowest.clone(),
            ninetieth_percentile: self.highest.clone(),
            x_labels: self.x_labels.clone(),
            y_min: 0.0,
            y_max,
        })
    }
}
