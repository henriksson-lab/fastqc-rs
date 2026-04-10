use crate::charts::{ChartData, Series};
use crate::config::FastQCConfig;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::utils::base_group::BaseGroup;

pub struct NContent {
    n_counts: Vec<u64>,
    not_n_counts: Vec<u64>,
    calculated: bool,
    percentages: Vec<f64>,
    x_categories: Vec<String>,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
}

impl NContent {
    pub fn new(config: ModuleConfig, fqc_config: FastQCConfig) -> Self {
        Self {
            n_counts: Vec::new(),
            not_n_counts: Vec::new(),
            calculated: false,
            percentages: Vec::new(),
            x_categories: Vec::new(),
            config,
            fqc_config,
        }
    }

    fn get_percentages(&mut self) {
        if self.n_counts.is_empty() {
            self.calculated = true;
            return;
        }
        let groups = BaseGroup::make_base_groups(self.n_counts.len(), &self.fqc_config);

        self.x_categories = Vec::with_capacity(groups.len());
        self.percentages = vec![0.0; groups.len()];

        for (i, group) in groups.iter().enumerate() {
            self.x_categories.push(group.to_string());

            let mut n_count: u64 = 0;
            let mut total: u64 = 0;

            for bp in (group.lower_count - 1)..group.upper_count {
                n_count += self.n_counts[bp];
                total += self.n_counts[bp];
                total += self.not_n_counts[bp];
            }

            self.percentages[i] = if total > 0 {
                100.0 * (n_count as f64 / total as f64)
            } else {
                0.0
            };
        }

        self.calculated = true;
    }
}

impl QCModule for NContent {
    fn name(&self) -> &str {
        "Per base N content"
    }

    fn description(&self) -> &str {
        "Shows the percentage of bases at each position which are not being called"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("n_content", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;
        let seq = sequence.sequence.as_bytes();

        if self.n_counts.len() < seq.len() {
            self.n_counts.resize(seq.len(), 0);
            self.not_n_counts.resize(seq.len(), 0);
        }

        for (i, &ch) in seq.iter().enumerate() {
            if ch == b'N' {
                self.n_counts[i] += 1;
            } else {
                self.not_n_counts[i] += 1;
            }
        }
    }

    fn reset(&mut self) {
        self.n_counts.clear();
        self.not_n_counts.clear();
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }
        let threshold = self.config.get_param("n_content", "error");
        self.percentages.iter().any(|&p| p > threshold)
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }
        let threshold = self.config.get_param("n_content", "warn");
        self.percentages.iter().any(|&p| p > threshold)
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.get_percentages();
        }

        buf.push_str("#Base\tN-Count\n");
        for i in 0..self.x_categories.len() {
            buf.push_str(&self.x_categories[i]);
            buf.push('\t');
            buf.push_str(&format_java_double(self.percentages[i]));
            buf.push('\n');
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.get_percentages();
        }
        if self.percentages.is_empty() {
            return None;
        }

        Some(ChartData::LineGraph {
            series: vec![Series {
                name: "%N".to_string(),
                data: self.percentages.clone(),
            }],
            x_labels: self.x_categories.clone(),
            x_axis_label: "Position in read (bp)".to_string(),
            title: "N content across all bases".to_string(),
            min_y: 0.0,
            max_y: 100.0,
        })
    }
}
