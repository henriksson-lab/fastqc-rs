use crate::charts::{ChartData, Series};
use crate::config::FastQCConfig;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::utils::base_group::BaseGroup;

const DEFAULT_ADAPTERS: &str = include_str!("../resources/adapter_list.txt");

struct Adapter {
    sequence: String,
    positions: Vec<u64>,
}

impl Adapter {
    fn new(sequence: &str) -> Self {
        Self {
            sequence: sequence.to_string(),
            positions: vec![0; 1],
        }
    }

    fn expand_length_to(&mut self, new_length: usize) {
        if new_length <= self.positions.len() {
            return;
        }
        let old_len = self.positions.len();
        let last_val = if old_len > 0 {
            self.positions[old_len - 1]
        } else {
            0
        };
        self.positions.resize(new_length, last_val);
    }

    fn increment_count(&mut self, position: usize) {
        self.positions[position] += 1;
    }

    fn reset(&mut self) {
        self.positions.clear();
    }
}

pub struct AdapterContent {
    adapters: Vec<Adapter>,
    labels: Vec<String>,
    longest_sequence: usize,
    longest_adapter: usize,
    total_count: u64,
    calculated: bool,
    enrichments: Vec<Vec<f64>>,
    x_labels: Vec<String>,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
}

impl AdapterContent {
    pub fn new(config: ModuleConfig, fqc_config: FastQCConfig) -> Self {
        let mut adapters = Vec::new();
        let mut labels = Vec::new();
        let mut longest_adapter = 0usize;

        let text = if let Some(ref adapter_file) = fqc_config.adapter_file {
            std::fs::read_to_string(adapter_file).unwrap_or_else(|e| {
                eprintln!("Warning: could not read adapter file {:?}: {}", adapter_file, e);
                DEFAULT_ADAPTERS.to_string()
            })
        } else {
            DEFAULT_ADAPTERS.to_string()
        };

        for line in text.lines() {
            let line = line.trim();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let sections: Vec<&str> = line.split('\t').filter(|s| !s.is_empty()).collect();
            if sections.len() != 2 {
                eprintln!(
                    "Expected 2 sections for contaminant line but got {} from {}",
                    sections.len(),
                    line
                );
                continue;
            }

            let name = sections[0].trim();
            let seq = sections[1].trim();
            if seq.len() > longest_adapter {
                longest_adapter = seq.len();
            }
            labels.push(name.to_string());
            adapters.push(Adapter::new(seq));
        }

        Self {
            adapters,
            labels,
            longest_sequence: 0,
            longest_adapter,
            total_count: 0,
            calculated: false,
            enrichments: Vec::new(),
            x_labels: Vec::new(),
            config,
            fqc_config,
        }
    }

    fn calculate_enrichment(&mut self) {
        if self.total_count == 0 {
            self.calculated = true;
            return;
        }
        let mut max_length = 0usize;
        for adapter in &self.adapters {
            if adapter.positions.len() > max_length {
                max_length = adapter.positions.len();
            }
        }

        let groups = BaseGroup::make_base_groups(max_length, &self.fqc_config);

        self.x_labels = Vec::with_capacity(groups.len());
        for group in &groups {
            self.x_labels.push(group.to_string());
        }

        self.enrichments = vec![vec![0.0; groups.len()]; self.adapters.len()];

        for (a, adapter) in self.adapters.iter().enumerate() {
            let positions = &adapter.positions;

            for (g, group) in groups.iter().enumerate() {
                let mut p = group.lower_count - 1;
                while p < group.upper_count && p < positions.len() {
                    self.enrichments[a][g] += (positions[p] as f64 * 100.0) / self.total_count as f64;
                    p += 1;
                }

                self.enrichments[a][g] /=
                    (group.upper_count - group.lower_count + 1) as f64;
            }
        }

        self.calculated = true;
    }
}

impl QCModule for AdapterContent {
    fn name(&self) -> &str {
        "Adapter Content"
    }

    fn description(&self) -> &str {
        "Searches for specific adapter sequences in a library"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("adapter", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;
        self.total_count += 1;

        let seq = &sequence.sequence;

        if seq.len() > self.longest_sequence
            && seq.len() as i64 - self.longest_adapter as i64 > 0
        {
            self.longest_sequence = seq.len();
            for adapter in self.adapters.iter_mut() {
                adapter.expand_length_to((self.longest_sequence - self.longest_adapter) + 1);
            }
        }

        for adapter in self.adapters.iter_mut() {
            if let Some(index) = seq.find(&adapter.sequence) {
                let max_pos = self.longest_sequence.saturating_sub(self.longest_adapter);
                for i in index..=max_pos {
                    if i < adapter.positions.len() {
                        adapter.increment_count(i);
                    }
                }
            }
        }
    }

    fn reset(&mut self) {
        self.calculated = false;
        self.total_count = 0;
        self.longest_sequence = 0;
        for adapter in self.adapters.iter_mut() {
            adapter.reset();
        }
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.calculate_enrichment();
        }

        for enrichment in &self.enrichments {
            for &val in enrichment {
                if val > self.config.get_param("adapter", "error") {
                    return true;
                }
            }
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        if self.longest_adapter > self.longest_sequence {
            return true;
        }

        if !self.calculated {
            self.calculate_enrichment();
        }

        for enrichment in &self.enrichments {
            for &val in enrichment {
                if val > self.config.get_param("adapter", "warn") {
                    return true;
                }
            }
        }
        false
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if self.longest_adapter > self.longest_sequence {
            // Can't analyse adapters as read length is too short
            return;
        }

        if !self.calculated {
            self.calculate_enrichment();
        }

        // Header: #Position\tAdapter1\tAdapter2\t...
        buf.push_str("#Position");
        for label in &self.labels {
            buf.push('\t');
            buf.push_str(label);
        }
        buf.push('\n');

        // Data rows
        if !self.enrichments.is_empty() {
            for r in 0..self.enrichments[0].len() {
                buf.push_str(&self.x_labels[r]);
                for a in 0..self.adapters.len() {
                    buf.push('\t');
                    buf.push_str(&format_java_double(self.enrichments[a][r]));
                }
                buf.push('\n');
            }
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.calculate_enrichment();
        }
        if self.enrichments.is_empty() || self.x_labels.is_empty() {
            return None;
        }

        let series: Vec<Series> = self.labels.iter().zip(self.enrichments.iter())
            .map(|(name, data)| Series {
                name: name.clone(),
                data: data.clone(),
            })
            .collect();

        Some(ChartData::LineGraph {
            series,
            x_labels: self.x_labels.clone(),
            x_axis_label: "Position in read (bp)".to_string(),
            title: "% Adapter".to_string(),
            min_y: 0.0,
            max_y: 100.0,
        })
    }
}
