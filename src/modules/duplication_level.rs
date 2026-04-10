use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use crate::charts::{ChartData, Series};
use crate::modules::module_config::ModuleConfig;
use crate::modules::overrepresented_seqs::SharedDuplicationData;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;

pub struct DuplicationLevel {
    shared: Arc<RwLock<SharedDuplicationData>>,
    total_percentages: Option<Vec<f64>>,
    percent_different_seqs: f64,
    labels: Vec<String>,
    config: ModuleConfig,
}

impl DuplicationLevel {
    pub fn new(config: ModuleConfig, shared: Arc<RwLock<SharedDuplicationData>>) -> Self {
        Self {
            shared,
            total_percentages: None,
            percent_different_seqs: 0.0,
            labels: Vec::new(),
            config,
        }
    }

    pub fn calculate_levels(&mut self) {
        if self.total_percentages.is_some() {
            return;
        }

        let mut total_percentages = vec![0.0_f64; 16];

        let data = self.shared.read().unwrap();

        if data.count == 0 || data.sequences.is_empty() {
            let mut tp = vec![0.0; 16];
            tp[0] = 100.0;
            self.total_percentages = Some(tp);
            self.percent_different_seqs = 100.0;
            self.labels = vec![
                "1".into(), "2".into(), "3".into(), "4".into(), "5".into(),
                "6".into(), "7".into(), "8".into(), "9".into(), ">10".into(),
                ">50".into(), ">100".into(), ">500".into(), ">1k".into(),
                ">5k".into(), ">10k+".into(),
            ];
            return;
        }

        // Collate counts: how many sequences have each duplication level
        let mut collated_counts: HashMap<u64, u64> = HashMap::new();
        for &this_count in data.sequences.values() {
            *collated_counts.entry(this_count).or_insert(0) += 1;
        }

        // Correct each count
        let mut corrected_counts: HashMap<u64, f64> = HashMap::new();
        for (&dup_level, &count) in &collated_counts {
            corrected_counts.insert(
                dup_level,
                get_corrected_count(
                    data.count_at_unique_limit,
                    data.count,
                    dup_level,
                    count,
                ),
            );
        }

        // Calculate raw and deduplicated proportions
        let mut dedup_total = 0.0_f64;
        let mut raw_total = 0.0_f64;

        for (&dup_level, &count) in &corrected_counts {
            dedup_total += count;
            raw_total += count * dup_level as f64;

            let temp_dup_slot = dup_level as i64 - 1;

            let dup_slot: usize = if !(0..=9999).contains(&temp_dup_slot) {
                15
            } else if temp_dup_slot > 4999 {
                14
            } else if temp_dup_slot > 999 {
                13
            } else if temp_dup_slot > 499 {
                12
            } else if temp_dup_slot > 99 {
                11
            } else if temp_dup_slot > 49 {
                10
            } else if temp_dup_slot > 9 {
                9
            } else {
                temp_dup_slot as usize
            };

            total_percentages[dup_slot] += count * dup_level as f64;
        }

        self.labels = vec![String::new(); 16];
        for i in 0..total_percentages.len() {
            if i < 9 {
                self.labels[i] = format!("{}", i + 1);
            } else if i == 9 {
                self.labels[i] = ">10".to_string();
            } else if i == 10 {
                self.labels[i] = ">50".to_string();
            } else if i == 11 {
                self.labels[i] = ">100".to_string();
            } else if i == 12 {
                self.labels[i] = ">500".to_string();
            } else if i == 13 {
                self.labels[i] = ">1k".to_string();
            } else if i == 14 {
                self.labels[i] = ">5k".to_string();
            } else if i == 15 {
                self.labels[i] = ">10k".to_string();
            }

            if raw_total > 0.0 {
                total_percentages[i] /= raw_total;
                total_percentages[i] *= 100.0;
            }
        }

        self.percent_different_seqs = if raw_total == 0.0 {
            100.0
        } else {
            (dedup_total / raw_total) * 100.0
        };

        self.total_percentages = Some(total_percentages);
    }
}

fn get_corrected_count(
    count_at_limit: u64,
    total_count: u64,
    duplication_level: u64,
    number_of_observations: u64,
) -> f64 {
    // Bail out early
    if count_at_limit == total_count {
        return number_of_observations as f64;
    }

    if total_count - number_of_observations < count_at_limit {
        return number_of_observations as f64;
    }

    let mut p_not_seeing_at_limit = 1.0_f64;

    let limit_of_caring =
        1.0 - (number_of_observations as f64 / (number_of_observations as f64 + 0.01));

    for i in 0..count_at_limit {
        p_not_seeing_at_limit *= ((total_count - i) as f64 - duplication_level as f64)
            / (total_count - i) as f64;

        if p_not_seeing_at_limit < limit_of_caring {
            p_not_seeing_at_limit = 0.0;
            break;
        }
    }

    let p_seeing_at_limit = 1.0 - p_not_seeing_at_limit;

    if p_seeing_at_limit > 0.0 {
        number_of_observations as f64 / p_seeing_at_limit
    } else {
        number_of_observations as f64
    }
}

impl QCModule for DuplicationLevel {
    fn name(&self) -> &str {
        "Sequence Duplication Levels"
    }

    fn description(&self) -> &str {
        "Plots the number of sequences which are duplicated to different levels"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        self.config.get_param("duplication", "ignore") > 0.0
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("duplication", "ignore") > 0.0
    }

    fn process_sequence(&mut self, _sequence: &Sequence) {
        // No-op: uses data from OverRepresentedSeqs
    }

    fn reset(&mut self) {
        self.total_percentages = None;
    }

    fn raises_error(&mut self) -> bool {
        self.calculate_levels();
        self.percent_different_seqs < self.config.get_param("duplication", "error")
    }

    fn raises_warning(&mut self) -> bool {
        self.calculate_levels();
        self.percent_different_seqs < self.config.get_param("duplication", "warn")
    }

    fn make_data_report(&mut self, buf: &mut String) {
        self.calculate_levels();

        let total_percentages = self.total_percentages.as_ref().unwrap();

        buf.push_str("#Total Deduplicated Percentage\t");
        buf.push_str(&format_java_double(self.percent_different_seqs));
        buf.push('\n');

        buf.push_str("#Duplication Level\tPercentage of total\n");
        for i in 0..self.labels.len() {
            buf.push_str(&self.labels[i]);
            if i == self.labels.len() - 1 {
                buf.push('+');
            }
            buf.push('\t');
            buf.push_str(&format_java_double(total_percentages[i]));
            buf.push('\n');
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        self.calculate_levels();
        let total_percentages = self.total_percentages.as_ref()?;

        // Build labels with "+" suffix on the last one
        let mut x_labels = self.labels.clone();
        if let Some(last) = x_labels.last_mut() {
            last.push('+');
        }

        Some(ChartData::LineGraph {
            series: vec![Series {
                name: "% Total sequences".to_string(),
                data: total_percentages.clone(),
            }],
            x_labels,
            x_axis_label: "Sequence Duplication Level".to_string(),
            title: "Percent of seqs remaining if deduplicated".to_string(),
            min_y: 0.0,
            max_y: 100.0,
        })
    }
}
