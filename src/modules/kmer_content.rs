use std::collections::HashMap;

use crate::charts::{ChartData, Series};
use crate::config::FastQCConfig;
use crate::modules::module_config::ModuleConfig;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::utils::base_group::BaseGroup;

use statrs::distribution::{Binomial, DiscreteCDF};

struct Kmer {
    sequence: String,
    count: u64,
    lowest_p_value: f32,
    obs_exp_positions: Option<Vec<f32>>,
    positions: Vec<u64>,
}

impl Kmer {
    fn new(sequence: &str, position: usize, seq_length: usize) -> Self {
        let mut positions = vec![0u64; seq_length];
        positions[position] += 1;
        Self {
            sequence: sequence.to_string(),
            count: 1,
            lowest_p_value: 0.0,
            obs_exp_positions: None,
            positions,
        }
    }

    fn increment_count(&mut self, position: usize) {
        self.count += 1;

        if position >= self.positions.len() {
            self.positions.resize(position + 1, 0);
        }
        self.positions[position] += 1;
    }
}

pub struct KmerContent {
    kmers: HashMap<String, Kmer>,
    longest_sequence: usize,
    total_kmer_counts: Vec<Vec<u64>>,
    skip_count: u64,
    min_kmer_size: usize,
    max_kmer_size: usize,
    calculated: bool,
    enriched_kmers: Option<Vec<EnrichedKmer>>,
    x_labels: Vec<String>,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
}

/// Flattened kmer data after calculation, ready for reporting.
struct EnrichedKmer {
    sequence: String,
    count: u64,
    p_value: f32,
    max_obs_exp: f32,
    max_position_label: String,
    obs_exp_positions: Vec<f32>,
}

impl KmerContent {
    pub fn new(config: ModuleConfig, fqc_config: FastQCConfig) -> Self {
        let kmer_size = fqc_config.kmer_size;
        Self {
            kmers: HashMap::with_capacity(4usize.pow(kmer_size as u32)),
            longest_sequence: 0,
            total_kmer_counts: Vec::new(),
            skip_count: 0,
            min_kmer_size: kmer_size,
            max_kmer_size: kmer_size,
            calculated: false,
            enriched_kmers: None,
            x_labels: Vec::new(),
            config,
            fqc_config,
        }
    }

    fn add_kmer_count(&mut self, position: usize, kmer_length: usize, kmer: &str) {
        if position >= self.total_kmer_counts.len() {
            let old_len = self.total_kmer_counts.len();
            self.total_kmer_counts
                .resize_with(position + 1, || vec![0u64; self.max_kmer_size]);
            // Ensure existing entries have the right length
            for i in old_len..self.total_kmer_counts.len() {
                if self.total_kmer_counts[i].len() < self.max_kmer_size {
                    self.total_kmer_counts[i].resize(self.max_kmer_size, 0);
                }
            }
        }

        if kmer.contains('N') {
            return;
        }

        self.total_kmer_counts[position][kmer_length - 1] += 1;
    }

    fn calculate_enrichment(&mut self) {
        let group_length = if self.longest_sequence >= self.min_kmer_size {
            (self.longest_sequence - self.min_kmer_size) + 1
        } else {
            0
        };

        if group_length == 0 {
            self.enriched_kmers = Some(Vec::new());
            self.x_labels = Vec::new();
            self.calculated = true;
            return;
        }

        let groups = BaseGroup::make_base_groups(group_length, &self.fqc_config);
        self.x_labels = groups.iter().map(|g| g.to_string()).collect();

        let mut uneven_kmers: Vec<(String, u64, f32, Vec<f32>)> = Vec::new();

        for kmer in self.kmers.values_mut() {
            let kmer_len = kmer.sequence.len();

            let mut total_kmer_count: u64 = 0;
            for i in 0..self.total_kmer_counts.len() {
                total_kmer_count += self.total_kmer_counts[i][kmer_len - 1];
            }

            let expected_proportion = kmer.count as f64 / total_kmer_count as f64;

            let mut obs_exp_positions = vec![0.0_f32; groups.len()];
            let mut binomial_p_values = vec![0.0_f32; groups.len()];

            let position_counts = &kmer.positions;

            for (g, group) in groups.iter().enumerate() {
                let mut total_group_count: u64 = 0;
                let mut total_group_hits: u64 = 0;

                let mut p = group.lower_count - 1;
                while p < group.upper_count && p < position_counts.len() {
                    total_group_count += self.total_kmer_counts[p][kmer_len - 1];
                    total_group_hits += position_counts[p];
                    p += 1;
                }

                let predicted = expected_proportion * total_group_count as f64;
                obs_exp_positions[g] = if predicted > 0.0 {
                    (total_group_hits as f64 / predicted) as f32
                } else {
                    0.0
                };

                // Binomial test
                if total_group_count > 0 {
                    let bd = Binomial::new(expected_proportion, total_group_count).unwrap();
                    if total_group_hits as f64 > predicted {
                        let cdf = bd.cdf(total_group_hits);
                        binomial_p_values[g] = ((1.0 - cdf) * 4.0_f64.powi(kmer_len as i32)) as f32;
                    } else {
                        binomial_p_values[g] = 1.0;
                    }
                } else {
                    binomial_p_values[g] = 1.0;
                }
            }

            kmer.obs_exp_positions = Some(obs_exp_positions.clone());

            let mut lowest_p_value: f32 = 1.0;
            for i in 0..binomial_p_values.len() {
                if binomial_p_values[i] < 0.01
                    && obs_exp_positions[i] > 5.0
                    && binomial_p_values[i] < lowest_p_value
                {
                    lowest_p_value = binomial_p_values[i];
                }
            }

            if lowest_p_value < 0.01 {
                kmer.lowest_p_value = lowest_p_value;
                uneven_kmers.push((
                    kmer.sequence.clone(),
                    kmer.count,
                    lowest_p_value,
                    obs_exp_positions,
                ));
            }
        }

        // Sort by max obs/exp descending
        // Need to compute max_obs_exp for sorting
        uneven_kmers.sort_by(|a, b| {
            let a_max = a.3.iter().cloned().fold(0.0_f32, f32::max);
            let b_max = b.3.iter().cloned().fold(0.0_f32, f32::max);
            b_max
                .partial_cmp(&a_max)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Truncate to top 20
        if uneven_kmers.len() > 20 {
            uneven_kmers.truncate(20);
        }

        // Build enriched kmers for reporting
        let mut enriched = Vec::with_capacity(uneven_kmers.len());
        for (seq, count, p_value, obs_exp) in &uneven_kmers {
            let max_obs_exp = obs_exp.iter().cloned().fold(0.0_f32, f32::max);

            // Find max position
            let mut max_val = 0.0_f32;
            let mut max_pos = 0usize;
            for (i, &v) in obs_exp.iter().enumerate() {
                if v > max_val {
                    max_val = v;
                    max_pos = i;
                }
            }

            let max_position_label = if max_pos < groups.len() {
                groups[max_pos].to_string()
            } else {
                "1".to_string()
            };

            enriched.push(EnrichedKmer {
                sequence: seq.clone(),
                count: *count,
                p_value: *p_value,
                max_obs_exp,
                max_position_label,
                obs_exp_positions: obs_exp.clone(),
            });
        }

        self.enriched_kmers = Some(enriched);
        self.kmers.clear();
        self.calculated = true;
    }
}

impl QCModule for KmerContent {
    fn name(&self) -> &str {
        "Kmer Content"
    }

    fn description(&self) -> &str {
        "Identifies short sequences which have uneven representation"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("kmer", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;

        self.skip_count += 1;
        if !self.skip_count.is_multiple_of(50) {
            return;
        }

        let seq = if sequence.sequence.len() > 500 {
            &sequence.sequence[..500]
        } else {
            &sequence.sequence
        };

        if seq.len() > self.longest_sequence {
            self.longest_sequence = seq.len();
        }

        for kmer_size in self.min_kmer_size..=self.max_kmer_size {
            if seq.len() < kmer_size {
                continue;
            }
            for i in 0..=(seq.len() - kmer_size) {
                let kmer = &seq[i..i + kmer_size];

                self.add_kmer_count(i, kmer_size, kmer);

                if kmer.contains('N') {
                    continue;
                }

                if let Some(existing) = self.kmers.get_mut(kmer) {
                    existing.increment_count(i);
                } else {
                    let seq_length = (seq.len() - kmer_size) + 1;
                    self.kmers
                        .insert(kmer.to_string(), Kmer::new(kmer, i, seq_length));
                }
            }
        }
    }

    fn reset(&mut self) {
        self.calculated = false;
        self.total_kmer_counts.clear();
        self.longest_sequence = 0;
        self.skip_count = 0;
        self.enriched_kmers = None;
        self.x_labels.clear();
        self.kmers.clear();
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.calculate_enrichment();
        }

        let enriched = self.enriched_kmers.as_ref().unwrap();
        if !enriched.is_empty() {
            let log_p = 0.0 - (enriched[0].p_value as f64).log10();
            return log_p > self.config.get_param("kmer", "error");
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.calculate_enrichment();
        }

        let enriched = self.enriched_kmers.as_ref().unwrap();
        if !enriched.is_empty() {
            let log_p = 0.0 - (enriched[0].p_value as f64).log10();
            return log_p > self.config.get_param("kmer", "warn");
        }
        false
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.calculate_enrichment();
        }

        let enriched = self.enriched_kmers.as_ref().unwrap();

        if !enriched.is_empty() {
            buf.push_str("#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position\n");
            for ek in enriched {
                buf.push_str(&ek.sequence);
                buf.push('\t');
                // Java: count * 5 (because 2% sampling = every 50th, so multiply by 5... wait)
                // Actually Java getValueAt case 1: kmers[rowIndex].count()*5
                buf.push_str(&(ek.count * 5).to_string());
                buf.push('\t');
                buf.push_str(&format_java_float(ek.p_value));
                buf.push('\t');
                buf.push_str(&format_java_float(ek.max_obs_exp));
                buf.push('\t');
                buf.push_str(&ek.max_position_label);
                buf.push('\n');
            }
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.calculate_enrichment();
        }

        let enriched = self.enriched_kmers.as_ref()?;
        if enriched.is_empty() || self.x_labels.is_empty() {
            return None;
        }

        // Take up to 6 top enriched kmers for the chart
        let top_n = enriched.len().min(6);
        let series: Vec<Series> = enriched[..top_n]
            .iter()
            .map(|ek| {
                // Convert obs/exp to log2 scale
                let data: Vec<f64> = ek
                    .obs_exp_positions
                    .iter()
                    .map(|&v| if v > 0.0 { (v as f64).log2() } else { 0.0 })
                    .collect();
                Series {
                    name: ek.sequence.clone(),
                    data,
                }
            })
            .collect();

        let max_y = series
            .iter()
            .flat_map(|s| s.data.iter())
            .cloned()
            .fold(0.0_f64, f64::max);
        let min_y = series
            .iter()
            .flat_map(|s| s.data.iter())
            .cloned()
            .fold(0.0_f64, f64::min);

        Some(ChartData::LineGraph {
            series,
            x_labels: self.x_labels.clone(),
            x_axis_label: "Position in read (bp)".to_string(),
            title: "Log2 Obs/Exp".to_string(),
            min_y,
            max_y,
        })
    }
}

/// Format f32 like Java Float.toString()
fn format_java_float(v: f32) -> String {
    // Java Float.toString() format
    format!("{}", v)
}
