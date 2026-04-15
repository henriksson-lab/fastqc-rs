use crate::charts::{ChartData, Series};
use crate::modules::gc_model::GCModel;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::statistics::NormalDistribution;

pub struct PerSequenceGCContent {
    gc_distribution: [f64; 101],
    theoretical_distribution: [f64; 101],
    x_categories: Vec<usize>,
    deviation_percent: f64,
    calculated: bool,
    cached_models: Vec<Option<GCModel>>,
    config: ModuleConfig,
}

impl PerSequenceGCContent {
    pub fn new(config: ModuleConfig) -> Self {
        Self {
            gc_distribution: [0.0; 101],
            theoretical_distribution: [0.0; 101],
            x_categories: Vec::new(),
            deviation_percent: 0.0,
            calculated: false,
            cached_models: vec![None; 200],
            config,
        }
    }

    fn truncate_sequence<'a>(&self, seq: &'a str) -> &'a str {
        if seq.len() > 1000 {
            let length = (seq.len() / 1000) * 1000;
            &seq[..length]
        } else if seq.len() > 100 {
            let length = (seq.len() / 100) * 100;
            &seq[..length]
        } else {
            seq
        }
    }

    fn calculate_distribution(&mut self) {
        let mut max = 0.0_f64;
        self.x_categories = Vec::with_capacity(101);
        let mut total_count = 0.0_f64;

        let mut first_mode: usize = 0;
        let mut mode_count = 0.0_f64;

        for i in 0..self.gc_distribution.len() {
            self.x_categories.push(i);
            total_count += self.gc_distribution[i];

            if self.gc_distribution[i] > mode_count {
                mode_count = self.gc_distribution[i];
                first_mode = i;
            }
            if self.gc_distribution[i] > max {
                max = self.gc_distribution[i];
            }
        }

        // Average over adjacent points which stay above 90% of the modal value
        let mut mode = 0.0_f64;
        let mut mode_duplicates = 0;

        let mut fell_off_top = true;

        for i in first_mode..self.gc_distribution.len() {
            if self.gc_distribution[i]
                > self.gc_distribution[first_mode] - (self.gc_distribution[first_mode] / 10.0)
            {
                mode += i as f64;
                mode_duplicates += 1;
            } else {
                fell_off_top = false;
                break;
            }
        }

        let mut fell_off_bottom = true;

        if first_mode > 0 {
            for i in (0..first_mode).rev() {
                if self.gc_distribution[i]
                    > self.gc_distribution[first_mode] - (self.gc_distribution[first_mode] / 10.0)
                {
                    mode += i as f64;
                    mode_duplicates += 1;
                } else {
                    fell_off_bottom = false;
                    break;
                }
            }
        }

        if fell_off_bottom || fell_off_top {
            mode = first_mode as f64;
        } else {
            mode /= mode_duplicates as f64;
        }

        // Calculate standard deviation
        let mut stdev = 0.0_f64;

        for i in 0..self.gc_distribution.len() {
            stdev += (i as f64 - mode).powi(2) * self.gc_distribution[i];
        }

        stdev /= total_count - 1.0;
        stdev = stdev.sqrt();

        let nd = NormalDistribution::new(mode, stdev);

        self.deviation_percent = 0.0;

        for i in 0..self.theoretical_distribution.len() {
            let probability = nd.pdf(i as f64);
            self.theoretical_distribution[i] = probability * total_count;

            self.deviation_percent +=
                (self.theoretical_distribution[i] - self.gc_distribution[i]).abs();
        }

        self.deviation_percent /= total_count;
        self.deviation_percent *= 100.0;

        self.calculated = true;
    }
}

impl QCModule for PerSequenceGCContent {
    fn name(&self) -> &str {
        "Per sequence GC content"
    }

    fn description(&self) -> &str {
        "Shows the distribution of GC contents for whole sequences"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("gc_sequence", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        let seq = self.truncate_sequence(&sequence.sequence);

        if seq.is_empty() {
            return;
        }

        let seq_bytes = seq.as_bytes();
        let seq_len = seq_bytes.len();

        let mut gc_count = 0usize;
        for &b in seq_bytes {
            if b == b'G' || b == b'C' {
                gc_count += 1;
            }
        }

        // Extend cached_models if needed
        if seq_len >= self.cached_models.len() {
            self.cached_models.resize_with(seq_len + 1, || None);
        }

        if self.cached_models[seq_len].is_none() {
            self.cached_models[seq_len] = Some(GCModel::new(seq_len));
        }

        let values = self.cached_models[seq_len]
            .as_ref()
            .unwrap()
            .get_model_values(gc_count);

        for v in values {
            self.gc_distribution[v.percentage] += v.increment;
        }
    }

    fn reset(&mut self) {
        self.gc_distribution = [0.0; 101];
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.calculate_distribution();
        }
        self.deviation_percent > self.config.get_param("gc_sequence", "error")
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.calculate_distribution();
        }
        self.deviation_percent > self.config.get_param("gc_sequence", "warn")
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.calculate_distribution();
        }

        buf.push_str("#GC Content\tCount\n");
        for i in 0..self.x_categories.len() {
            buf.push_str(&self.x_categories[i].to_string());
            buf.push('\t');
            buf.push_str(&format_java_double(self.gc_distribution[i]));
            buf.push('\n');
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.calculate_distribution();
        }

        let x_labels: Vec<String> = (0..=100).map(|i| i.to_string()).collect();
        let max_y = self
            .gc_distribution
            .iter()
            .chain(self.theoretical_distribution.iter())
            .cloned()
            .fold(0.0_f64, f64::max);

        Some(ChartData::LineGraph {
            series: vec![
                Series {
                    name: "GC count per read".to_string(),
                    data: self.gc_distribution.to_vec(),
                },
                Series {
                    name: "Theoretical Distribution".to_string(),
                    data: self.theoretical_distribution.to_vec(),
                },
            ],
            x_labels,
            x_axis_label: "Mean GC content (%)".to_string(),
            title: "GC distribution over all sequences".to_string(),
            min_y: 0.0,
            max_y,
        })
    }
}
