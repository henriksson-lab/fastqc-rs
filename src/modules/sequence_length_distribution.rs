use crate::config::FastQCConfig;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;

pub struct SequenceLengthDistribution {
    length_counts: Vec<u64>,
    graph_counts: Vec<f64>,
    x_categories: Vec<String>,
    calculated: bool,
    config: ModuleConfig,
    nogroup: bool,
}

impl SequenceLengthDistribution {
    pub fn new(config: ModuleConfig, fqc_config: &FastQCConfig) -> Self {
        Self {
            length_counts: Vec::new(),
            graph_counts: Vec::new(),
            x_categories: Vec::new(),
            calculated: false,
            config,
            nogroup: fqc_config.nogroup,
        }
    }

    /// Get size distribution parameters [start, interval].
    /// Exact port of Java SequenceLengthDistribution.getSizeDistribution().
    fn get_size_distribution(&self, min: i32, max: i32) -> [i32; 2] {
        if self.nogroup {
            return [min, 1];
        }

        let mut base: i32 = 1;
        while base > (max - min) {
            base /= 10;
        }

        let interval;
        let divisions = [1, 2, 5];

        'outer: loop {
            for &d in &divisions {
                let tester = base * d;
                if (max - min) / tester <= 50 {
                    interval = tester;
                    break 'outer;
                }
            }
            base *= 10;
        }

        let basic_division = min / interval;
        let starting = basic_division * interval;

        [starting, interval]
    }

    fn calculate_distribution(&mut self) {
        let mut max_len: i32 = 0;
        let mut min_len: i32 = -1;

        for (i, &count) in self.length_counts.iter().enumerate() {
            if count > 0 {
                if min_len < 0 {
                    min_len = i as i32;
                }
                max_len = i as i32;
            }
        }

        if min_len < 0 {
            min_len = 0;
        }

        // Add padding
        if min_len > 0 {
            min_len -= 1;
        }
        max_len += 1;

        let start_and_interval = self.get_size_distribution(min_len, max_len);
        let start = start_and_interval[0];
        let interval = start_and_interval[1];

        // Count categories
        let mut categories = 0;
        let mut current_value = start;
        while current_value <= max_len {
            categories += 1;
            current_value += interval;
        }

        self.graph_counts = vec![0.0; categories];
        self.x_categories = Vec::with_capacity(categories);

        for i in 0..categories {
            let min_value = start + interval * i as i32;
            let mut max_value = start + interval * (i as i32 + 1) - 1;
            if max_value > max_len {
                max_value = max_len;
            }

            for bp in min_value..=max_value {
                if (bp as usize) < self.length_counts.len() {
                    self.graph_counts[i] += self.length_counts[bp as usize] as f64;
                }
            }

            if interval == 1 {
                self.x_categories.push(min_value.to_string());
            } else {
                self.x_categories
                    .push(format!("{}-{}", min_value, max_value));
            }
        }

        self.calculated = true;
    }
}

impl QCModule for SequenceLengthDistribution {
    fn name(&self) -> &str {
        "Sequence Length Distribution"
    }

    fn description(&self) -> &str {
        "Shows the distribution of sequence length over all sequences"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("sequence_length", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;
        let seq_len = sequence.sequence.len();

        if seq_len + 2 > self.length_counts.len() {
            self.length_counts.resize(seq_len + 2, 0);
        }

        self.length_counts[seq_len] += 1;
    }

    fn reset(&mut self) {
        self.length_counts.clear();
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.calculate_distribution();
        }
        if self.config.get_param("sequence_length", "error") == 0.0 {
            return false;
        }
        if !self.length_counts.is_empty() && self.length_counts[0] > 0 {
            return true;
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.calculate_distribution();
        }
        if self.config.get_param("sequence_length", "warn") == 0.0 {
            return false;
        }
        let mut seen_length = false;
        for &count in &self.length_counts {
            if count > 0 {
                if seen_length {
                    return true;
                }
                seen_length = true;
            }
        }
        false
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.calculate_distribution();
        }

        buf.push_str("#Length\tCount\n");
        for i in 0..self.x_categories.len() {
            if (i == 0 || i == self.x_categories.len() - 1) && self.graph_counts[i] == 0.0 {
                continue;
            }
            buf.push_str(&self.x_categories[i]);
            buf.push('\t');
            buf.push_str(&format_java_double(self.graph_counts[i]));
            buf.push('\n');
        }
    }
}
