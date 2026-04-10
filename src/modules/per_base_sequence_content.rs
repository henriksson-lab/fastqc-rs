use crate::charts::{ChartData, Series};
use crate::config::FastQCConfig;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::utils::base_group::BaseGroup;

pub struct PerBaseSequenceContent {
    g_counts: Vec<u64>,
    a_counts: Vec<u64>,
    c_counts: Vec<u64>,
    t_counts: Vec<u64>,
    /// Percentages in order: [t_percent, c_percent, a_percent, g_percent]
    percentages: Option<[Vec<f64>; 4]>,
    x_categories: Vec<String>,
    calculated: bool,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
}

impl PerBaseSequenceContent {
    pub fn new(config: ModuleConfig, fqc_config: FastQCConfig) -> Self {
        Self {
            g_counts: Vec::new(),
            a_counts: Vec::new(),
            c_counts: Vec::new(),
            t_counts: Vec::new(),
            percentages: None,
            x_categories: Vec::new(),
            calculated: false,
            config,
            fqc_config,
        }
    }

    fn get_percentages(&mut self) {
        if self.g_counts.is_empty() {
            self.percentages = Some([vec![], vec![], vec![], vec![]]);
            self.calculated = true;
            return;
        }
        let groups = BaseGroup::make_base_groups(self.g_counts.len(), &self.fqc_config);

        self.x_categories = Vec::with_capacity(groups.len());

        let mut g_percent = vec![0.0; groups.len()];
        let mut a_percent = vec![0.0; groups.len()];
        let mut t_percent = vec![0.0; groups.len()];
        let mut c_percent = vec![0.0; groups.len()];

        for (i, group) in groups.iter().enumerate() {
            self.x_categories.push(group.to_string());

            let mut g_count: u64 = 0;
            let mut a_count: u64 = 0;
            let mut t_count: u64 = 0;
            let mut c_count: u64 = 0;
            let mut total: u64 = 0;

            // Java: for (int bp=groups[i].lowerCount()-1;bp<groups[i].upperCount();bp++)
            for bp in (group.lower_count - 1)..group.upper_count {
                total += self.g_counts[bp];
                total += self.c_counts[bp];
                total += self.a_counts[bp];
                total += self.t_counts[bp];

                a_count += self.a_counts[bp];
                t_count += self.t_counts[bp];
                c_count += self.c_counts[bp];
                g_count += self.g_counts[bp];
            }

            if total > 0 {
                g_percent[i] = (g_count as f64 / total as f64) * 100.0;
                a_percent[i] = (a_count as f64 / total as f64) * 100.0;
                t_percent[i] = (t_count as f64 / total as f64) * 100.0;
                c_percent[i] = (c_count as f64 / total as f64) * 100.0;
            }
        }

        // Order: T, C, A, G (matching Java)
        self.percentages = Some([t_percent, c_percent, a_percent, g_percent]);
        self.calculated = true;
    }
}

impl QCModule for PerBaseSequenceContent {
    fn name(&self) -> &str {
        "Per base sequence content"
    }

    fn description(&self) -> &str {
        "Shows the relative amounts of each base at each position in a sequencing run"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("sequence", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;
        let seq = sequence.sequence.as_bytes();

        if self.g_counts.len() < seq.len() {
            self.g_counts.resize(seq.len(), 0);
            self.a_counts.resize(seq.len(), 0);
            self.c_counts.resize(seq.len(), 0);
            self.t_counts.resize(seq.len(), 0);
        }

        for (i, &ch) in seq.iter().enumerate() {
            match ch {
                b'G' => self.g_counts[i] += 1,
                b'A' => self.a_counts[i] += 1,
                b'T' => self.t_counts[i] += 1,
                b'C' => self.c_counts[i] += 1,
                _ => {}
            }
        }
    }

    fn reset(&mut self) {
        self.g_counts.clear();
        self.a_counts.clear();
        self.t_counts.clear();
        self.c_counts.clear();
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }

        let pct = self.percentages.as_ref().unwrap();
        // Percentages in order TCAG
        for i in 0..pct[0].len() {
            let gc_diff = (pct[1][i] - pct[3][i]).abs();
            let at_diff = (pct[0][i] - pct[2][i]).abs();
            if gc_diff > self.config.get_param("sequence", "error")
                || at_diff > self.config.get_param("sequence", "error")
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

        let pct = self.percentages.as_ref().unwrap();
        for i in 0..pct[0].len() {
            let gc_diff = (pct[1][i] - pct[3][i]).abs();
            let at_diff = (pct[0][i] - pct[2][i]).abs();
            if gc_diff > self.config.get_param("sequence", "warn")
                || at_diff > self.config.get_param("sequence", "warn")
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

        let pct = self.percentages.as_ref().unwrap();

        buf.push_str("#Base\tG\tA\tT\tC\n");
        for i in 0..self.x_categories.len() {
            buf.push_str(&self.x_categories[i]);
            buf.push('\t');
            // Java order in report: G, A, T, C  (indices 3, 2, 0, 1)
            buf.push_str(&format_java_double(pct[3][i]));
            buf.push('\t');
            buf.push_str(&format_java_double(pct[2][i]));
            buf.push('\t');
            buf.push_str(&format_java_double(pct[0][i]));
            buf.push('\t');
            buf.push_str(&format_java_double(pct[1][i]));
            buf.push('\n');
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.get_percentages();
        }
        let pct = self.percentages.as_ref()?;
        if pct[0].is_empty() {
            return None;
        }

        // Order: %T, %C, %A, %G (indices 0, 1, 2, 3 in the percentages array)
        Some(ChartData::LineGraph {
            series: vec![
                Series { name: "%T".to_string(), data: pct[0].clone() },
                Series { name: "%C".to_string(), data: pct[1].clone() },
                Series { name: "%A".to_string(), data: pct[2].clone() },
                Series { name: "%G".to_string(), data: pct[3].clone() },
            ],
            x_labels: self.x_categories.clone(),
            x_axis_label: "Position in read (bp)".to_string(),
            title: "Sequence content across all bases".to_string(),
            min_y: 0.0,
            max_y: 100.0,
        })
    }
}
