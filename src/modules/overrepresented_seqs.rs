use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use crate::config::FastQCConfig;
use crate::contaminant::ContaminantFinder;
use crate::modules::module_config::ModuleConfig;
use crate::modules::QCModule;
use crate::sequence::Sequence;

/// Shared state between OverRepresentedSeqs and DuplicationLevel modules.
pub struct SharedDuplicationData {
    pub sequences: HashMap<String, u64>,
    pub count: u64,
    pub unique_sequence_count: usize,
    pub count_at_unique_limit: u64,
}

impl Default for SharedDuplicationData {
    fn default() -> Self {
        Self::new()
    }
}

impl SharedDuplicationData {
    /// Construct empty shared data.
    pub fn new() -> Self {
        Self {
            sequences: HashMap::new(),
            count: 0,
            unique_sequence_count: 0,
            count_at_unique_limit: 0,
        }
    }
}

/// A single overrepresented sequence with its frequency and best contaminant match.
struct OverrepresentedSeq {
    seq: String,
    count: u64,
    percentage: f64,
    contaminant_hit: String,
}

/// Records unique-sequence counts (up to OBSERVATION_CUTOFF distinct sequences) and
/// produces a list of sequences appearing in more than a configurable fraction of reads.
pub struct OverRepresentedSeqs {
    shared: Arc<RwLock<SharedDuplicationData>>,
    overrepresented_seqs: Vec<OverrepresentedSeq>,
    calculated: bool,
    frozen: bool,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
    contaminant_finder: ContaminantFinder,
}

const OBSERVATION_CUTOFF: usize = 100000;

impl OverRepresentedSeqs {
    /// Construct the module wired to its shared sequence-count table and a contaminant DB.
    pub fn new(
        config: ModuleConfig,
        fqc_config: FastQCConfig,
        shared: Arc<RwLock<SharedDuplicationData>>,
    ) -> Self {
        let contaminant_finder = ContaminantFinder::new(fqc_config.contaminant_file.as_deref());

        Self {
            shared,
            overrepresented_seqs: Vec::new(),
            calculated: false,
            frozen: false,
            config,
            fqc_config,
            contaminant_finder,
        }
    }

    /// Filter the shared sequence table down to entries above the warn threshold and
    /// annotate each with its best contaminant match.
    fn get_overrepresented_seqs(&mut self) {
        let data = self.shared.read().unwrap();
        let count = data.count;

        let mut keepers = Vec::new();
        let warn_threshold = self.config.get_param("overrepresented", "warn");

        for (seq, &seq_count) in &data.sequences {
            let percentage = (seq_count as f64 / count as f64) * 100.0;
            if percentage > warn_threshold {
                let hit = self.contaminant_finder.find_contaminant_hit(seq);
                let hit_str = match hit {
                    Some(h) => h.to_string(),
                    None => "No Hit".to_string(),
                };
                keepers.push(OverrepresentedSeq {
                    seq: seq.clone(),
                    count: seq_count,
                    percentage,
                    contaminant_hit: hit_str,
                });
            }
        }

        // Sort by count descending
        keepers.sort_by(|a, b| b.count.cmp(&a.count));

        self.overrepresented_seqs = keepers;
        self.calculated = true;
    }
}

impl QCModule for OverRepresentedSeqs {
    fn name(&self) -> &str {
        "Overrepresented sequences"
    }

    fn description(&self) -> &str {
        "Identifies sequences which are overrepresented in the set"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.config.get_param("overrepresented", "ignore") > 0.0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        self.calculated = false;

        let mut data = self.shared.write().unwrap();
        data.count += 1;

        let mut seq = sequence.sequence.clone();

        if self.fqc_config.dup_length != 0 && seq.len() > self.fqc_config.dup_length {
            seq.truncate(self.fqc_config.dup_length);
        } else if self.fqc_config.dup_length == 0 && seq.len() > 50 {
            seq.truncate(50);
        }

        if let Some(entry) = data.sequences.get_mut(&seq) {
            *entry += 1;
            if !self.frozen {
                data.count_at_unique_limit = data.count;
            }
        } else if !self.frozen {
            data.sequences.insert(seq, 1);
            data.unique_sequence_count += 1;
            data.count_at_unique_limit = data.count;
            if data.unique_sequence_count == OBSERVATION_CUTOFF {
                self.frozen = true;
            }
        }
    }

    fn reset(&mut self) {
        let mut data = self.shared.write().unwrap();
        data.count = 0;
        data.sequences.clear();
        data.unique_sequence_count = 0;
        data.count_at_unique_limit = 0;
        self.frozen = false;
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.get_overrepresented_seqs();
        }
        if !self.overrepresented_seqs.is_empty()
            && self.overrepresented_seqs[0].percentage
                > self.config.get_param("overrepresented", "error")
        {
            return true;
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.get_overrepresented_seqs();
        }
        !self.overrepresented_seqs.is_empty()
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.get_overrepresented_seqs();
        }

        if self.overrepresented_seqs.is_empty() {
            // Java: no header line when there are no overrepresented sequences
        } else {
            buf.push_str("#Sequence\tCount\tPercentage\tPossible Source\n");
            for os in &self.overrepresented_seqs {
                buf.push_str(&os.seq);
                buf.push('\t');
                buf.push_str(&os.count.to_string());
                buf.push('\t');
                // Java: Math.round(percentage * 100.0) / 100.0
                let rounded = (os.percentage * 100.0).round() / 100.0;
                buf.push_str(&format_java_rounded(rounded));
                buf.push('\t');
                buf.push_str(&os.contaminant_hit);
                buf.push('\n');
            }
        }
    }
}

/// Format a rounded percentage like Java does.
/// Java's Math.round(x*100.0)/100.0 then toString() gives e.g. "100.0", "0.12"
fn format_java_rounded(v: f64) -> String {
    if v.fract() == 0.0 && v.is_finite() {
        format!("{:.1}", v)
    } else {
        format!("{}", v)
    }
}
