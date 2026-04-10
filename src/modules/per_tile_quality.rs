use std::collections::HashMap;

use crate::charts::ChartData;
use crate::config::FastQCConfig;
use crate::encoding::phred::PhredEncoding;
use crate::modules::module_config::ModuleConfig;
use crate::modules::per_sequence_quality::format_java_double;
use crate::modules::QCModule;
use crate::sequence::Sequence;
use crate::utils::base_group::BaseGroup;
use crate::utils::quality_count::QualityCount;

pub struct PerTileQuality {
    per_tile_quality_counts: HashMap<i32, Vec<QualityCount>>,
    current_length: usize,
    means: Vec<Vec<f64>>,
    x_labels: Vec<String>,
    tiles: Vec<i32>,
    max_deviation: f64,
    calculated: bool,
    total_count: u64,
    split_position: i32,
    ignore_in_report: bool,
    config: ModuleConfig,
    fqc_config: FastQCConfig,
}

impl PerTileQuality {
    pub fn new(config: ModuleConfig, fqc_config: FastQCConfig) -> Self {
        Self {
            per_tile_quality_counts: HashMap::new(),
            current_length: 0,
            means: Vec::new(),
            x_labels: Vec::new(),
            tiles: Vec::new(),
            max_deviation: 0.0,
            calculated: false,
            total_count: 0,
            split_position: -1,
            ignore_in_report: false,
            config,
            fqc_config,
        }
    }

    fn calculate_offsets(&self) -> (char, char) {
        let mut min_char: char = '\0';
        let mut max_char: char = '\0';

        for quality_counts in self.per_tile_quality_counts.values() {
            for qc in quality_counts {
                if min_char == '\0' {
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

        let groups = BaseGroup::make_base_groups(self.current_length, &self.fqc_config);

        let mut tile_numbers: Vec<i32> = self.per_tile_quality_counts.keys().copied().collect();
        tile_numbers.sort();

        self.tiles = tile_numbers.clone();

        self.means = vec![vec![0.0; groups.len()]; tile_numbers.len()];
        self.x_labels = Vec::with_capacity(groups.len());

        for (t, &tile_num) in tile_numbers.iter().enumerate() {
            for (i, group) in groups.iter().enumerate() {
                if t == 0 {
                    self.x_labels.push(group.to_string());
                }

                let min_base = group.lower_count;
                let max_base = group.upper_count;
                self.means[t][i] =
                    self.get_mean(tile_num, min_base, max_base, encoding.offset);
            }
        }

        // Normalize: subtract column averages
        let mut max_deviation = 0.0_f64;
        let num_tiles = tile_numbers.len();

        let mut average_qualities_per_group = vec![0.0_f64; groups.len()];

        for t in 0..num_tiles {
            for i in 0..groups.len() {
                average_qualities_per_group[i] += self.means[t][i];
            }
        }

        for avg in average_qualities_per_group.iter_mut() {
            *avg /= num_tiles as f64;
        }

        for i in 0..groups.len() {
            for t in 0..num_tiles {
                self.means[t][i] -= average_qualities_per_group[i];
                if self.means[t][i].abs() > max_deviation {
                    max_deviation = self.means[t][i].abs();
                }
            }
        }

        self.max_deviation = max_deviation;
        self.calculated = true;
    }

    fn get_mean(&self, tile: i32, min_bp: usize, max_bp: usize, offset: i32) -> f64 {
        let mut count = 0;
        let mut total = 0.0;

        let quality_counts = &self.per_tile_quality_counts[&tile];

        for i in (min_bp - 1)..max_bp {
            if quality_counts[i].get_total_count() > 0 {
                count += 1;
                total += quality_counts[i].get_mean(offset);
            }
        }

        if count > 0 {
            total / count as f64
        } else {
            0.0
        }
    }
}

impl QCModule for PerTileQuality {
    fn name(&self) -> &str {
        "Per tile sequence quality"
    }

    fn description(&self) -> &str {
        "Shows the perl tile Quality scores of all bases at a given position in a sequencing run"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore_in_report
            || self.config.get_param("tile", "ignore") > 0.0
            || self.current_length == 0
    }

    fn process_sequence(&mut self, sequence: &Sequence) {
        // Check if we can skip
        if self.total_count == 0 && self.config.get_param("tile", "ignore") > 0.0 {
            self.ignore_in_report = true;
        }

        if self.ignore_in_report {
            return;
        }

        if sequence.quality.is_empty() {
            return;
        }

        self.calculated = false;

        self.total_count += 1;
        if self.total_count > 10000 && !self.total_count.is_multiple_of(10) {
            return;
        }

        // Parse tile ID from sequence header
        let split_id: Vec<&str> = sequence.id.split(':').collect();

        let tile: i32;

        if self.split_position >= 0 {
            let pos = self.split_position as usize;
            if split_id.len() <= pos {
                self.ignore_in_report = true;
                return;
            }
            match split_id[pos].parse::<i32>() {
                Ok(v) => tile = v,
                Err(_) => {
                    self.ignore_in_report = true;
                    return;
                }
            }
        } else if split_id.len() >= 7 {
            self.split_position = 4;
            match split_id[4].parse::<i32>() {
                Ok(v) => tile = v,
                Err(_) => {
                    self.ignore_in_report = true;
                    return;
                }
            }
        } else if split_id.len() >= 5 {
            self.split_position = 2;
            match split_id[2].parse::<i32>() {
                Ok(v) => tile = v,
                Err(_) => {
                    self.ignore_in_report = true;
                    return;
                }
            }
        } else {
            self.ignore_in_report = true;
            return;
        }

        let qual_len = sequence.quality.len();

        if self.current_length < qual_len {
            // Extend existing tile arrays
            for qc_vec in self.per_tile_quality_counts.values_mut() {
                qc_vec.resize_with(qual_len, QualityCount::new);
            }
            self.current_length = qual_len;
        }

        if !self.per_tile_quality_counts.contains_key(&tile) {
            if self.per_tile_quality_counts.len() > 2500 {
                eprintln!("Too many tiles (>2500) so giving up trying to do per-tile qualities since we're probably parsing the file wrongly");
                self.ignore_in_report = true;
                self.per_tile_quality_counts.clear();
                return;
            }

            let mut quality_counts = Vec::with_capacity(self.current_length);
            for _ in 0..self.current_length {
                quality_counts.push(QualityCount::new());
            }
            self.per_tile_quality_counts.insert(tile, quality_counts);
        }

        let quality_counts = self.per_tile_quality_counts.get_mut(&tile).unwrap();
        for (i, ch) in sequence.quality.chars().enumerate() {
            quality_counts[i].add_value(ch);
        }
    }

    fn reset(&mut self) {
        self.total_count = 0;
        self.per_tile_quality_counts.clear();
    }

    fn raises_error(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }
        self.max_deviation > self.config.get_param("tile", "error")
    }

    fn raises_warning(&mut self) -> bool {
        if !self.calculated {
            self.get_percentages();
        }
        self.max_deviation > self.config.get_param("tile", "warn")
    }

    fn make_data_report(&mut self, buf: &mut String) {
        if !self.calculated {
            self.get_percentages();
        }

        buf.push_str("#Tile\tBase\tMean\n");

        for (t, &tile) in self.tiles.iter().enumerate() {
            for i in 0..self.means[t].len() {
                buf.push_str(&tile.to_string());
                buf.push('\t');
                buf.push_str(&self.x_labels[i]);
                buf.push('\t');
                buf.push_str(&format_java_double(self.means[t][i]));
                buf.push('\n');
            }
        }
    }

    fn chart_data(&mut self) -> Option<ChartData> {
        if !self.calculated {
            self.get_percentages();
        }
        if self.tiles.is_empty() {
            return None;
        }

        Some(ChartData::TileHeatmap {
            tile_ids: self.tiles.clone(),
            x_labels: self.x_labels.clone(),
            deviations: self.means.clone(),
            max_deviation: self.max_deviation,
        })
    }
}
