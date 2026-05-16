use crate::config::FastQCConfig;

/// A group of base positions in a read, used for binning.
/// Early positions get individual groups, later positions are binned together.
#[derive(Debug, Clone)]
pub struct BaseGroup {
    pub lower_count: usize,
    pub upper_count: usize,
}

impl BaseGroup {
    /// Create a `BaseGroup` covering positions `lower_count..=upper_count` (inclusive, 1-based).
    pub fn new(lower_count: usize, upper_count: usize) -> Self {
        Self {
            lower_count,
            upper_count,
        }
    }

    /// True if `value` falls within this group's inclusive range.
    pub fn contains_value(&self, value: usize) -> bool {
        value >= self.lower_count && value <= self.upper_count
    }

    /// Create base groups matching the Java BaseGroup.makeBaseGroups() logic exactly.
    pub fn make_base_groups(mut max_length: usize, config: &FastQCConfig) -> Vec<BaseGroup> {
        // They might have set a fixed max length
        if config.min_length > max_length {
            max_length = config.min_length;
        }

        if config.nogroup {
            Self::make_ungrouped_groups(max_length)
        } else if config.expgroup {
            Self::make_exponential_base_groups(max_length)
        } else {
            Self::make_linear_base_groups(max_length)
        }
    }

    /// One group per base position (no binning), used when `--nogroup` is set.
    pub fn make_ungrouped_groups(max_length: usize) -> Vec<BaseGroup> {
        let mut groups = Vec::new();
        let interval = 1;
        let mut starting_base = 1;

        while starting_base <= max_length {
            let mut end_base = starting_base + interval - 1;
            if end_base > max_length {
                end_base = max_length;
            }
            groups.push(BaseGroup::new(starting_base, end_base));
            starting_base += interval;
        }

        groups
    }

    /// Exponentially-widening bins (interval grows at 10/50/100/500/1000 bp thresholds).
    /// Matches Java's `makeExponentialBaseGroups` exactly; used when `--expgroup` is set.
    pub fn make_exponential_base_groups(max_length: usize) -> Vec<BaseGroup> {
        let mut groups = Vec::new();
        let mut starting_base: usize = 1;
        let mut interval: usize = 1;

        while starting_base <= max_length {
            let mut end_base = starting_base + interval - 1;
            if end_base > max_length {
                end_base = max_length;
            }
            groups.push(BaseGroup::new(starting_base, end_base));
            starting_base += interval;

            // Increase interval at specific thresholds - matches Java exactly
            if starting_base == 10 && max_length > 75 {
                interval = 5;
            }
            if starting_base == 50 && max_length > 200 {
                interval = 10;
            }
            if starting_base == 100 && max_length > 300 {
                interval = 50;
            }
            if starting_base == 500 && max_length > 1000 {
                interval = 100;
            }
            if starting_base == 1000 && max_length > 2000 {
                interval = 500;
            }
        }

        groups
    }

    /// Pick the smallest interval from {2,5,10} * 10^k that keeps the group count under 75.
    fn get_linear_interval(length: usize) -> usize {
        let base_values = [2, 5, 10];
        let mut multiplier: usize = 1;

        loop {
            for &b in &base_values {
                let interval = b * multiplier;
                let mut group_count = 9 + (length - 9) / interval;
                if !(length - 9).is_multiple_of(interval) {
                    group_count += 1;
                }
                if group_count < 75 {
                    return interval;
                }
            }
            multiplier *= 10;
            if multiplier == 10_000_000 {
                panic!(
                    "Couldn't find a sensible interval grouping for length '{}'",
                    length
                );
            }
        }
    }

    /// Default grouping: first 9 positions individual, then fixed-width bins chosen
    /// to keep the total group count under 75. Reads <=75 bp are returned ungrouped.
    pub fn make_linear_base_groups(max_length: usize) -> Vec<BaseGroup> {
        // For lengths below 75bp we just return everything
        if max_length <= 75 {
            return Self::make_ungrouped_groups(max_length);
        }

        let interval = Self::get_linear_interval(max_length);
        let mut groups = Vec::new();
        let mut starting_base: usize = 1;

        while starting_base <= max_length {
            let mut end_base = starting_base + interval - 1;

            if starting_base < 10 {
                end_base = starting_base;
            }

            if starting_base == 10 && interval > 10 {
                end_base = interval - 1;
            }

            if end_base > max_length {
                end_base = max_length;
            }

            groups.push(BaseGroup::new(starting_base, end_base));

            if starting_base < 10 {
                starting_base += 1;
            } else if starting_base == 10 && interval > 10 {
                starting_base = interval;
            } else {
                starting_base += interval;
            }
        }

        groups
    }
}

impl std::fmt::Display for BaseGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.lower_count == self.upper_count {
            write!(f, "{}", self.lower_count)
        } else {
            write!(f, "{}-{}", self.lower_count, self.upper_count)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ungrouped_16bp() {
        let groups = BaseGroup::make_ungrouped_groups(16);
        assert_eq!(groups.len(), 16);
        assert_eq!(groups[0].lower_count, 1);
        assert_eq!(groups[0].upper_count, 1);
        assert_eq!(groups[15].lower_count, 16);
        assert_eq!(groups[15].upper_count, 16);
    }

    #[test]
    fn test_linear_short() {
        // <= 75bp returns ungrouped
        let groups = BaseGroup::make_linear_base_groups(50);
        assert_eq!(groups.len(), 50);
    }

    #[test]
    fn test_linear_long() {
        let groups = BaseGroup::make_linear_base_groups(150);
        assert!(groups.len() < 75);
        // First 9 should be individual
        for i in 0..9 {
            assert_eq!(groups[i].lower_count, i + 1);
            assert_eq!(groups[i].upper_count, i + 1);
        }
    }
}
