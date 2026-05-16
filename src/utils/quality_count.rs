/// Histogram of quality characters at a single position.
/// Mirrors Java QualityCount exactly: 150 slots indexed by ASCII value.
#[derive(Debug, Clone)]
pub struct QualityCount {
    actual_counts: [u64; 150],
    total_counts: u64,
}

impl QualityCount {
    /// Create an empty histogram with all 150 ASCII slots zeroed.
    pub fn new() -> Self {
        Self {
            actual_counts: [0u64; 150],
            total_counts: 0,
        }
    }

    /// Record one observation of quality character `c`. Panics if its ASCII value exceeds 149.
    pub fn add_value(&mut self, c: char) {
        let idx = c as usize;
        if idx >= self.actual_counts.len() {
            panic!(
                "Got character {} as a quality value which has ASCII {} which is higher than we can cope with",
                c, idx
            );
        }
        self.total_counts += 1;
        self.actual_counts[idx] += 1;
    }

    /// Total number of quality observations recorded.
    pub fn get_total_count(&self) -> u64 {
        self.total_counts
    }

    /// Lowest observed quality character (NUL if no observations).
    pub fn get_min_char(&self) -> char {
        for i in 0..self.actual_counts.len() {
            if self.actual_counts[i] > 0 {
                return char::from(i as u8);
            }
        }
        // Java returns (char)1000 which is out of ASCII range
        char::from(0)
    }

    /// Highest observed quality character (NUL if no observations).
    pub fn get_max_char(&self) -> char {
        for i in (0..self.actual_counts.len()).rev() {
            if self.actual_counts[i] > 0 {
                return char::from(i as u8);
            }
        }
        char::from(0)
    }

    /// Calculate mean quality score adjusted by offset.
    /// Matches Java: total / count where total = sum(count[i] * (i - offset))
    pub fn get_mean(&self, offset: i32) -> f64 {
        let mut total: i64 = 0;
        let mut count: u64 = 0;
        let offset_usize = offset as usize;
        for i in offset_usize..self.actual_counts.len() {
            total += self.actual_counts[i] as i64 * (i as i64 - offset as i64);
            count += self.actual_counts[i];
        }
        total as f64 / count as f64
    }

    /// Calculate percentile quality score adjusted by offset.
    /// Matches Java behavior exactly including NaN for single-count cases.
    pub fn get_percentile(&self, offset: i32, percentile: i32) -> f64 {
        let mut total = self.total_counts;
        total *= percentile as u64;
        total /= 100;

        let mut count: u64 = 0;
        let offset_usize = offset as usize;
        for i in offset_usize..self.actual_counts.len() {
            count += self.actual_counts[i];
            if count >= total {
                return (i as i32 - offset) as f64;
            }
        }
        -1.0
    }
}

impl Default for QualityCount {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_value() {
        let mut qc = QualityCount::new();
        qc.add_value('I'); // ASCII 73
        assert_eq!(qc.get_total_count(), 1);
        assert_eq!(qc.get_min_char(), 'I');
        assert_eq!(qc.get_max_char(), 'I');
        assert_eq!(qc.get_mean(64), 9.0); // 73 - 64 = 9
    }

    #[test]
    fn test_percentile_single_value() {
        let mut qc = QualityCount::new();
        qc.add_value('I'); // ASCII 73
                           // total = 1 * 10 / 100 = 0 (integer division)
                           // Loop starts at i=offset=64, count starts at 0
                           // actualCounts[64] = 0, count = 0, 0 >= 0 => returns (64-64) = 0.0
                           // This matches Java behavior exactly.
        let p = qc.get_percentile(64, 10);
        assert_eq!(p, 0.0);

        // With offset 33 (Sanger): same thing, returns (33-33) = 0.0
        let p2 = qc.get_percentile(33, 10);
        assert_eq!(p2, 0.0);

        // For percentile 100: total = 1*100/100 = 1
        // Loop counts until count >= 1: at i=73, count=1, returns (73-64)=9.0
        let p3 = qc.get_percentile(64, 100);
        assert_eq!(p3, 9.0);
    }
}
