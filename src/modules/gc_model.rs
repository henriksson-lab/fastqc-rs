use crate::encoding::phred::java_round;

/// A single percentage bucket with its weight for GC distribution.
#[derive(Debug, Clone)]
pub struct GCModelValue {
    pub percentage: usize,
    pub increment: f64,
}

/// Model for mapping GC base count to percentage bins for a given read length.
/// Exactly mirrors Java GCModel behavior including java_round for rounding.
#[derive(Debug, Clone)]
pub struct GCModel {
    pub read_length: usize,
    models: Vec<Vec<GCModelValue>>,
}

impl GCModel {
    /// Precompute, for each possible GC base count in a read of length `read_length`,
    /// the set of percentage bins it contributes to along with fractional weights.
    pub fn new(read_length: usize) -> Self {
        let mut claiming_counts = vec![0i32; 101];

        // First pass: count how many GC counts claim each percentage bucket
        for pos in 0..=read_length {
            let mut low_count = pos as f64 - 0.5;
            let mut high_count = pos as f64 + 0.5;

            if low_count < 0.0 {
                low_count = 0.0;
            }
            if high_count < 0.0 {
                high_count = 0.0;
            }
            if high_count > read_length as f64 {
                high_count = read_length as f64;
            }
            if low_count > read_length as f64 {
                low_count = read_length as f64;
            }

            let low_percentage = java_round((low_count * 100.0) / read_length as f64) as usize;
            let high_percentage = java_round((high_count * 100.0) / read_length as f64) as usize;

            for item in claiming_counts
                .iter_mut()
                .take(high_percentage + 1)
                .skip(low_percentage)
            {
                *item += 1;
            }
        }

        // Second pass: create model values with weights
        let mut models = Vec::with_capacity(read_length + 1);
        for pos in 0..=read_length {
            let mut low_count = pos as f64 - 0.5;
            let mut high_count = pos as f64 + 0.5;

            if low_count < 0.0 {
                low_count = 0.0;
            }
            if high_count < 0.0 {
                high_count = 0.0;
            }
            if high_count > read_length as f64 {
                high_count = read_length as f64;
            }
            if low_count > read_length as f64 {
                low_count = read_length as f64;
            }

            let low_percentage = java_round((low_count * 100.0) / read_length as f64) as usize;
            let high_percentage = java_round((high_count * 100.0) / read_length as f64) as usize;

            let mut model_values = Vec::with_capacity(high_percentage - low_percentage + 1);
            for (p, &claim_count) in claiming_counts
                .iter()
                .enumerate()
                .take(high_percentage + 1)
                .skip(low_percentage)
            {
                model_values.push(GCModelValue {
                    percentage: p,
                    increment: 1.0 / claim_count as f64,
                });
            }
            models.push(model_values);
        }

        Self {
            read_length,
            models,
        }
    }

    /// Return the percentage-bin contributions for a read with `gc_count` G/C bases.
    pub fn get_model_values(&self, gc_count: usize) -> &[GCModelValue] {
        &self.models[gc_count]
    }
}
