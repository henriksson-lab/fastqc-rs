/// Normal distribution for GC content theoretical curve fitting.
pub struct NormalDistribution {
    pub mean: f64,
    pub stdev: f64,
}

impl NormalDistribution {
    pub fn new(mean: f64, stdev: f64) -> Self {
        Self { mean, stdev }
    }

    /// Probability density function value at x.
    pub fn pdf(&self, x: f64) -> f64 {
        let exponent = -((x - self.mean).powi(2)) / (2.0 * self.stdev * self.stdev);
        (1.0 / (self.stdev * (2.0 * std::f64::consts::PI).sqrt())) * exponent.exp()
    }
}
