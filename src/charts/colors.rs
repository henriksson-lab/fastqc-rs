use plotters::style::RGBColor;

/// Tol colorblind-friendly palette used for line graph series.
/// Matches Java LineGraph.java exactly.
pub const TOL_PALETTE: [RGBColor; 8] = [
    RGBColor(136, 34, 85),   // Purple/Mauve
    RGBColor(51, 34, 136),   // Dark blue
    RGBColor(17, 119, 51),   // Green
    RGBColor(221, 204, 119), // Mustard/Gold
    RGBColor(68, 170, 153),  // Teal
    RGBColor(170, 68, 153),  // Magenta
    RGBColor(204, 102, 119), // Rose
    RGBColor(136, 204, 238), // Light blue
];

// Quality zone background colors (alternating light/dark for visual clarity)
/// Background fill for Q28+ ("good") quality zone in box plots.
pub const QUALITY_GOOD: RGBColor = RGBColor(195, 230, 195);
/// Alternating darker stripe variant of `QUALITY_GOOD`.
pub const QUALITY_GOOD_DARK: RGBColor = RGBColor(175, 230, 175);
/// Background fill for Q20-Q28 ("warning") quality zone.
pub const QUALITY_BAD: RGBColor = RGBColor(230, 220, 195);
/// Alternating darker stripe variant of `QUALITY_BAD`.
pub const QUALITY_BAD_DARK: RGBColor = RGBColor(230, 215, 175);
/// Background fill for Q<20 ("bad") quality zone.
pub const QUALITY_UGLY: RGBColor = RGBColor(230, 195, 195);
/// Alternating darker stripe variant of `QUALITY_UGLY`.
pub const QUALITY_UGLY_DARK: RGBColor = RGBColor(230, 175, 175);

// Box plot element colors
/// Fill color for the IQR box in quality box plots (yellow, matches Java).
pub const BOX_FILL: RGBColor = RGBColor(240, 240, 0);
/// Color of the median line drawn inside each IQR box.
pub const MEDIAN_COLOR: RGBColor = RGBColor(200, 0, 0);
/// Color of the mean line connecting per-position means.
pub const MEAN_COLOR: RGBColor = RGBColor(0, 0, 200);

// Chart background colors
/// Alternating column stripe color used as a subtle background guide.
pub const STRIPE_COLOR: RGBColor = RGBColor(230, 230, 230);
/// Color used for horizontal y-axis grid lines.
pub const GRID_COLOR: RGBColor = RGBColor(180, 180, 180);

/// Hot-cold color gradient for tile quality heatmaps.
/// Blue (cold/negative) → Green (neutral) → Red (hot/positive).
/// Uses sqrt scaling to emphasize deviations.
pub struct HotColdGradient {
    colors: Vec<RGBColor>,
}

impl HotColdGradient {
    /// Creates a palette of 100 pre-cached colours from which the closest match
    /// will be selected to return for future queries. Pre-caching saves on the
    /// overhead of generating a lot of new objects on each lookup.
    pub fn new() -> Self {
        let mut colors = Vec::with_capacity(100);
        for i in 0..100 {
            let fraction = i as f64 / 99.0;
            let r;
            let g;
            let b;
            if fraction < 0.5 {
                // Blue → Green
                let f = fraction * 2.0;
                r = 0;
                g = (255.0 * f) as u8;
                b = (255.0 * (1.0 - f)) as u8;
            } else {
                // Green → Red
                let f = (fraction - 0.5) * 2.0;
                r = (255.0 * f) as u8;
                g = (255.0 * (1.0 - f)) as u8;
                b = 0;
            }
            colors.push(RGBColor(r, g, b));
        }
        Self { colors }
    }

    /// Gets a colour from the gradient for `value`, where `min` and `max` define
    /// the range of the gradient. Applies sqrt scaling to emphasize extremes.
    pub fn get_color(&self, value: f64, min: f64, max: f64) -> RGBColor {
        if max <= min {
            return self.colors[50]; // neutral green
        }
        let range = max - min;
        let fraction = ((value - min) / range).clamp(0.0, 1.0);
        // Apply sqrt transformation to emphasize deviations
        let scaled = fraction.sqrt();
        let index = (scaled * 99.0) as usize;
        self.colors[index.min(99)]
    }
}

impl Default for HotColdGradient {
    fn default() -> Self {
        Self::new()
    }
}
