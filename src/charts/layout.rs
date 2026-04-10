/// Calculate chart dimensions from data point count.
/// Matches Java: max(800, data_points * 15) x 600.
pub fn chart_dimensions(data_points: usize) -> (u32, u32) {
    let width = 800u32.max(data_points as u32 * 15);
    (width, 600)
}

/// Find the optimal y-axis interval for grid lines.
/// Port of Java's findOptimalYInterval logic.
pub fn find_optimal_y_interval(max_value: f64) -> f64 {
    let bases = [1.0, 2.0, 2.5, 5.0];
    let mut multiplier = 1.0;

    if max_value <= 0.0 {
        return 1.0;
    }

    // Scale down if needed
    while multiplier * 10.0 < max_value {
        multiplier *= 10.0;
    }

    for &base in &bases {
        let interval = base * multiplier;
        let divisions = max_value / interval;
        if divisions <= 10.0 {
            return interval;
        }
    }

    // Fallback
    max_value / 5.0
}

/// Chart layout margins in pixels.
pub const MARGIN_TOP: i32 = 40;
pub const MARGIN_BOTTOM: i32 = 40;
pub const MARGIN_RIGHT: i32 = 10;
pub const MARGIN_LEFT_MIN: i32 = 50;

/// Standard font sizes.
pub const TITLE_FONT_SIZE: f64 = 14.0;
pub const LABEL_FONT_SIZE: f64 = 11.0;
pub const LEGEND_FONT_SIZE: f64 = 11.0;
