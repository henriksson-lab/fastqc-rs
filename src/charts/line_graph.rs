use plotters::prelude::*;

use super::colors::{GRID_COLOR, STRIPE_COLOR, TOL_PALETTE};
use super::layout;
use super::ChartData;

/// Render a line graph matching Java FastQC's LineGraph.paint().
pub fn render<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    data: &ChartData,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let (series, x_labels, x_axis_label, title, min_y, max_y) = match data {
        ChartData::LineGraph {
            series,
            x_labels,
            x_axis_label,
            title,
            min_y,
            max_y,
        } => (series, x_labels, x_axis_label, title, *min_y, *max_y),
        _ => return Err("Expected LineGraph data".into()),
    };

    let (w, h) = root.dim_in_pixel();
    let w = w as i32;
    let h = h as i32;

    // White background
    root.fill(&WHITE)?;

    let n = x_labels.len();
    if n == 0 {
        return Ok(());
    }

    // Layout calculations
    let x_offset = layout::MARGIN_LEFT_MIN;
    let y_top = layout::MARGIN_TOP;
    let y_bottom = h - layout::MARGIN_BOTTOM;
    let x_right = w - layout::MARGIN_RIGHT;
    let plot_width = x_right - x_offset;
    let plot_height = y_bottom - y_top;

    // Draw title
    root.draw(&Text::new(
        title.clone(),
        (w / 2, 15),
        ("sans-serif", layout::TITLE_FONT_SIZE)
            .into_font()
            .color(&BLACK),
    ))?;

    // Draw alternating vertical stripes
    let bar_width = plot_width as f64 / n as f64;
    for i in (1..n).step_by(2) {
        let x1 = x_offset + (i as f64 * bar_width) as i32;
        let x2 = x_offset + ((i + 1) as f64 * bar_width) as i32;
        root.draw(&Rectangle::new(
            [(x1, y_top), (x2.min(x_right), y_bottom)],
            STRIPE_COLOR.filled(),
        ))?;
    }

    // Y-axis grid lines
    let y_interval = layout::find_optimal_y_interval(max_y);
    let mut y_val = y_interval;
    while y_val < max_y {
        let y_px = y_bottom - ((y_val - min_y) / (max_y - min_y) * plot_height as f64) as i32;
        if y_px > y_top && y_px < y_bottom {
            root.draw(&PathElement::new(
                vec![(x_offset, y_px), (x_right, y_px)],
                GRID_COLOR.stroke_width(1),
            ))?;
            // Y-axis label
            root.draw(&Text::new(
                format_y_label(y_val),
                (x_offset - 5, y_px - 5),
                ("sans-serif", layout::LABEL_FONT_SIZE)
                    .into_font()
                    .color(&BLACK),
            ))?;
        }
        y_val += y_interval;
    }

    // Draw axes
    root.draw(&PathElement::new(
        vec![(x_offset, y_top), (x_offset, y_bottom)],
        BLACK.stroke_width(1),
    ))?;
    root.draw(&PathElement::new(
        vec![(x_offset, y_bottom), (x_right, y_bottom)],
        BLACK.stroke_width(1),
    ))?;

    // X-axis labels (with collision avoidance)
    let max_label_chars = 6;
    let skip = if n > 50 { n / 25 } else { 1 };
    for (i, label) in x_labels.iter().enumerate() {
        if i % skip != 0 && i != n - 1 {
            continue;
        }
        let x_px = x_offset + ((i as f64 + 0.5) * bar_width) as i32;
        let display = if label.len() > max_label_chars {
            &label[..max_label_chars]
        } else {
            label
        };
        root.draw(&Text::new(
            display.to_string(),
            (x_px - 10, y_bottom + 5),
            ("sans-serif", 9.0).into_font().color(&BLACK),
        ))?;
    }

    // X-axis title
    root.draw(&Text::new(
        x_axis_label.clone(),
        (w / 2 - 50, h - 10),
        ("sans-serif", layout::LABEL_FONT_SIZE)
            .into_font()
            .color(&BLACK),
    ))?;

    // Draw data series
    for (s_idx, s) in series.iter().enumerate() {
        let color = TOL_PALETTE[s_idx % TOL_PALETTE.len()];
        let points: Vec<(i32, i32)> = s
            .data
            .iter()
            .enumerate()
            .filter(|(_, &v)| v.is_finite())
            .map(|(i, &v)| {
                let x_px = x_offset + ((i as f64 + 0.5) * bar_width) as i32;
                let y_frac = (v - min_y) / (max_y - min_y);
                let y_px = y_bottom - (y_frac * plot_height as f64) as i32;
                (x_px, y_px)
            })
            .collect();

        if points.len() >= 2 {
            root.draw(&PathElement::new(points, color.stroke_width(2)))?;
        }
    }

    // Draw legend
    if !series.is_empty() {
        let legend_x = x_right - 180;
        let legend_y = y_top + 5;
        let line_height = 16;
        let legend_h = series.len() as i32 * line_height + 10;

        // Legend background
        root.draw(&Rectangle::new(
            [(legend_x, legend_y), (legend_x + 170, legend_y + legend_h)],
            WHITE.filled(),
        ))?;
        root.draw(&Rectangle::new(
            [(legend_x, legend_y), (legend_x + 170, legend_y + legend_h)],
            RGBColor(200, 200, 200).stroke_width(1),
        ))?;

        for (i, s) in series.iter().enumerate() {
            let color = TOL_PALETTE[i % TOL_PALETTE.len()];
            let y = legend_y + 8 + i as i32 * line_height;
            // Color line
            root.draw(&PathElement::new(
                vec![(legend_x + 5, y), (legend_x + 25, y)],
                color.stroke_width(2),
            ))?;
            // Label
            root.draw(&Text::new(
                s.name.clone(),
                (legend_x + 30, y - 5),
                ("sans-serif", 10.0).into_font().color(&BLACK),
            ))?;
        }
    }

    Ok(())
}

/// Format a y-axis tick: integer if whole, otherwise one decimal place.
fn format_y_label(val: f64) -> String {
    if val == val.floor() {
        format!("{}", val as i64)
    } else {
        format!("{:.1}", val)
    }
}
