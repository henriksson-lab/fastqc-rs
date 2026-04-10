use plotters::prelude::*;

use super::colors::HotColdGradient;
use super::layout;
use super::ChartData;

/// Render a tile heatmap matching Java FastQC's TileGraph.paint().
pub fn render<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    data: &ChartData,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let (tile_ids, x_labels, deviations, max_deviation) = match data {
        ChartData::TileHeatmap {
            tile_ids,
            x_labels,
            deviations,
            max_deviation,
        } => (tile_ids, x_labels, deviations, *max_deviation),
        _ => return Err("Expected TileHeatmap data".into()),
    };

    let (w, h) = root.dim_in_pixel();
    let w = w as i32;
    let h = h as i32;

    root.fill(&WHITE)?;

    let n_cols = x_labels.len();
    let n_rows = tile_ids.len();
    if n_cols == 0 || n_rows == 0 {
        return Ok(());
    }

    let x_offset = layout::MARGIN_LEFT_MIN + 10; // extra space for tile ID labels
    let y_top = layout::MARGIN_TOP;
    let y_bottom = h - layout::MARGIN_BOTTOM;
    let x_right = w - layout::MARGIN_RIGHT;

    let cell_width = (x_right - x_offset) as f64 / n_cols as f64;
    let cell_height = (y_bottom - y_top) as f64 / n_rows as f64;

    let gradient = HotColdGradient::new();

    // Draw title
    root.draw(&Text::new(
        "Quality per tile".to_string(),
        (w / 2 - 50, 15),
        ("sans-serif", layout::TITLE_FONT_SIZE)
            .into_font()
            .color(&BLACK),
    ))?;

    // Draw heatmap cells
    for (row, tile_devs) in deviations.iter().enumerate() {
        for (col, &dev) in tile_devs.iter().enumerate() {
            let x1 = x_offset + (col as f64 * cell_width) as i32;
            let y1 = y_top + (row as f64 * cell_height) as i32;
            let x2 = x_offset + ((col + 1) as f64 * cell_width) as i32;
            let y2 = y_top + ((row + 1) as f64 * cell_height) as i32;

            // Map deviation to color: 0 = green, negative = blue, positive = red
            let color = gradient.get_color(
                (dev + max_deviation) / 2.0,
                0.0,
                max_deviation,
            );
            root.draw(&Rectangle::new([(x1, y1), (x2, y2)], color.filled()))?;
        }
    }

    // Y-axis: tile ID labels (skip if too dense)
    let label_skip = if n_rows > 50 { n_rows / 25 } else { 1 };
    for (row, &tile_id) in tile_ids.iter().enumerate() {
        if row % label_skip != 0 {
            continue;
        }
        let y = y_top + ((row as f64 + 0.5) * cell_height) as i32;
        root.draw(&Text::new(
            tile_id.to_string(),
            (5, y - 5),
            ("sans-serif", 8.0).into_font().color(&BLACK),
        ))?;
    }

    // X-axis labels
    let x_skip = if n_cols > 50 { n_cols / 25 } else { 1 };
    for (i, label) in x_labels.iter().enumerate() {
        if i % x_skip != 0 && i != n_cols - 1 {
            continue;
        }
        let x = x_offset + ((i as f64 + 0.5) * cell_width) as i32;
        root.draw(&Text::new(
            label.clone(),
            (x - 10, y_bottom + 5),
            ("sans-serif", 9.0).into_font().color(&BLACK),
        ))?;
    }

    // X-axis title
    root.draw(&Text::new(
        "Position in read (bp)".to_string(),
        (w / 2 - 50, h - 10),
        ("sans-serif", layout::LABEL_FONT_SIZE)
            .into_font()
            .color(&BLACK),
    ))?;

    // Axes
    root.draw(&PathElement::new(
        vec![(x_offset, y_top), (x_offset, y_bottom)],
        BLACK.stroke_width(1),
    ))?;
    root.draw(&PathElement::new(
        vec![(x_offset, y_bottom), (x_right, y_bottom)],
        BLACK.stroke_width(1),
    ))?;

    Ok(())
}
