use plotters::prelude::*;

use super::colors::*;
use super::layout;
use super::ChartData;

/// Render a quality box plot matching Java FastQC's QualityBoxPlot.paint().
pub fn render<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    data: &ChartData,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let (means, medians, lower_q, upper_q, tenth, ninetieth, x_labels, y_min, y_max) = match data {
        ChartData::QualityBoxPlot {
            means,
            medians,
            lower_quartile,
            upper_quartile,
            tenth_percentile,
            ninetieth_percentile,
            x_labels,
            y_min,
            y_max,
        } => (
            means,
            medians,
            lower_quartile,
            upper_quartile,
            tenth_percentile,
            ninetieth_percentile,
            x_labels,
            *y_min,
            *y_max,
        ),
        _ => return Err("Expected QualityBoxPlot data".into()),
    };

    let (w, h) = root.dim_in_pixel();
    let w = w as i32;
    let h = h as i32;

    root.fill(&WHITE)?;

    let n = x_labels.len();
    if n == 0 {
        return Ok(());
    }

    let x_offset = layout::MARGIN_LEFT_MIN;
    let y_top = layout::MARGIN_TOP;
    let y_bottom = h - layout::MARGIN_BOTTOM;
    let x_right = w - layout::MARGIN_RIGHT;
    let plot_height = (y_bottom - y_top) as f64;
    let bar_width = (x_right - x_offset) as f64 / n as f64;

    // Helper: map quality value to y pixel
    let y_px = |val: f64| -> i32 {
        let frac = (val - y_min) / (y_max - y_min);
        y_bottom - (frac * plot_height) as i32
    };

    // Draw title
    root.draw(&Text::new(
        "Quality scores across all bases".to_string(),
        (w / 2 - 100, 15),
        ("sans-serif", layout::TITLE_FONT_SIZE)
            .into_font()
            .color(&BLACK),
    ))?;

    // Draw quality zone backgrounds
    let q28_y = y_px(28.0);
    let q20_y = y_px(20.0);
    // Green zone (Q28+)
    if q28_y > y_top {
        root.draw(&Rectangle::new(
            [(x_offset, y_top), (x_right, q28_y.min(y_bottom))],
            QUALITY_GOOD.filled(),
        ))?;
    }
    // Amber zone (Q20-28)
    if q20_y > y_top {
        let zone_top = q28_y.max(y_top);
        root.draw(&Rectangle::new(
            [(x_offset, zone_top), (x_right, q20_y.min(y_bottom))],
            QUALITY_BAD.filled(),
        ))?;
    }
    // Red zone (Q<20)
    {
        let zone_top = q20_y.max(y_top);
        root.draw(&Rectangle::new(
            [(x_offset, zone_top), (x_right, y_bottom)],
            QUALITY_UGLY.filled(),
        ))?;
    }

    // Draw alternating darker stripes in quality zones
    for i in (1..n).step_by(2) {
        let x1 = x_offset + (i as f64 * bar_width) as i32;
        let x2 = x_offset + ((i + 1) as f64 * bar_width) as i32;
        // Green stripe
        if q28_y > y_top {
            root.draw(&Rectangle::new(
                [(x1, y_top), (x2.min(x_right), q28_y.min(y_bottom))],
                QUALITY_GOOD_DARK.filled(),
            ))?;
        }
        // Amber stripe
        if q20_y > y_top {
            let zt = q28_y.max(y_top);
            root.draw(&Rectangle::new(
                [(x1, zt), (x2.min(x_right), q20_y.min(y_bottom))],
                QUALITY_BAD_DARK.filled(),
            ))?;
        }
        // Red stripe
        {
            let zt = q20_y.max(y_top);
            root.draw(&Rectangle::new(
                [(x1, zt), (x2.min(x_right), y_bottom)],
                QUALITY_UGLY_DARK.filled(),
            ))?;
        }
    }

    // Draw box plots for each position
    let box_inset = (bar_width * 0.15) as i32;
    for i in 0..n {
        let x_center = x_offset + ((i as f64 + 0.5) * bar_width) as i32;
        let x_left = x_offset + (i as f64 * bar_width) as i32 + box_inset;
        let x_right_box = x_offset + ((i + 1) as f64 * bar_width) as i32 - box_inset;

        if !medians[i].is_finite() {
            continue;
        }

        let med_y = y_px(medians[i]);
        let lq_y = y_px(lower_q[i]);
        let uq_y = y_px(upper_q[i]);
        let low_y = y_px(tenth[i]);
        let high_y = y_px(ninetieth[i]);

        // Whisker lines (10th to 25th, 75th to 90th)
        root.draw(&PathElement::new(
            vec![(x_center, low_y), (x_center, lq_y)],
            BLACK.stroke_width(1),
        ))?;
        root.draw(&PathElement::new(
            vec![(x_center, uq_y), (x_center, high_y)],
            BLACK.stroke_width(1),
        ))?;

        // Whisker caps
        let cap = (bar_width * 0.25) as i32;
        root.draw(&PathElement::new(
            vec![(x_center - cap, low_y), (x_center + cap, low_y)],
            BLACK.stroke_width(1),
        ))?;
        root.draw(&PathElement::new(
            vec![(x_center - cap, high_y), (x_center + cap, high_y)],
            BLACK.stroke_width(1),
        ))?;

        // IQR box (yellow fill)
        root.draw(&Rectangle::new(
            [(x_left, uq_y), (x_right_box, lq_y)],
            BOX_FILL.filled(),
        ))?;
        // Box outline
        root.draw(&Rectangle::new(
            [(x_left, uq_y), (x_right_box, lq_y)],
            BLACK.stroke_width(1),
        ))?;

        // Median line (red)
        root.draw(&PathElement::new(
            vec![(x_left, med_y), (x_right_box, med_y)],
            MEDIAN_COLOR.stroke_width(2),
        ))?;
    }

    // Mean line (blue, connecting all positions)
    let mean_points: Vec<(i32, i32)> = means
        .iter()
        .enumerate()
        .filter(|(_, &v)| v.is_finite())
        .map(|(i, &v)| {
            let x = x_offset + ((i as f64 + 0.5) * bar_width) as i32;
            (x, y_px(v))
        })
        .collect();
    if mean_points.len() >= 2 {
        root.draw(&PathElement::new(mean_points, MEAN_COLOR.stroke_width(1)))?;
    }

    // Axes
    root.draw(&PathElement::new(
        vec![(x_offset, y_top), (x_offset, y_bottom)],
        BLACK.stroke_width(1),
    ))?;
    root.draw(&PathElement::new(
        vec![(x_offset, y_bottom), (w - layout::MARGIN_RIGHT, y_bottom)],
        BLACK.stroke_width(1),
    ))?;

    // Y-axis labels
    let y_interval = layout::find_optimal_y_interval(y_max);
    let mut yv = 0.0;
    while yv <= y_max {
        let yp = y_px(yv);
        if yp >= y_top && yp <= y_bottom {
            root.draw(&Text::new(
                format!("{}", yv as i32),
                (x_offset - 25, yp - 5),
                ("sans-serif", layout::LABEL_FONT_SIZE)
                    .into_font()
                    .color(&BLACK),
            ))?;
        }
        yv += y_interval;
    }

    // X-axis labels
    let skip = if n > 50 { n / 25 } else { 1 };
    for (i, label) in x_labels.iter().enumerate() {
        if i % skip != 0 && i != n - 1 {
            continue;
        }
        let x = x_offset + ((i as f64 + 0.5) * bar_width) as i32;
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

    Ok(())
}
