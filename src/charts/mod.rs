pub mod colors;
pub mod layout;
pub mod line_graph;
pub mod quality_box_plot;
pub mod tile_graph;

/// A named data series for line graphs.
#[derive(Debug, Clone)]
pub struct Series {
    pub name: String,
    pub data: Vec<f64>,
}

/// Chart data extracted from a QC module for rendering.
#[derive(Debug, Clone)]
pub enum ChartData {
    /// Line graph with one or more data series.
    LineGraph {
        series: Vec<Series>,
        x_labels: Vec<String>,
        x_axis_label: String,
        title: String,
        min_y: f64,
        max_y: f64,
    },
    /// Box-and-whisker plot for per-base quality scores.
    QualityBoxPlot {
        means: Vec<f64>,
        medians: Vec<f64>,
        lower_quartile: Vec<f64>,
        upper_quartile: Vec<f64>,
        tenth_percentile: Vec<f64>,
        ninetieth_percentile: Vec<f64>,
        x_labels: Vec<String>,
        y_min: f64,
        y_max: f64,
    },
    /// Heatmap for per-tile quality deviations.
    TileHeatmap {
        tile_ids: Vec<i32>,
        x_labels: Vec<String>,
        /// deviations[tile_index][position_group_index]
        deviations: Vec<Vec<f64>>,
        max_deviation: f64,
    },
}

impl ChartData {
    /// Return the expected image filename for this chart type (matches Java output).
    pub fn image_filename(&self, module_name: &str) -> String {
        match module_name {
            "Per base sequence quality" => "per_base_quality.png".into(),
            "Per tile sequence quality" => "per_tile_quality.png".into(),
            "Per sequence quality scores" => "per_sequence_quality.png".into(),
            "Per base sequence content" => "per_base_sequence_content.png".into(),
            "Per sequence GC content" => "per_sequence_gc_content.png".into(),
            "Per base N content" => "per_base_n_content.png".into(),
            "Sequence Length Distribution" => "sequence_length_distribution.png".into(),
            "Sequence Duplication Levels" => "duplication_levels.png".into(),
            "Adapter Content" => "adapter_content.png".into(),
            "Kmer Content" => "kmer_profiles.png".into(),
            _ => format!("{}.png", module_name.to_lowercase().replace(' ', "_")),
        }
    }

    /// Calculate chart dimensions (width, height) matching Java behavior.
    pub fn dimensions(&self) -> (u32, u32) {
        let data_points = match self {
            ChartData::LineGraph { x_labels, .. } => x_labels.len(),
            ChartData::QualityBoxPlot { x_labels, .. } => x_labels.len(),
            ChartData::TileHeatmap { x_labels, .. } => x_labels.len(),
        };
        layout::chart_dimensions(data_points)
    }
}

/// Render a chart to PNG bytes.
pub fn render_chart_to_png(data: &ChartData) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let (w, h) = data.dimensions();
    // RGB buffer
    let mut buf = vec![0u8; (w * h * 3) as usize];
    {
        let backend = BitMapBackend::with_buffer(&mut buf, (w, h));
        let root = backend.into_drawing_area();
        match data {
            ChartData::LineGraph { .. } => line_graph::render(&root, data)?,
            ChartData::QualityBoxPlot { .. } => quality_box_plot::render(&root, data)?,
            ChartData::TileHeatmap { .. } => tile_graph::render(&root, data)?,
        }
        root.present()?;
    }
    // Encode RGB buffer to PNG
    encode_rgb_to_png(&buf, w, h)
}

/// Render a chart to SVG bytes.
pub fn render_chart_to_svg(data: &ChartData) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let (w, h) = data.dimensions();
    let mut svg = String::new();
    {
        let backend = SVGBackend::with_string(&mut svg, (w, h));
        let root = backend.into_drawing_area();
        match data {
            ChartData::LineGraph { .. } => line_graph::render(&root, data)?,
            ChartData::QualityBoxPlot { .. } => quality_box_plot::render(&root, data)?,
            ChartData::TileHeatmap { .. } => tile_graph::render(&root, data)?,
        }
        root.present()?;
    }

    Ok(svg.into_bytes())
}

fn encode_rgb_to_png(rgb: &[u8], w: u32, h: u32) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    use image::{ImageBuffer, Rgb};
    let img: ImageBuffer<Rgb<u8>, _> =
        ImageBuffer::from_raw(w, h, rgb.to_vec()).ok_or("Failed to create image buffer")?;
    let mut png_bytes = Vec::new();
    let encoder = image::codecs::png::PngEncoder::new(&mut png_bytes);
    image::ImageEncoder::write_image(encoder, img.as_raw(), w, h, image::ExtendedColorType::Rgb8)?;
    Ok(png_bytes)
}
