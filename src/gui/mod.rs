use std::path::PathBuf;

use iced::widget::{
    button, center, column, container, horizontal_rule, image, row, scrollable, text,
};
use iced::{Element, Length, Subscription, Task, Theme};

use crate::charts::{self, ChartData};
use crate::config::FastQCConfig;
use crate::modules::QCStatus;

/// Parsed table data from a module's TSV report.
#[derive(Debug, Clone)]
pub struct TableData {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl TableData {
    pub fn from_tsv(tsv: &str) -> Self {
        let mut headers = Vec::new();
        let mut rows = Vec::new();
        for line in tsv.lines() {
            if line.is_empty() {
                continue;
            }
            if let Some(header_line) = line.strip_prefix('#') {
                headers = header_line.split('\t').map(|s| s.to_string()).collect();
            } else {
                rows.push(line.split('\t').map(|s| s.to_string()).collect());
            }
        }
        TableData { headers, rows }
    }
}

/// Pre-computed result for a single module, safe to use from GUI thread.
#[derive(Debug, Clone)]
pub struct ModuleResult {
    pub name: String,
    pub status: QCStatus,
    pub table: TableData,
    pub chart_data: Option<ChartData>,
}

/// Status icon handles, pre-rendered at startup.
struct StatusIcons {
    pass: iced::widget::image::Handle,
    warn: iced::widget::image::Handle,
    fail: iced::widget::image::Handle,
}

/// Render a filled circle with a symbol as a PNG and return an iced image handle.
fn render_status_icon(size: u32, bg: [u8; 3], symbol: &str) -> iced::widget::image::Handle {
    use ::image::{ImageBuffer, Rgba, RgbaImage};

    let mut img: RgbaImage = ImageBuffer::new(size, size);
    let center = size as f64 / 2.0;
    let radius = center - 1.0;

    // Draw filled circle
    for y in 0..size {
        for x in 0..size {
            let dx = x as f64 - center;
            let dy = y as f64 - center;
            let dist = (dx * dx + dy * dy).sqrt();
            if dist <= radius {
                img.put_pixel(x, y, Rgba([bg[0], bg[1], bg[2], 255]));
            } else if dist <= radius + 1.0 {
                // Anti-alias edge
                let alpha = ((radius + 1.0 - dist) * 255.0) as u8;
                img.put_pixel(x, y, Rgba([bg[0], bg[1], bg[2], alpha]));
            }
        }
    }

    let white = Rgba([255, 255, 255, 255]);
    let cx = center as i32;
    let cy = center as i32;
    let s = (size as f64 * 0.22) as i32; // symbol half-size

    match symbol {
        "check" => {
            // Draw a checkmark: two lines forming a V shape
            // Short leg: from (-s, 0) to (-s/3, s*2/3)
            // Long leg: from (-s/3, s*2/3) to (s, -s*2/3)
            draw_thick_line(&mut img, cx - s, cy, cx - s / 3, cy + s * 2 / 3, white, 2);
            draw_thick_line(
                &mut img,
                cx - s / 3,
                cy + s * 2 / 3,
                cx + s,
                cy - s * 2 / 3,
                white,
                2,
            );
        }
        "exclaim" => {
            // Draw exclamation mark: vertical bar + dot
            let bar_top = cy - s;
            let bar_bottom = cy + s / 3;
            draw_thick_line(&mut img, cx, bar_top, cx, bar_bottom, white, 2);
            // Dot
            let dot_y = cy + s;
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let px = (cx + dx) as u32;
                    let py = (dot_y + dy) as u32;
                    if px < size && py < size {
                        img.put_pixel(px, py, white);
                    }
                }
            }
        }
        "cross" => {
            // Draw an X
            draw_thick_line(&mut img, cx - s, cy - s, cx + s, cy + s, white, 2);
            draw_thick_line(&mut img, cx + s, cy - s, cx - s, cy + s, white, 2);
        }
        _ => {}
    }

    // Encode to PNG bytes
    let mut buf = std::io::Cursor::new(Vec::new());
    ::image::write_buffer_with_format(
        &mut buf,
        img.as_raw(),
        size,
        size,
        ::image::ColorType::Rgba8,
        ::image::ImageFormat::Png,
    )
    .expect("PNG encode");

    iced::widget::image::Handle::from_bytes(buf.into_inner())
}

fn draw_thick_line(
    img: &mut ::image::RgbaImage,
    x0: i32,
    y0: i32,
    x1: i32,
    y1: i32,
    color: ::image::Rgba<u8>,
    thickness: i32,
) {
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    let mut x = x0;
    let mut y = y0;
    let (w, h) = img.dimensions();

    loop {
        for ty in -(thickness / 2)..=(thickness / 2) {
            for tx in -(thickness / 2)..=(thickness / 2) {
                let px = (x + tx) as u32;
                let py = (y + ty) as u32;
                if px < w && py < h {
                    img.put_pixel(px, py, color);
                }
            }
        }
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            err += dx;
            y += sy;
        }
    }
}

impl StatusIcons {
    fn new() -> Self {
        let size = 32; // Render at 2x for retina, display at 16px
        Self {
            pass: render_status_icon(size, [0x22, 0xbb, 0x22], "check"), // green
            warn: render_status_icon(size, [0xe6, 0xa8, 0x00], "exclaim"), // orange/yellow
            fail: render_status_icon(size, [0xdd, 0x22, 0x22], "cross"), // red
        }
    }

    fn get(&self, status: &QCStatus) -> iced::widget::image::Handle {
        match status {
            QCStatus::Pass => self.pass.clone(),
            QCStatus::Warn => self.warn.clone(),
            QCStatus::Fail => self.fail.clone(),
        }
    }
}

/// Application view state.
enum AppState {
    Empty,
    Analyzing(String),
    Loaded {
        file_name: String,
        modules: Vec<ModuleResult>,
        selected: usize,
        chart_cache: Vec<Option<iced::widget::image::Handle>>,
    },
    Error(String),
}

/// Main GUI application state.
pub struct FastQCApp {
    state: AppState,
    icons: StatusIcons,
}

#[derive(Debug, Clone)]
pub enum Message {
    SelectModule(usize),
    FileDropped(PathBuf),
    AnalysisComplete(Result<(String, Vec<ModuleResult>), String>),
}

fn analyze_file(path: PathBuf) -> Result<(String, Vec<ModuleResult>), String> {
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };

    let mut modules_raw = crate::create_modules(&config);
    let format = crate::detect_format(&path, &config);
    let only_mapped = format == "bam_mapped" || format == "sam_mapped";

    let mut run = || -> Result<(), Box<dyn std::error::Error>> {
        match format.as_str() {
            "bam" | "bam_mapped" => {
                let bam = crate::sequence::bam_file::BamFileReader::open(&path, only_mapped)?;
                let sequences =
                    bam.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
                crate::analysis::run_analysis(
                    sequences,
                    &mut modules_raw,
                    true,
                    config.min_length,
                )?;
            }
            "sam" | "sam_mapped" => {
                let sam = crate::sequence::bam_file::SamFileReader::open(&path, only_mapped)?;
                let sequences =
                    sam.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
                crate::analysis::run_analysis(
                    sequences,
                    &mut modules_raw,
                    true,
                    config.min_length,
                )?;
            }
            _ => {
                let fq = crate::sequence::fastq_file::FastQFile::open(
                    &path,
                    config.casava,
                    config.nofilter,
                )?;
                let sequences =
                    fq.map(|r| r.map_err(|e| Box::new(e) as Box<dyn std::error::Error>));
                crate::analysis::run_analysis(
                    sequences,
                    &mut modules_raw,
                    true,
                    config.min_length,
                )?;
            }
        }
        Ok(())
    };

    run().map_err(|e| e.to_string())?;

    let mut results = Vec::new();
    for module in modules_raw.iter_mut() {
        if module.ignore_in_report() {
            continue;
        }
        let mut data_text = String::new();
        module.make_data_report(&mut data_text);
        let table = TableData::from_tsv(&data_text);
        results.push(ModuleResult {
            name: module.name().to_string(),
            status: module.status(),
            table,
            chart_data: module.chart_data(),
        });
    }

    Ok((file_name, results))
}

impl FastQCApp {
    fn new(initial_file: Option<PathBuf>) -> (Self, Task<Message>) {
        let icons = StatusIcons::new();
        match initial_file {
            Some(path) => {
                let file_name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_else(|| "unknown".to_string());
                let task =
                    Task::perform(async move { analyze_file(path) }, Message::AnalysisComplete);
                (
                    Self {
                        state: AppState::Analyzing(file_name),
                        icons,
                    },
                    task,
                )
            }
            None => (
                Self {
                    state: AppState::Empty,
                    icons,
                },
                Task::none(),
            ),
        }
    }

    fn title(&self) -> String {
        match &self.state {
            AppState::Loaded { file_name, .. } => format!("FastQC - {}", file_name),
            AppState::Analyzing(name) => format!("FastQC - Analyzing {}...", name),
            _ => "FastQC".to_string(),
        }
    }

    fn update(&mut self, message: Message) -> Task<Message> {
        match message {
            Message::SelectModule(idx) => {
                if let AppState::Loaded {
                    modules,
                    selected,
                    chart_cache,
                    ..
                } = &mut self.state
                {
                    *selected = idx;
                    if chart_cache[idx].is_none() {
                        if let Some(ref cd) = modules[idx].chart_data {
                            if let Ok(png) = charts::render_chart_to_png(cd) {
                                chart_cache[idx] =
                                    Some(iced::widget::image::Handle::from_bytes(png));
                            }
                        }
                    }
                }
            }
            Message::FileDropped(path) => {
                let file_name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_else(|| "unknown".to_string());
                self.state = AppState::Analyzing(file_name);
                return Task::perform(async move { analyze_file(path) }, Message::AnalysisComplete);
            }
            Message::AnalysisComplete(result) => match result {
                Ok((file_name, modules)) => {
                    let n = modules.len();
                    let mut chart_cache = vec![None; n];
                    if !modules.is_empty() {
                        if let Some(ref cd) = modules[0].chart_data {
                            if let Ok(png) = charts::render_chart_to_png(cd) {
                                chart_cache[0] = Some(iced::widget::image::Handle::from_bytes(png));
                            }
                        }
                    }
                    self.state = AppState::Loaded {
                        file_name,
                        modules,
                        selected: 0,
                        chart_cache,
                    };
                }
                Err(e) => {
                    self.state = AppState::Error(e);
                }
            },
        }
        Task::none()
    }

    fn subscription(&self) -> Subscription<Message> {
        iced::event::listen_with(|event, _status, _window| {
            if let iced::event::Event::Window(iced::window::Event::FileDropped(path)) = event {
                Some(Message::FileDropped(path))
            } else {
                None
            }
        })
    }

    fn view(&self) -> Element<'_, Message> {
        match &self.state {
            AppState::Empty => center(
                text("Drop FASTQ file here")
                    .size(28)
                    .color(iced::Color::from_rgb(0.5, 0.5, 0.5)),
            )
            .width(Length::Fill)
            .height(Length::Fill)
            .into(),

            AppState::Analyzing(name) => center(
                column![
                    text(format!("Analyzing {}...", name)).size(24),
                    text("Please wait")
                        .size(16)
                        .color(iced::Color::from_rgb(0.5, 0.5, 0.5)),
                ]
                .spacing(10)
                .align_x(iced::Alignment::Center),
            )
            .width(Length::Fill)
            .height(Length::Fill)
            .into(),

            AppState::Error(msg) => center(
                column![
                    text("Error")
                        .size(24)
                        .color(iced::Color::from_rgb(0.8, 0.2, 0.2)),
                    text(msg).size(16),
                    text("Drop another file to try again")
                        .size(14)
                        .color(iced::Color::from_rgb(0.5, 0.5, 0.5)),
                ]
                .spacing(10)
                .align_x(iced::Alignment::Center),
            )
            .width(Length::Fill)
            .height(Length::Fill)
            .into(),

            AppState::Loaded {
                modules,
                selected,
                chart_cache,
                ..
            } => self.view_loaded(modules, *selected, chart_cache),
        }
    }

    fn view_loaded<'a>(
        &'a self,
        modules: &'a [ModuleResult],
        selected: usize,
        chart_cache: &'a [Option<iced::widget::image::Handle>],
    ) -> Element<'a, Message> {
        // Sidebar
        let sidebar_content: Vec<Element<Message>> = modules
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let icon_handle = self.icons.get(&m.status);
                let icon = image(icon_handle).width(16).height(16);

                let label = text(&m.name).size(13);

                let item = button(
                    row![icon, label]
                        .spacing(6)
                        .align_y(iced::Alignment::Center),
                )
                .on_press(Message::SelectModule(i))
                .width(Length::Fill)
                .padding([4, 8]);

                let bg = if i == selected {
                    iced::Color::from_rgb(0.82, 0.87, 0.95) // light blue selection
                } else {
                    iced::Color::WHITE
                };

                container(item)
                    .style(move |_theme: &Theme| container::Style {
                        background: Some(iced::Background::Color(bg)),
                        ..Default::default()
                    })
                    .into()
            })
            .collect();

        let sidebar = container(scrollable(
            column(sidebar_content).spacing(1).width(Length::Fill),
        ))
        .width(220)
        .height(Length::Fill)
        .style(|_theme: &Theme| container::Style {
            background: Some(iced::Background::Color(iced::Color::from_rgb(
                0.95, 0.95, 0.95,
            ))),
            border: iced::Border {
                color: iced::Color::from_rgb(0.8, 0.8, 0.8),
                width: 1.0,
                radius: 0.0.into(),
            },
            ..Default::default()
        });

        // Main panel
        let main_content = if modules.is_empty() {
            column![text("No modules loaded").size(20)]
        } else {
            let module = &modules[selected];
            let mut content = column![text(&module.name).size(20), horizontal_rule(1),].spacing(10);

            if let Some(ref handle) = chart_cache[selected] {
                content = content.push(
                    image(handle.clone())
                        .width(Length::Fill)
                        .height(Length::Shrink),
                );
            }

            // Data table
            let table = &module.table;
            let mut table_col = column![].spacing(0);

            if !table.headers.is_empty() {
                let header_row: Vec<Element<Message>> = table
                    .headers
                    .iter()
                    .map(|h| {
                        container(text(h).size(11).font(iced::Font::MONOSPACE))
                            .width(Length::FillPortion(1))
                            .padding(2)
                            .style(|_theme: &Theme| container::Style {
                                background: Some(iced::Background::Color(iced::Color::from_rgb(
                                    0.9, 0.9, 0.9,
                                ))),
                                ..Default::default()
                            })
                            .into()
                    })
                    .collect();
                table_col = table_col.push(row(header_row));
            }

            for (i, data_row) in table.rows.iter().enumerate() {
                let bg = if i % 2 == 0 {
                    iced::Color::WHITE
                } else {
                    iced::Color::from_rgb(0.96, 0.96, 0.96)
                };
                let cells: Vec<Element<Message>> = data_row
                    .iter()
                    .map(|cell| {
                        container(text(cell).size(11).font(iced::Font::MONOSPACE))
                            .width(Length::FillPortion(1))
                            .padding(2)
                            .style(move |_theme: &Theme| container::Style {
                                background: Some(iced::Background::Color(bg)),
                                ..Default::default()
                            })
                            .into()
                    })
                    .collect();
                table_col = table_col.push(row(cells));
            }

            content = content.push(scrollable(table_col).height(Length::Fill));
            content
        };

        let main_panel = container(main_content.padding(10))
            .width(Length::Fill)
            .height(Length::Fill);

        row![sidebar, main_panel].into()
    }
}

/// Launch the GUI, optionally with an initial file to analyze.
pub fn run_gui(file_path: Option<&std::path::Path>) -> Result<(), Box<dyn std::error::Error>> {
    let initial = file_path.map(|p| p.to_path_buf());

    iced::application(FastQCApp::title, FastQCApp::update, FastQCApp::view)
        .subscription(FastQCApp::subscription)
        .window_size((1024.0, 700.0))
        .run_with(move || FastQCApp::new(initial))?;

    Ok(())
}
