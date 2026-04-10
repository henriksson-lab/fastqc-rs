use iced::widget::{button, column, container, horizontal_rule, image, row, scrollable, text};
use iced::{Element, Length, Task, Theme};

use crate::charts::{self, ChartData};
use crate::config::FastQCConfig;
use crate::modules::QCStatus;

/// Pre-computed result for a single module, safe to use from GUI thread.
#[derive(Debug, Clone)]
pub struct ModuleResult {
    pub name: String,
    pub status: QCStatus,
    pub data_text: String,
    pub chart_data: Option<ChartData>,
}

/// Main GUI application state.
pub struct FastQCApp {
    modules: Vec<ModuleResult>,
    selected: usize,
    chart_cache: Vec<Option<iced::widget::image::Handle>>,
    file_name: String,
}

#[derive(Debug, Clone)]
pub enum Message {
    SelectModule(usize),
    OpenFile,
}

impl FastQCApp {
    pub fn new(modules: Vec<ModuleResult>, file_name: String) -> (Self, Task<Message>) {
        let n = modules.len();
        let mut chart_cache = vec![None; n];

        // Pre-render first module's chart
        if !modules.is_empty() {
            if let Some(ref cd) = modules[0].chart_data {
                if let Ok(png) = charts::render_chart_to_png(cd) {
                    chart_cache[0] = Some(iced::widget::image::Handle::from_bytes(png));
                }
            }
        }

        (
            Self {
                modules,
                selected: 0,
                chart_cache,
                file_name,
            },
            Task::none(),
        )
    }

    pub fn title(&self) -> String {
        format!("FastQC - {}", self.file_name)
    }

    pub fn update(&mut self, message: Message) -> Task<Message> {
        match message {
            Message::SelectModule(idx) => {
                self.selected = idx;
                // Render chart if not cached
                if self.chart_cache[idx].is_none() {
                    if let Some(ref cd) = self.modules[idx].chart_data {
                        if let Ok(png) = charts::render_chart_to_png(cd) {
                            self.chart_cache[idx] =
                                Some(iced::widget::image::Handle::from_bytes(png));
                        }
                    }
                }
            }
            Message::OpenFile => {
                // TODO: file open dialog
            }
        }
        Task::none()
    }

    pub fn view(&self) -> Element<'_, Message> {
        // Sidebar: module list
        let sidebar_content: Vec<Element<Message>> = self
            .modules
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let status_char = match m.status {
                    QCStatus::Pass => "✓",
                    QCStatus::Warn => "⚠",
                    QCStatus::Fail => "✗",
                };
                let label = format!("{} {}", status_char, m.name);
                let btn = button(text(label).size(13))
                    .on_press(Message::SelectModule(i))
                    .width(Length::Fill);

                if i == self.selected {
                    container(btn)
                        .style(|_theme: &Theme| container::Style {
                            background: Some(iced::Background::Color(iced::Color::from_rgb(
                                0.85, 0.85, 0.85,
                            ))),
                            ..Default::default()
                        })
                        .into()
                } else {
                    container(btn).into()
                }
            })
            .collect();

        let sidebar = container(scrollable(
            column(sidebar_content).spacing(2).width(Length::Fill),
        ))
        .width(250)
        .height(Length::Fill);

        // Main panel: chart + data
        let main_content = if self.modules.is_empty() {
            column![text("No modules loaded").size(20)]
        } else {
            let module = &self.modules[self.selected];
            let mut content = column![
                text(&module.name).size(20),
                horizontal_rule(1),
            ]
            .spacing(10);

            // Chart image
            if let Some(ref handle) = self.chart_cache[self.selected] {
                content = content.push(
                    image(handle.clone())
                        .width(Length::Fill)
                        .height(Length::Shrink),
                );
            }

            // Data table as text
            content = content.push(
                scrollable(text(&module.data_text).size(11).font(iced::Font::MONOSPACE))
                    .height(Length::Fill),
            );

            content
        };

        let main_panel = container(main_content.padding(10))
            .width(Length::Fill)
            .height(Length::Fill);

        row![sidebar, main_panel].into()
    }
}

/// Run QC analysis on a file and launch the GUI.
pub fn run_gui(file_path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };

    let file_name = file_path
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    eprintln!("Analyzing {}...", file_name);

    // Run analysis
    let mut modules_raw = crate::create_modules(&config);
    let format = crate::detect_format(file_path, &config);
    let only_mapped = format == "bam_mapped" || format == "sam_mapped";

    match format.as_str() {
        "bam" | "bam_mapped" => {
            let bam =
                crate::sequence::bam_file::BamFileReader::open(file_path, only_mapped)?;
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
            let sam =
                crate::sequence::bam_file::SamFileReader::open(file_path, only_mapped)?;
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
                file_path,
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

    // Extract module results for GUI
    let mut results: Vec<ModuleResult> = Vec::new();
    for module in modules_raw.iter_mut() {
        if module.ignore_in_report() {
            continue;
        }
        let mut data_text = String::new();
        module.make_data_report(&mut data_text);
        results.push(ModuleResult {
            name: module.name().to_string(),
            status: module.status(),
            data_text,
            chart_data: module.chart_data(),
        });
    }

    eprintln!("Analysis complete. Launching GUI...");

    // Launch iced application
    iced::application(
        FastQCApp::title,
        FastQCApp::update,
        FastQCApp::view,
    )
    .window_size((1024.0, 700.0))
    .run_with(move || FastQCApp::new(results, file_name))?;

    Ok(())
}
