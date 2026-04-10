use crate::modules::{QCModule, QCStatus};

// Embed icon PNGs as base64 at compile time
const FASTQC_ICON_PNG: &[u8] = include_bytes!("../resources/icons/fastqc_icon.png");
const TICK_PNG: &[u8] = include_bytes!("../resources/icons/tick.png");
const WARNING_PNG: &[u8] = include_bytes!("../resources/icons/warning.png");
const ERROR_PNG: &[u8] = include_bytes!("../resources/icons/error.png");

// Embed CSS template at compile time
const HEADER_CSS: &str = include_str!("../resources/header_template.html");

/// Generate the fastqc_data.txt content from all modules.
/// Mirrors Java HTMLReportArchive text data generation.
pub fn generate_data_report(
    modules: &mut [Box<dyn QCModule>],
    version: &str,
) -> String {
    let mut buf = String::new();

    // Header
    buf.push_str("##FastQC\t");
    buf.push_str(version);
    buf.push('\n');

    for module in modules.iter_mut() {
        if module.ignore_in_report() {
            continue;
        }

        // Module header with status
        let status = module.status();
        buf.push_str(">>");
        buf.push_str(module.name());
        buf.push('\t');
        buf.push_str(&status.to_string());
        buf.push('\n');

        // Module data
        module.make_data_report(&mut buf);

        // Module footer
        buf.push_str(">>END_MODULE\n");
    }

    buf
}

/// Encode bytes as base64 string.
fn to_base64(data: &[u8]) -> String {
    use base64::Engine;
    base64::engine::general_purpose::STANDARD.encode(data)
}

/// Return the base64 data URI for a status icon.
fn status_icon_base64(status: QCStatus) -> String {
    let png_data = match status {
        QCStatus::Pass => TICK_PNG,
        QCStatus::Warn => WARNING_PNG,
        QCStatus::Fail => ERROR_PNG,
    };
    format!("data:image/png;base64,{}", to_base64(png_data))
}

/// Return the alt text for a status icon.
fn status_alt(status: QCStatus) -> &'static str {
    match status {
        QCStatus::Pass => "[PASS]",
        QCStatus::Warn => "[WARNING]",
        QCStatus::Fail => "[FAIL]",
    }
}

/// Convert a module's text data report into an HTML table.
/// Lines starting with '#' are treated as header rows (column names),
/// other lines are tab-separated data rows.
fn module_data_to_html_table(data_text: &str) -> String {
    let mut buf = String::new();
    let mut in_thead = false;
    let mut in_tbody = false;

    buf.push_str("<table>\n");

    for line in data_text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('#') {
            // Header row: strip leading '#' and split by tab
            let header = line.trim_start_matches('#');
            let cols: Vec<&str> = header.split('\t').collect();

            if !in_thead {
                buf.push_str("     <thead>\n");
                in_thead = true;
            }
            buf.push_str("      <tr>\n");
            for col in &cols {
                buf.push_str("       <th>");
                buf.push_str(col.trim());
                buf.push_str("</th>\n");
            }
            buf.push_str("      </tr>\n");
        } else {
            // Close thead if open, open tbody
            if in_thead {
                buf.push_str("     </thead>\n");
                in_thead = false;
            }
            if !in_tbody {
                buf.push_str("     <tbody>\n");
                in_tbody = true;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            buf.push_str("      <tr>\n");
            for col in &cols {
                buf.push_str("       <td>");
                buf.push_str(col);
                buf.push_str("</td>\n");
            }
            buf.push_str("      </tr>\n");
        }
    }

    if in_thead {
        buf.push_str("     </thead>\n");
    }
    if in_tbody {
        buf.push_str("     </tbody>\n");
    }
    buf.push_str("    </table>\n");

    buf
}

/// Generate the fastqc_report.html content from all modules.
/// Mirrors the Java FastQC HTML report format.
pub fn generate_html_report(
    modules: &mut [Box<dyn QCModule>],
    file_name: &str,
    version: &str,
) -> String {
    let mut buf = String::new();
    let fastqc_icon_b64 = to_base64(FASTQC_ICON_PNG);

    // Collect module info: (index_in_report, name, status)
    let mut module_info: Vec<(usize, String, QCStatus)> = Vec::new();
    let mut idx = 0usize;
    for module in modules.iter_mut() {
        if module.ignore_in_report() {
            continue;
        }
        module_info.push((idx, module.name().to_string(), module.status()));
        idx += 1;
    }

    // HTML head
    buf.push_str("<!DOCTYPE html>\n<html>\n <head>\n");
    buf.push_str("  <title>");
    buf.push_str(file_name);
    buf.push_str(" FastQC Report</title>\n");
    buf.push_str("  <style type=\"text/css\">\n");
    buf.push_str(HEADER_CSS);
    buf.push_str("</style>\n");
    buf.push_str(" </head>\n");

    // Body
    buf.push_str(" <body>\n");

    // Header div
    buf.push_str("  <div class=\"header\">\n");
    buf.push_str("   <div id=\"header_title\">\n");
    buf.push_str("    <img src=\"data:image/png;base64,");
    buf.push_str(&fastqc_icon_b64);
    buf.push_str("\" alt=\"FastQC\">FastQC Report\n");
    buf.push_str("   </div>\n");
    buf.push_str("   <div id=\"header_filename\">\n");
    buf.push_str("    ");
    buf.push_str(file_name);
    buf.push('\n');
    buf.push_str("   </div>\n");
    buf.push_str("  </div>\n");

    // Summary sidebar
    buf.push_str("  <div class=\"summary\">\n");
    buf.push_str("   <h2>Summary</h2>\n");
    buf.push_str("   <ul>\n");
    for (mod_idx, name, status) in &module_info {
        let icon_uri = status_icon_base64(*status);
        let alt = status_alt(*status);
        buf.push_str("    <li><img src=\"");
        buf.push_str(&icon_uri);
        buf.push_str("\" alt=\"");
        buf.push_str(alt);
        buf.push_str("\"><a href=\"#M");
        buf.push_str(&mod_idx.to_string());
        buf.push_str("\">");
        buf.push_str(name);
        buf.push_str("</a></li>\n");
    }
    buf.push_str("   </ul>\n");
    buf.push_str("  </div>\n");

    // Main div with module sections
    buf.push_str("  <div class=\"main\">\n");

    let mut report_idx = 0usize;
    for module in modules.iter_mut() {
        if module.ignore_in_report() {
            continue;
        }

        let status = module.status();
        let icon_uri = status_icon_base64(status);
        let alt = status_alt(status);
        let name = module.name().to_string();

        buf.push_str("   <div class=\"module\">\n");
        buf.push_str("    <h2 id=\"M");
        buf.push_str(&report_idx.to_string());
        buf.push_str("\"><img src=\"");
        buf.push_str(&icon_uri);
        buf.push_str("\" alt=\"");
        buf.push_str(alt);
        buf.push_str("\">");
        buf.push_str(&name);
        buf.push_str("</h2>\n");

        // Render chart image if available
        if let Some(chart_data) = module.chart_data() {
            match crate::charts::render_chart_to_png(&chart_data) {
                Ok(png_bytes) => {
                    let b64 = to_base64(&png_bytes);
                    buf.push_str("    <p><img class=\"indented\" src=\"data:image/png;base64,");
                    buf.push_str(&b64);
                    buf.push_str("\" alt=\"");
                    buf.push_str(&name);
                    buf.push_str("\"></p>\n");
                }
                Err(e) => {
                    eprintln!("Warning: failed to render chart for {}: {}", name, e);
                }
            }
        }

        // Generate text data and convert to HTML table
        let mut data_text = String::new();
        module.make_data_report(&mut data_text);
        let table_html = module_data_to_html_table(&data_text);
        buf.push_str("    ");
        buf.push_str(&table_html);

        buf.push_str("   </div>\n");
        report_idx += 1;
    }

    buf.push_str("  </div>\n");

    // Footer
    buf.push_str("  <div class=\"footer\">\n");
    buf.push_str("   Produced by <a href=\"http://www.bioinformatics.babraham.ac.uk/projects/fastqc/\">FastQC</a> (version ");
    buf.push_str(version);
    buf.push_str(")\n");
    buf.push_str("  </div>\n");

    buf.push_str(" </body>\n</html>\n");

    buf
}
