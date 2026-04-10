// Hide console window on Windows when running as GUI
#![cfg_attr(
    all(target_os = "windows", feature = "gui"),
    windows_subsystem = "windows"
)]

#[cfg(feature = "gui")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.is_empty() {
        eprintln!("Usage: fastqc-rs-gui <input.fastq>");
        std::process::exit(1);
    }

    let path = std::path::PathBuf::from(&args[0]);
    if !path.exists() {
        eprintln!("Error: file not found: {}", path.display());
        std::process::exit(1);
    }

    fastqc_rs::gui::run_gui(&path)
}

#[cfg(not(feature = "gui"))]
fn main() {
    eprintln!("This binary requires the 'gui' feature. Build with: cargo build --features gui");
    std::process::exit(1);
}
