// Hide console window on Windows when running as GUI
#![cfg_attr(
    all(target_os = "windows", feature = "gui"),
    windows_subsystem = "windows"
)]

/// Entry point for the `fastqc-rs-gui` binary; opens the iced GUI, optionally with a file argument.
#[cfg(feature = "gui")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    let file_path = if args.is_empty() {
        None
    } else {
        let path = std::path::PathBuf::from(&args[0]);
        if !path.exists() {
            eprintln!("Error: file not found: {}", path.display());
            std::process::exit(1);
        }
        Some(path)
    };

    fastqc_rs::gui::run_gui(file_path.as_deref())
}

#[cfg(not(feature = "gui"))]
fn main() {
    eprintln!("This binary requires the 'gui' feature. Build with: cargo build --features gui");
    std::process::exit(1);
}
