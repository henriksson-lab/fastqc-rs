fn main() {
    // Embed Windows icon and version info into the GUI executable
    #[cfg(target_os = "windows")]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("package/fastqc-rs.ico");
        res.set("ProductName", "fastqc-compliant-rs");
        res.set("FileDescription", "FastQC Quality Control Tool");
        res.set("LegalCopyright", "GPL-3.0");
        if let Err(e) = res.compile() {
            eprintln!("Warning: failed to embed Windows resources: {}", e);
        }
    }
}
