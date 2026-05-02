fn main() {
    std::process::exit(fastqc_rs::cli::run_cli_from(std::env::args()));
}
