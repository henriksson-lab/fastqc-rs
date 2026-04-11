# fastqc-rs

A Rust port of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — a quality control tool for high-throughput sequencing data.

`fastqc-rs` produces output identical to the original Java FastQC, verified against golden test files. It supports FASTQ (plain, gzip, bzip2) and BAM/SAM input formats.

## Installation

### Pre-built binaries

Download from the [Releases](../../releases) page:
- **Windows**: `fastqc-rs-windows.zip` — portable, no installer needed. Extract and run `fastqc-rs.exe` (CLI) or `fastqc-rs-gui.exe` (GUI).
- **macOS**: `fastqc-rs-macos.zip` — extract `FastQC.app` and drag to Applications.
- **Linux**: `fastqc-rs-linux.tar.gz` — extract and run.

### From source

```sh
git clone <repo-url>
cd fastqc-rs
cargo build --release                               # CLI only
cargo build --release --features gui                 # CLI + GUI
```

For best performance, compile with native CPU optimizations:

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### GUI

The GUI is an optional feature. Build with:

```sh
cargo build --release --features gui
# Run: target/release/fastqc-rs-gui input.fastq
```

### Packaging for distribution

```sh
# Windows portable exe (cross-compile from Linux)
bash package/build-windows.sh

# macOS .app bundle
bash package/build-macos.sh
```

## CLI Usage

```sh
# Basic usage
fastqc-rs input.fastq

# Multiple files
fastqc-rs sample1.fastq sample2.fastq.gz sample3.bam

# Specify output directory
fastqc-rs --outdir results/ input.fastq

# Suppress progress output
fastqc-rs --quiet input.fastq

# Use multiple threads for parallel file processing
fastqc-rs --threads 4 *.fastq.gz

# Force input format
fastqc-rs --format bam input.bam

# Custom adapter/contaminant lists
fastqc-rs --adapters my_adapters.txt --contaminants my_contaminants.txt input.fastq

# All options
fastqc-rs --help
```

### Output Files

For each input file (e.g., `sample.fastq`), three output files are generated:
- `sample_fastqc_data.txt` — tab-separated text report
- `sample_fastqc_report.html` — self-contained HTML report with embedded CSS and icons
- `sample_fastqc.zip` — ZIP archive containing both reports

### CLI Options

| Option | Description |
|--------|-------------|
| `--outdir`, `-o` | Output directory |
| `--threads`, `-t` | Number of parallel threads (default: 1) |
| `--contaminants`, `-c` | Custom contaminant list file |
| `--adapters`, `-a` | Custom adapter list file |
| `--limits`, `-l` | Custom warn/error thresholds file |
| `--kmers`, `-k` | K-mer size (2-10, default: 7) |
| `--quiet`, `-q` | Suppress progress output |
| `--nogroup` | Disable base position grouping for reads >50bp |
| `--expgroup` | Use exponential base grouping |
| `--format`, `-f` | Force format: fastq, bam, sam, bam_mapped, sam_mapped |
| `--casava` | CASAVA mode: group files and filter |
| `--nano` | Nanopore mode |
| `--nofilter` | Don't filter CASAVA reads |
| `--extract` | Unzip output archive |
| `--min_length` | Minimum sequence length |
| `--dup_length` | Truncation length for duplication detection (default: 50) |
| `--svg` | Generate SVG images instead of PNG |

## Library Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
fastqc-rs = "0.1"
```

### Run QC on a file

```rust
use fastqc_rs::config::FastQCConfig;
use fastqc_rs::{FastQCRunner, FastQCReport, ModuleResult};
use fastqc_rs::modules::QCStatus;
use std::path::Path;

let config = FastQCConfig::default();
let runner = FastQCRunner::new(config);
let report: FastQCReport = runner.run_file(Path::new("input.fastq")).unwrap();

// Access raw reports
println!("{}", report.data_report);   // fastqc_data.txt content
println!("{}", report.html_report);   // HTML with embedded charts

// Access structured per-module results
for module in &report.modules {
    println!("{}: {}", module.name, module.status);
    // module.data_text  — tab-separated data for this module
    // module.chart_data — ChartData enum for rendering charts
}
```

### Run QC on in-memory sequences

No file I/O required — pass sequences directly:

```rust
use fastqc_rs::config::FastQCConfig;
use fastqc_rs::sequence::Sequence;
use fastqc_rs::FastQCRunner;

let sequences = vec![
    Sequence::new(
        "sample.fastq".into(),
        "ACGTACGTACGT".into(),
        "IIIIIIIIIIII".into(),
        "@read1".into(),
    ),
    Sequence::new(
        "sample.fastq".into(),
        "GCTAGCTAGCTA".into(),
        "HHHHHHHHHHHH".into(),
        "@read2".into(),
    ),
];

let config = FastQCConfig { quiet: true, ..Default::default() };
let runner = FastQCRunner::new(config);
let report = runner.run_sequences(sequences.into_iter()).unwrap();

// Check pass/warn/fail per module
for m in &report.modules {
    match m.status {
        fastqc_rs::modules::QCStatus::Pass => println!("  PASS  {}", m.name),
        fastqc_rs::modules::QCStatus::Warn => println!("  WARN  {}", m.name),
        fastqc_rs::modules::QCStatus::Fail => println!("  FAIL  {}", m.name),
    }
}
```

### Render charts programmatically

```rust
use fastqc_rs::charts;

let report = runner.run_file(Path::new("input.fastq")).unwrap();

for module in &report.modules {
    if let Some(ref chart_data) = module.chart_data {
        let png_bytes = charts::render_chart_to_png(chart_data).unwrap();
        std::fs::write(format!("{}.png", module.name), &png_bytes).unwrap();
    }
}
```

### Key types

| Type | Description |
|------|-------------|
| `FastQCRunner` | Main entry point — holds config, runs analysis |
| `FastQCReport` | Analysis result: `data_report`, `html_report`, `modules` |
| `ModuleResult` | Per-module result: `name`, `status`, `data_text`, `chart_data` |
| `FastQCConfig` | Configuration (mirrors all CLI options) |
| `Sequence` | A single read: `id`, `sequence`, `quality`, `file_name` |
| `QCStatus` | `Pass`, `Warn`, or `Fail` |
| `ChartData` | `LineGraph`, `QualityBoxPlot`, or `TileHeatmap` — for rendering |

## QC Modules

All 12 standard FastQC modules are implemented:

| Module | Description |
|--------|-------------|
| Basic Statistics | Read counts, lengths, GC content |
| Per Base Sequence Quality | Quality score box plots per position |
| Per Tile Sequence Quality | Quality variation across flow cell tiles |
| Per Sequence Quality Scores | Distribution of mean quality per read |
| Per Base Sequence Content | A/T/G/C frequency per position |
| Per Sequence GC Content | GC% distribution with theoretical normal |
| Per Base N Content | Unknown base (N) frequency per position |
| Sequence Length Distribution | Read length histogram |
| Sequence Duplication Levels | Duplicate sequence analysis |
| Overrepresented Sequences | High-frequency sequences with contaminant lookup |
| Adapter Content | Adapter contamination per position |
| Kmer Content | K-mer enrichment analysis |

## Benchmarks

### Real-world dataset (1.47M reads x 126bp, 119MB gzipped FASTQ)

| Implementation | Wall time | Memory |
|----------------|-----------|--------|
| Java FastQC (JDK 19) | 13.2s | 512 MB |
| **fastqc-rs** (native CPU) | **11.4s** | **57 MB** |

### Small dataset (100k reads x 150bp, 31 MB FASTQ)

| Implementation | Wall time | Speedup |
|----------------|-----------|---------|
| Java FastQC (JDK 19) | 3.23s | 1.0x |
| **fastqc-rs** (native CPU) | **0.94s** | **3.4x** |

### Multi-file (4 x 100k reads, parallel)

| Threads | Wall time | Speedup |
|---------|-----------|---------|
| 1 thread | 3.59s | 1.0x |
| **4 threads** | **1.01s** | **3.6x** |

- Hardware: Linux x86_64
- Rust: compiled with `RUSTFLAGS="-C target-cpu=native"`, LTO enabled
- Java: OpenJDK 19, `-Xmx512m`, headless mode
- Output verified identical between both implementations

## Correctness

Output is validated against Java FastQC's approved golden test files:
- `minimal.fastq` (1 read, 16bp) — character-for-character match
- `complex.fastq` (5 reads, 16bp) — character-for-character match
- `benchmark_100k.fastq` (100k reads, 150bp) — character-for-character match

## Supported Input Formats

| Format | Extension | Status |
|--------|-----------|--------|
| FASTQ | `.fastq`, `.fq` | Supported |
| Gzipped FASTQ | `.fastq.gz`, `.fq.gz` | Supported |
| Bzip2 FASTQ | `.fastq.bz2`, `.fq.bz2` | Supported |
| BAM | `.bam`, `.ubam` | Supported |
| SAM | `.sam` | Supported |
| Colorspace | SOLiD format | Supported |

## License

GPL-3.0, matching the original FastQC license.
