# fastqc-compliant-rs

A Rust port of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — a quality control tool for high-throughput sequencing data.

`fastqc-compliant-rs` targets output compatibility with the original Java FastQC, verified against golden test files for core FASTQ reports. It supports FASTQ (plain, gzip, bzip2), BAM/SAM, and legacy Nanopore Fast5/HDF5 input formats. The project scope is now full original FastQC coverage: remaining Java FastQC behaviors such as `fastqc.fo` archive generation are tracked as compatibility work rather than intentionally excluded.

## This is an LLM-mediated faithful (hopefully) translation, not the original code! 

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this manifesto](https://rewrites.bio/) but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back into the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information
* **If you are the author of the original code and wish to move to Rust, you can obtain ownership of this repository and crate**. Until then, our commitment is to offer an as-faithful-as-possible translation of a snapshot of your code. If we find serious bugs, we will report them to you. Otherwise we will just replicate them, to ensure comparability across studies that claim to use package XYZ v.666. Think of this like a fancy Ubuntu .deb-package of your software - that is how we treat it

This blurb might be out of date. Go to [this page](https://github.com/henriksson-lab/rustification) for the latest information and further information about how we approach translation

## Installation

### Pre-built binaries

Download from the [Releases](../../releases) page:
- **Windows**: `fastqc-compliant-rs-windows.zip` — portable, no installer needed. Extract and run `fastqc-compliant-rs.exe` (CLI) or `fastqc-compliant-rs-gui.exe` (GUI).
- **macOS**: `fastqc-compliant-rs-macos.zip` — extract `FastQC.app` and drag to Applications.
- **Linux**: `fastqc-compliant-rs-linux.tar.gz` — extract and run.

(these might not be generated yet; but you can produce them yourself from source)

### From source

```sh
git clone https://github.com/henriksson-lab/fastqc-compliant-rs.git
cd fastqc-compliant-rs
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
# Run: target/release/fastqc-compliant-rs-gui input.fastq
```

### Packaging for distribution

```sh
# Windows portable exe (cross-compile from Linux)
bash package/build-windows.sh

# macOS .app bundle
bash package/build-macos.sh
```

### Release checklist

Before publishing a release, run:

```sh
cargo fmt --check
cargo clippy --all-targets --all-features
cargo test
```

Also run the Java FastQC golden comparisons for representative FASTQ fixtures, and do a manual smoke run on representative FASTQ, BAM, SAM, and Fast5 inputs where fixtures are available. Do not claim complete FastQC replacement status in release notes until the parity matrix below has no unsupported items.

Fast5/HDF5 support uses the crates.io `hdf5-pure-rust` dependency, not native `libhdf5`.

## Performance

Local single-threaded timing on a 123,888,712 byte gzipped FASTQ (`Undetermined_S0_L001_R1_001.fastq.gz`) showed `fastqc-compliant-rs` completing the same HTML+ZIP, no-extract workflow about 1.30x faster than Java FastQC v0.11.9.

Benchmark environment: Linux 6.8, Intel Xeon Gold 6138, release build from this repository, 1 warmup run followed by 3 measured runs per tool.

| Tool | Command shape | Mean time | Range |
|------|---------------|-----------|-------|
| `fastqc-compliant-rs 0.4.0` | `target/release/fastqc-compliant-rs --threads 1 --quiet --noextract --outdir <tmp> Undetermined_S0_L001_R1_001.fastq.gz` | 9.64s | 9.59-9.69s |
| `FastQC v0.11.9` | `fastqc --threads 1 --quiet --noextract --outdir <tmp> Undetermined_S0_L001_R1_001.fastq.gz` | 12.55s | 12.21-13.20s |

These numbers are machine- and input-dependent; use them as a snapshot, not a general guarantee.

## CLI Usage

```sh
# Basic usage
fastqc-compliant-rs input.fastq

# Multiple files
fastqc-compliant-rs sample1.fastq sample2.fastq.gz sample3.bam

# Specify output directory
fastqc-compliant-rs --outdir results/ input.fastq

# Suppress progress output
fastqc-compliant-rs --quiet input.fastq

# Use multiple threads for parallel file processing
fastqc-compliant-rs --threads 4 *.fastq.gz

# Force input format
fastqc-compliant-rs --format bam input.bam

# Custom adapter/contaminant lists
fastqc-compliant-rs --adapters my_adapters.txt --contaminants my_contaminants.txt input.fastq

# All options
fastqc-compliant-rs --help
```

### Output Files

For each input file (e.g., `sample.fastq`), FastQC-style output files are generated:
- `sample_fastqc.html` — top-level HTML report; by default it references `Icons/` and `Images/` assets like FastQC's archived report
- `sample_fastqc.zip` — ZIP archive containing `sample_fastqc/fastqc_data.txt`, `sample_fastqc/fastqc_report.html`, `sample_fastqc/summary.txt`, `Icons/`, and `Images/`

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
| `--noextract` | Do not unzip output archive |
| `--delete` | Delete zip after extraction |
| `--memory` | Accepted for FastQC CLI compatibility; validated but otherwise ignored |
| `--dir`, `-d` | Accepted for FastQC CLI compatibility; validates temp directory |
| `--java`, `-j` | Accepted for FastQC CLI compatibility; ignored |
| `--min_length` | Minimum base-grouping length for report modules, matching FastQC; does not filter reads |
| `--dup_length` | Truncation length for duplication detection (default: 50) |
| `--svg` | Generate SVG images instead of PNG |
| `--embed-images` | Embed icons and charts as base64 data URIs in HTML instead of referencing `Icons/` and `Images/` |

## Library Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
fastqc-compliant-rs = "0.4"
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
println!("{}", report.html_report);   // HTML report; references Icons/ and Images/ by default
println!("{}", report.summary_report); // summary.txt content

// Access structured per-module results
for module in &report.modules {
    println!("{}: {}", module.name, module.status);
    // module.data_text  — tab-separated data for this module
    // module.chart_data — ChartData enum for rendering charts
}
```

### Write FastQC-compatible report artifacts

The library returns the text reports and rendered chart images directly. Use `report::write_fastqc_archive` to create the same FastQC-compatible zip layout as the CLI, or persist selected artifacts yourself.

```rust
use fastqc_rs::config::FastQCConfig;
use fastqc_rs::report;
use fastqc_rs::FastQCRunner;
use std::path::Path;

let runner = FastQCRunner::new(FastQCConfig::default());
let report = runner.run_file(Path::new("input.fastq")).unwrap();

report::write_fastqc_archive(
    "sample_fastqc.zip",
    "sample",
    &report.data_report,
    &report.html_report,
    &report.summary_report,
    &report.chart_images,
).unwrap();

std::fs::write("fastqc_data.txt", &report.data_report).unwrap();
std::fs::write("fastqc_report.html", &report.html_report).unwrap();
std::fs::write("summary.txt", &report.summary_report).unwrap();

for image in &report.chart_images {
    std::fs::write(&image.filename, &image.bytes).unwrap();
    println!("wrote {} ({})", image.filename, image.mime_type);
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

### Render SVG charts

```rust
use fastqc_rs::charts;

let report = runner.run_file(Path::new("input.fastq")).unwrap();

for module in &report.modules {
    if let Some(ref chart_data) = module.chart_data {
        let svg_bytes = charts::render_chart_to_svg(chart_data).unwrap();
        std::fs::write(format!("{}.svg", module.name), &svg_bytes).unwrap();
    }
}
```

### Key types

| Type | Description |
|------|-------------|
| `FastQCRunner` | Main entry point — holds config, runs analysis |
| `FastQCReport` | Analysis result: `data_report`, `html_report`, `summary_report`, `chart_images`, `modules` |
| `ModuleResult` | Per-module result: `name`, `status`, `data_text`, `chart_data` |
| `ChartImage` | Rendered archive image: `filename`, `mime_type`, `bytes` |
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
| **fastqc-compliant-rs** (native CPU) | **11.4s** | **57 MB** |

### Small dataset (100k reads x 150bp, 31 MB FASTQ)

| Implementation | Wall time | Speedup |
|----------------|-----------|---------|
| Java FastQC (JDK 19) | 3.23s | 1.0x |
| **fastqc-compliant-rs** (native CPU) | **0.94s** | **3.4x** |

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
| Nanopore Fast5 | `.fast5` | Supported for legacy single-read and multi-read layouts with embedded FastQ datasets |

## FastQC Parity

`fastqc-compliant-rs` targets FastQC CLI/report parity for non-GUI workflows. The GUI is intentionally different, and chart rendering may differ visually, but the text report and archive layout are kept close to Java FastQC.

| Area | Implemented | Tested | Notes |
|------|-------------|--------|-------|
| 12 standard modules | Yes | Golden FASTQ + focused edge tests | Module order matches Java FastQC |
| FASTQ | Yes | Yes | Plain, gzip, multi-member gzip, bzip2, colorspace |
| BAM | Yes | Partial | Needs tiny BAM parity fixtures |
| SAM | Yes | Yes | Autodetect, forced format, mapped-only soft clipping, reverse-complement handling |
| CASAVA grouping | Yes | Yes | Invalid names become singleton groups |
| Nanopore grouping | Yes | Yes | Directory scan/grouping plus legacy Fast5/HDF5 extraction |
| CLI flags | Mostly | Yes | `--memory`, `--dir`, and `--java` are compatibility-only |
| `fastqc_data.txt` | Yes | Golden FASTQ | More Java golden fixtures planned |
| `summary.txt` | Yes | Yes | Stored in zip and extracted folder |
| HTML report | Yes | Structure tests | References archived `Icons/` and `Images/` by default; `--embed-images` keeps standalone data URIs |
| Zip archive | Mostly | Yes | `fastqc.fo` generation is planned |
| SVG output | Yes | Yes | HTML references SVG chart images by default and archive stores `.svg` chart images |

## Migrating From FastQC

Most CLI usage maps directly:

```sh
# Java FastQC
fastqc --outdir results --threads 4 sample.fastq.gz

# fastqc-compliant-rs
fastqc-compliant-rs --outdir results --threads 4 sample.fastq.gz
```

Known differences:

- `fastqc-compliant-rs` writes `sample_fastqc.html` and `sample_fastqc.zip`, matching FastQC's top-level report naming. The text report lives inside the zip as `sample_fastqc/fastqc_data.txt`.
- `--memory`, `--dir`, and `--java` are accepted for CLI compatibility. They are validated where applicable, but Rust does not use Java heap or temp image files.
- HTML reports reference `Icons/` and `Images/` assets by default. Use `--embed-images` to keep the previous standalone base64 data URI behavior.
- `--nano` handles filename scanning/grouping and legacy `.fast5` HDF5 extraction for the same embedded FastQ dataset paths used by Java FastQC.
- `fastqc.fo` generation is planned but not implemented yet.
- Chart pixels may differ from Java FastQC, but chart data and `fastqc_data.txt` are the compatibility target.

## License

GPL-3.0, matching the original FastQC license.
