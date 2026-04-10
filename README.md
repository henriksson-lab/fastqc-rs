# fastqc-rs

A Rust port of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — a quality control tool for high-throughput sequencing data.

`fastqc-rs` produces output identical to the original Java FastQC, verified against golden test files. It supports FASTQ (plain, gzip, bzip2) and BAM/SAM input formats.

## Installation

```sh
cargo install fastqc-rs
```

Or build from source:

```sh
git clone <repo-url>
cd fastqc-rs
cargo build --release
```

For best performance, compile with native CPU optimizations:

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release
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

`fastqc-rs` can be used as a library for programmatic QC analysis:

```rust
use fastqc_rs::config::FastQCConfig;
use fastqc_rs::{FastQCRunner, FastQCReport};
use std::path::Path;

// Run QC on a file
let config = FastQCConfig::default();
let runner = FastQCRunner::new(config);
let report: FastQCReport = runner.run_file(Path::new("input.fastq")).unwrap();

println!("{}", report.data_report);  // fastqc_data.txt content
println!("{}", report.html_report);  // HTML report content
```

### In-memory sequences

```rust
use fastqc_rs::config::FastQCConfig;
use fastqc_rs::sequence::Sequence;
use fastqc_rs::FastQCRunner;

let sequences = vec![
    Sequence::new(
        "test.fastq".to_string(),
        "ACGTACGTACGT".to_string(),
        "IIIIIIIIIIII".to_string(),
        "@read1".to_string(),
    ),
];

let config = FastQCConfig { quiet: true, ..Default::default() };
let runner = FastQCRunner::new(config);
let report = runner.run_sequences(sequences.into_iter()).unwrap();
println!("{}", report.data_report);
```

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

### Single file (100k reads x 150bp, 31 MB FASTQ)

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
