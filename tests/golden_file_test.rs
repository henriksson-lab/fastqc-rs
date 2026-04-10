use std::path::PathBuf;

use fastqc_rs::config::FastQCConfig;
use fastqc_rs::run_fastqc_on_file;

fn compare_with_golden(test_name: &str) {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fastq_path = base.join(format!("fastqc/test/data/{}.fastq", test_name));
    let golden_path = base.join(format!(
        "fastqc/test/integration/FileContentsTest_{}_fastqc_data.approved.txt",
        test_name
    ));

    assert!(
        fastq_path.exists(),
        "Test file not found: {}",
        fastq_path.display()
    );
    assert!(
        golden_path.exists(),
        "Golden file not found: {}",
        golden_path.display()
    );

    let config = FastQCConfig {
        quiet: true,
        ..Default::default()
    };

    let result = run_fastqc_on_file(&fastq_path, &config).expect("Analysis failed");

    let golden = std::fs::read_to_string(&golden_path).expect("Failed to read golden file");

    // Replace version in our output to match golden file version
    let report = result.data_report.replacen(
        &format!("##FastQC\t{}", env!("CARGO_PKG_VERSION")),
        "##FastQC\t0.12.1",
        1,
    );

    if report != golden {
        // Find first difference for debugging
        let report_lines: Vec<&str> = report.lines().collect();
        let golden_lines: Vec<&str> = golden.lines().collect();
        for (i, (r, g)) in report_lines.iter().zip(golden_lines.iter()).enumerate() {
            if r != g {
                panic!(
                    "Mismatch at line {} for {}:\n  got:      '{}'\n  expected: '{}'",
                    i + 1,
                    test_name,
                    r,
                    g
                );
            }
        }
        if report_lines.len() != golden_lines.len() {
            panic!(
                "Line count mismatch for {}: got {}, expected {}",
                test_name,
                report_lines.len(),
                golden_lines.len()
            );
        }
    }
}

#[test]
fn test_minimal_fastq_matches_golden() {
    compare_with_golden("minimal");
}

#[test]
fn test_complex_fastq_matches_golden() {
    compare_with_golden("complex");
}
