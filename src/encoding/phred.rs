/// Phred quality encoding scheme.
#[derive(Debug, Clone)]
pub struct PhredEncoding {
    pub name: String,
    pub offset: i32,
}

const SANGER_ENCODING_OFFSET: i32 = 33;
const ILLUMINA_1_3_ENCODING_OFFSET: i32 = 64;

impl PhredEncoding {
    /// Detect the Phred encoding from the lowest quality character observed.
    /// Matches Java PhredEncoding.getFastQEncodingOffset() exactly.
    pub fn get_fastq_encoding_offset(lowest_char: char) -> Result<PhredEncoding, String> {
        let val = lowest_char as i32;
        if val < 33 {
            Err(format!(
                "No known encodings with chars < 33 (Yours was '{}' with value {})",
                lowest_char, val
            ))
        } else if val < 64 {
            Ok(PhredEncoding {
                name: "Sanger / Illumina 1.9".to_string(),
                offset: SANGER_ENCODING_OFFSET,
            })
        } else if val == ILLUMINA_1_3_ENCODING_OFFSET + 1 {
            Ok(PhredEncoding {
                name: "Illumina 1.3".to_string(),
                offset: ILLUMINA_1_3_ENCODING_OFFSET,
            })
        } else if val <= 126 {
            Ok(PhredEncoding {
                name: "Illumina 1.5".to_string(),
                offset: ILLUMINA_1_3_ENCODING_OFFSET,
            })
        } else {
            Err(format!(
                "No known encodings with chars > 126 (Yours was {} with value {})",
                lowest_char, val
            ))
        }
    }

    /// Convert a Sanger-scale Phred score to an error probability (10^(-Q/10)).
    pub fn convert_sanger_phred_to_probability(phred: i32) -> f64 {
        10.0_f64.powf(phred as f64 / -10.0)
    }

    /// Convert an old Illumina (Solexa) Phred score to an error probability.
    pub fn convert_old_illumina_phred_to_probability(phred: i32) -> f64 {
        10.0_f64.powf((phred as f64 / (phred as f64 + 1.0)) / -10.0)
    }

    /// Uses Java-compatible rounding (round half up).
    pub fn convert_probability_to_sanger_phred(p: f64) -> i32 {
        java_round(-10.0 * p.log10()) as i32
    }

    /// Uses Java-compatible rounding (round half up).
    pub fn convert_probability_to_old_illumina_phred(p: f64) -> i32 {
        // Note: Java original has a bug: Math.log10(p/1-p) which is (p/1)-p = p-p = 0
        // We reproduce the exact Java behavior: -10d*Math.log10(p/1-p)
        java_round(-10.0 * (p / 1.0 - p).log10()) as i32
    }
}

impl std::fmt::Display for PhredEncoding {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

/// Java-compatible Math.round(): rounds half up for positive numbers.
/// Rust f64::round() uses round-half-to-even (banker's rounding).
pub fn java_round(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sanger_encoding() {
        let enc = PhredEncoding::get_fastq_encoding_offset('!').unwrap(); // ASCII 33
        assert_eq!(enc.offset, 33);
        assert!(enc.name.contains("Sanger"));
    }

    #[test]
    fn test_illumina_1_5_encoding() {
        let enc = PhredEncoding::get_fastq_encoding_offset('I').unwrap(); // ASCII 73
        assert_eq!(enc.offset, 64);
        assert!(enc.name.contains("Illumina 1.5"));
    }

    #[test]
    fn test_illumina_1_3_encoding() {
        let enc = PhredEncoding::get_fastq_encoding_offset('A').unwrap(); // ASCII 65 = 64+1
        assert_eq!(enc.offset, 64);
        assert!(enc.name.contains("Illumina 1.3"));
    }

    #[test]
    fn test_java_round() {
        assert_eq!(java_round(2.5), 3);
        assert_eq!(java_round(3.5), 4);
        assert_eq!(java_round(-0.5), 0);
    }
}
