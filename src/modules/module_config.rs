use std::collections::HashMap;
use std::path::Path;

/// Module threshold configuration, mirroring Java ModuleConfig exactly.
/// Reads limits.txt to determine warn/error thresholds.
pub struct ModuleConfig {
    parameters: HashMap<String, f64>,
}

const DEFAULT_LIMITS: &str = include_str!("../resources/limits.txt");

impl ModuleConfig {
    pub fn new(limits_file: Option<&Path>) -> Self {
        let mut params = HashMap::new();

        // Set defaults matching Java ModuleConfig.readParams()
        params.insert("duplication:warn".to_string(), 70.0);
        params.insert("duplication:error".to_string(), 50.0);
        params.insert("kmer:warn".to_string(), 2.0);
        params.insert("kmer:error".to_string(), 5.0);
        params.insert("n_content:warn".to_string(), 5.0);
        params.insert("n_content:error".to_string(), 20.0);
        params.insert("overrepresented:warn".to_string(), 0.1);
        params.insert("overrepresented:error".to_string(), 1.0);
        params.insert("quality_base_lower:warn".to_string(), 10.0);
        params.insert("quality_base_lower:error".to_string(), 5.0);
        params.insert("quality_base_median:warn".to_string(), 25.0);
        params.insert("quality_base_median:error".to_string(), 20.0);
        params.insert("sequence:warn".to_string(), 10.0);
        params.insert("sequence:error".to_string(), 20.0);
        params.insert("gc_sequence:warn".to_string(), 15.0);
        params.insert("gc_sequence:error".to_string(), 30.0);
        params.insert("quality_sequence:warn".to_string(), 20.0);
        params.insert("quality_sequence:error".to_string(), 27.0);
        params.insert("tile:warn".to_string(), 5.0);
        params.insert("tile:error".to_string(), 10.0);
        params.insert("sequence_length:warn".to_string(), 1.0);
        params.insert("sequence_length:error".to_string(), 1.0);
        params.insert("adapter:warn".to_string(), 5.0);
        params.insert("adapter:error".to_string(), 10.0);

        // Ignore flags (default 0 = not ignored)
        params.insert("duplication:ignore".to_string(), 0.0);
        params.insert("kmer:ignore".to_string(), 0.0);
        params.insert("n_content:ignore".to_string(), 0.0);
        params.insert("overrepresented:ignore".to_string(), 0.0);
        params.insert("quality_base:ignore".to_string(), 0.0);
        params.insert("sequence:ignore".to_string(), 0.0);
        params.insert("gc_sequence:ignore".to_string(), 0.0);
        params.insert("quality_sequence:ignore".to_string(), 0.0);
        params.insert("tile:ignore".to_string(), 0.0);
        params.insert("sequence_length:ignore".to_string(), 0.0);
        params.insert("adapter:ignore".to_string(), 0.0);

        // Override with file contents
        let text = if let Some(path) = limits_file {
            std::fs::read_to_string(path).unwrap_or_else(|e| {
                eprintln!("Warning: could not read limits file {:?}: {}", path, e);
                DEFAULT_LIMITS.to_string()
            })
        } else {
            DEFAULT_LIMITS.to_string()
        };

        Self::parse_limits(&text, &mut params);

        Self { parameters: params }
    }

    fn parse_limits(text: &str, params: &mut HashMap<String, f64>) {
        for line in text.lines() {
            let line = line.trim();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let sections: Vec<&str> = line.split_whitespace().collect();
            if sections.len() != 3 {
                eprintln!(
                    "Config line '{}' didn't contain the 3 required sections",
                    line
                );
                continue;
            }

            if sections[1] != "warn" && sections[1] != "error" && sections[1] != "ignore" {
                eprintln!(
                    "Second config field must be error, warn or ignore, not '{}'",
                    sections[1]
                );
                continue;
            }

            let value: f64 = match sections[2].parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("Value {} didn't look like a number", sections[2]);
                    continue;
                }
            };

            let key = format!("{}:{}", sections[0], sections[1]);
            params.insert(key, value);
        }
    }

    /// Get a parameter value. Panics if key not found (matching Java behavior).
    pub fn get_param(&self, module: &str, level: &str) -> f64 {
        let key = format!("{}:{}", module, level);
        *self.parameters.get(&key).unwrap_or_else(|| {
            panic!("No key called {} in the config data", key);
        })
    }
}

impl Default for ModuleConfig {
    fn default() -> Self {
        Self::new(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_defaults() {
        let config = ModuleConfig::default();
        assert_eq!(config.get_param("duplication", "warn"), 70.0);
        assert_eq!(config.get_param("quality_base_lower", "error"), 5.0);
        assert_eq!(config.get_param("overrepresented", "warn"), 0.1);
    }
}
