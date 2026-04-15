use crate::modules::QCModule;
use crate::sequence::Sequence;

/// Run all modules on a sequence iterator.
/// Mirrors Java AnalysisRunner.run() logic.
pub fn run_analysis(
    sequences: impl Iterator<Item = Result<Sequence, Box<dyn std::error::Error>>>,
    modules: &mut [Box<dyn QCModule>],
    quiet: bool,
    _min_length: usize,
) -> Result<u64, Box<dyn std::error::Error>> {
    let mut count: u64 = 0;

    for result in sequences {
        let seq = result?;

        for module in modules.iter_mut() {
            if seq.is_filtered && module.ignore_filtered_sequences() {
                continue;
            }
            module.process_sequence(&seq);
        }

        count += 1;

        if !quiet && count.is_multiple_of(100_000) {
            eprintln!("Processed {} sequences", count);
        }
    }

    Ok(count)
}
