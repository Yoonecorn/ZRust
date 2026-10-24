//! Deterministic finite-input DrugZ reference profile, qualified on macOS ARM64.
mod design;
mod error;
mod fold_change;
mod gene;
mod input;
mod normalization;
mod ordering;
mod output;
mod provenance;
mod statistics;
mod variance;
mod zscore;
pub use error::DrugzError;

/// Ordered sample lists and finite-reference scientific options.
#[derive(Clone, Debug)]
pub struct AnalysisConfig {
    pub control: Vec<String>,
    pub treatment: Vec<String>,
    pub exclude: Vec<String>,
    pub pseudocount: u64,
    pub minobs: usize,
    pub half_window_size: usize,
    pub unpaired: bool,
}
impl AnalysisConfig {
    /// Validate parameters and arm structure before opening input.
    pub fn validate(&self) -> Result<(), DrugzError> {
        design::validate(self)
    }
}
impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            control: Vec::new(),
            treatment: Vec::new(),
            exclude: Vec::new(),
            pseudocount: 5,
            minobs: 1,
            half_window_size: 500,
            unpaired: false,
        }
    }
}
/// Complete, finite results. Serialization and publication occur only after analysis.
#[derive(Debug)]
pub struct AnalysisResult {
    rows: Vec<input::Row>,
    normalized: Vec<Vec<f64>>,
    comparisons: Vec<fold_change::Comparison>,
    genes: Vec<gene::Gene>,
    config: AnalysisConfig,
    threads: usize,
}
impl AnalysisResult {
    pub fn profile(&self) -> &'static str {
        provenance::PROFILE
    }
    pub fn upstream_commit(&self) -> &'static str {
        provenance::UPSTREAM
    }
    pub fn gene_output(&self) -> String {
        output::genes(self)
    }
    pub fn guide_output(&self) -> String {
        output::guides(self)
    }
    /// Publish guides first, genes last; see the documented rollback limits.
    pub fn publish(
        &self,
        gene: &std::path::Path,
        guide: Option<&std::path::Path>,
    ) -> Result<(), DrugzError> {
        output::publish(self, gene, guide)
    }
    /// Number of retained guides, including guides excluded by gene-level minobs.
    pub fn guide_count(&self) -> usize {
        self.normalized.len()
    }
}
/// Validate UTF-8 TSV bytes and calculate the accepted finite reference profile.
pub fn analyze(input: &[u8], config: &AnalysisConfig) -> Result<AnalysisResult, DrugzError> {
    design::validate(config)?;
    let rows = input::parse(input, config)?;
    let names: Vec<_> = config.control.iter().chain(&config.treatment).collect();
    let normalized = normalization::normalize(&rows, &names)?;
    let mut comparisons = Vec::new();
    for replicate in 0..if config.unpaired {
        1
    } else {
        config.control.len()
    } {
        let mut comparison = fold_change::calculate(&rows, &normalized, config, replicate)?;
        variance::local(&mut comparison, config.half_window_size)?;
        zscore::calculate(&mut comparison, &rows)?;
        comparisons.push(comparison);
    }
    let mut genes = gene::aggregate(&rows, &comparisons, config.minobs)?;
    statistics::rank(&mut genes)?;
    Ok(AnalysisResult {
        rows,
        normalized,
        comparisons,
        genes,
        config: config.clone(),
        threads: 1,
    })
}
#[cfg(test)]
mod tests;

/// Calculate the unchanged serial profile, then format outputs with at most `threads` workers.
/// Arithmetic, validation, and publication retain their serial operation order.
pub fn analyze_with_threads(
    input: &[u8],
    config: &AnalysisConfig,
    threads: usize,
) -> Result<AnalysisResult, DrugzError> {
    if threads == 0 {
        return Err(DrugzError::Parameter("threads must be >= 1".into()));
    }
    let mut result = analyze(input, config)?;
    result.threads = threads;
    Ok(result)
}
