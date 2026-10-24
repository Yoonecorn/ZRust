use crate::{DrugzError, error::finite, gene::Gene};
pub(crate) fn sf(z: f64) -> f64 {
    0.5 * libm::erfc(z / std::f64::consts::SQRT_2)
}
pub(crate) fn rank(genes: &mut [Gene]) -> Result<(), DrugzError> {
    let n = genes.len();
    if genes.iter().all(|g| g.pre == genes[0].pre) {
        return Err(DrugzError::Numeric(
            "cross-gene standardization: identical scores (zero variance)".into(),
        ));
    }
    let mean = genes.iter().map(|g| g.pre).sum::<f64>() / n as f64;
    let sd = (genes.iter().map(|g| (g.pre - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
    if !sd.is_finite() || sd <= 0.0 {
        return Err(DrugzError::Numeric(
            "cross-gene standardization: zero or nonfinite population SD".into(),
        ));
    }
    for g in genes.iter_mut() {
        g.norm = finite(
            (g.pre - mean) / sd,
            format!("gene {:?}: nonfinite standardized Z", g.name),
        )?;
        g.p = [sf(-g.norm), sf(g.norm)];
    }
    for direction in 0..2 {
        let mut order: Vec<_> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let primary = if direction == 0 {
                genes[a]
                    .norm
                    .partial_cmp(&genes[b].norm)
                    .expect("finite normalized Z")
            } else {
                genes[b]
                    .norm
                    .partial_cmp(&genes[a].norm)
                    .expect("finite normalized Z")
            };
            primary.then(genes[a].name.as_bytes().cmp(genes[b].name.as_bytes()))
        });
        let mut minimum = f64::INFINITY;
        for (position, &i) in order.iter().enumerate().rev() {
            let rank = position + 1;
            let p = finite(
                genes[i].p[direction],
                format!("gene {:?}: nonfinite p-value", genes[i].name),
            )?;
            minimum = minimum.min(p / (rank as f64 / n as f64));
            finite(minimum, format!("gene {:?}: nonfinite FDR", genes[i].name))?;
            if !(0.0..=1.0).contains(&minimum) {
                return Err(DrugzError::Numeric(format!(
                    "gene {:?}: FDR outside probability domain",
                    genes[i].name
                )));
            }
            genes[i].rank[direction] = rank;
            genes[i].fdr[direction] = minimum;
        }
    }
    genes.sort_by(|a, b| {
        a.norm
            .partial_cmp(&b.norm)
            .expect("finite normalized Z")
            .then(a.name.as_bytes().cmp(b.name.as_bytes()))
    });
    Ok(())
}
