use crate::{AnalysisConfig, DrugzError, error::finite, input::Row};
#[derive(Debug)]
pub(crate) struct Comparison {
    pub order: Vec<usize>,
    pub raw_control: Vec<f64>,
    pub raw_treatment: Vec<f64>,
    pub norm_control: Vec<f64>,
    pub norm_treatment: Vec<f64>,
    pub fc: Vec<f64>,
    pub sd: Vec<f64>,
    pub z: Vec<f64>,
}
pub(crate) fn calculate(
    rows: &[Row],
    norm: &[Vec<f64>],
    c: &AnalysisConfig,
    rep: usize,
) -> Result<Comparison, DrugzError> {
    debug_assert!(rows.windows(2).all(|p| p[0].position < p[1].position));
    let mut result = Comparison {
        order: (0..rows.len()).collect(),
        raw_control: Vec::new(),
        raw_treatment: Vec::new(),
        norm_control: Vec::new(),
        norm_treatment: Vec::new(),
        fc: Vec::new(),
        sd: Vec::new(),
        z: Vec::new(),
    };
    let nc = c.control.len();
    for (row, normalized) in rows.iter().zip(norm) {
        let (rc, rt, cn, tn) = if c.unpaired {
            (
                mean(&row.counts[..nc]),
                mean(&row.counts[nc..]),
                mean(&normalized[..nc]),
                mean(&normalized[nc..]),
            )
        } else {
            (
                row.counts[rep],
                row.counts[nc + rep],
                normalized[rep],
                normalized[nc + rep],
            )
        };
        let control = cn + c.pseudocount as f64;
        let treatment = tn + c.pseudocount as f64;
        if control <= 0.0 || treatment <= 0.0 {
            return Err(DrugzError::Numeric(format!(
                "fold change comparison {rep}, guide {:?}, gene {:?}: logarithm inputs must be positive",
                row.id, row.gene
            )));
        }
        for value in [rc, rt, cn, tn, control, treatment] {
            if !value.is_finite() {
                finite(
                    value,
                    format!(
                        "fold change comparison {rep}, guide {:?}: nonfinite mean or pseudocount addition",
                        row.id
                    ),
                )?;
            }
        }
        let fc = finite(
            (treatment / control).log2(),
            format!(
                "fold change comparison {rep}, guide {:?}: nonfinite log2 ratio",
                row.id
            ),
        )?;
        result.raw_control.push(rc);
        result.raw_treatment.push(rt);
        result.norm_control.push(cn);
        result.norm_treatment.push(tn);
        result.fc.push(fc);
    }
    result.order = crate::ordering::descending(rows.len(), |a, b| {
        // Preserve integer-key distinctions above 2^53 without changing f64 arithmetic.
        if !c.unpaired
            && let (Some(x), Some(y)) = (rows[a].ordering_counts[rep], rows[b].ordering_counts[rep])
        {
            return x < y;
        }
        result.raw_control[a] < result.raw_control[b]
    });
    Ok(result)
}
pub(crate) fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}
