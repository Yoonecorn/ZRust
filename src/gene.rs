use crate::{DrugzError, error::finite, fold_change::Comparison, input::Row};
use std::collections::BTreeMap;
#[derive(Debug)]
pub(crate) struct Gene {
    pub name: String,
    pub sum: f64,
    pub count: usize,
    pub pre: f64,
    pub norm: f64,
    pub p: [f64; 2],
    pub fdr: [f64; 2],
    pub rank: [usize; 2],
}
pub(crate) fn aggregate(
    rows: &[Row],
    comparisons: &[Comparison],
    minobs: usize,
) -> Result<Vec<Gene>, DrugzError> {
    let mut genes: BTreeMap<&str, (f64, usize)> = BTreeMap::new();
    for &i in &comparisons[0].order {
        let entry = genes.entry(&rows[i].gene).or_default();
        for comparison in comparisons {
            let z = comparison.z[i];
            entry.0 += z;
            if z != 0.0 {
                entry.1 = entry.1.checked_add(1).ok_or_else(|| {
                    DrugzError::Numeric(format!(
                        "gene {:?}: observation count overflow",
                        rows[i].gene
                    ))
                })?;
            }
        }
    }
    let mut result = Vec::new();
    for (name, (sum, count)) in genes {
        if count < minobs {
            continue;
        }
        finite(sum, format!("gene {name:?}: nonfinite sumZ"))?;
        let pre = finite(
            sum / (count as f64).sqrt(),
            format!("gene {name:?}: nonfinite pre-standardization score"),
        )?;
        result.push(Gene {
            name: name.into(),
            sum,
            count,
            pre,
            norm: 0.0,
            p: [0.0; 2],
            fdr: [0.0; 2],
            rank: [0; 2],
        });
    }
    if result.len() < 2 {
        return Err(DrugzError::Numeric(format!(
            "cross-gene standardization: fewer than two retained genes after minobs={minobs}"
        )));
    }
    Ok(result)
}
