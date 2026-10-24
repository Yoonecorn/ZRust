use crate::{DrugzError, error::finite, input::Row};
pub(crate) fn normalize(rows: &[Row], names: &[&String]) -> Result<Vec<Vec<f64>>, DrugzError> {
    let mut totals = vec![0.0; names.len()];
    for row in rows {
        for (total, count) in totals.iter_mut().zip(&row.counts) {
            *total += count;
        }
    }
    for (name, total) in names.iter().zip(&totals) {
        if !total.is_finite() || *total <= 0.0 {
            return Err(DrugzError::Numeric(format!(
                "normalization sample {name:?}: total must be finite and positive, got {total}"
            )));
        }
    }
    rows.iter()
        .map(|row| {
            row.counts
                .iter()
                .zip(&totals)
                .enumerate()
                .map(|(j, (count, total))| {
                    finite(
                        1e7 * count / total,
                        format!(
                            "normalization guide {:?}, sample {:?}: nonfinite normalized count",
                            row.id, names[j]
                        ),
                    )
                })
                .collect()
        })
        .collect()
}
