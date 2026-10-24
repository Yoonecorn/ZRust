use crate::{DrugzError, error::finite, fold_change::Comparison, input::Row};
pub(crate) fn calculate(comparison: &mut Comparison, rows: &[Row]) -> Result<(), DrugzError> {
    comparison.z = comparison
        .fc
        .iter()
        .zip(&comparison.sd)
        .zip(rows)
        .map(|((fc, sd), row)| {
            finite(
                fc / sd,
                format!(
                    "guide Z, guide {:?}, gene {:?}: nonfinite score",
                    row.id, row.gene
                ),
            )
        })
        .collect::<Result<_, _>>()?;
    Ok(())
}
