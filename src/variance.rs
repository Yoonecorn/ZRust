use crate::{
    DrugzError,
    error::finite,
    fold_change::{Comparison, mean},
};
pub(crate) fn sample_sd(values: &[f64]) -> f64 {
    let center = mean(values);
    (values.iter().map(|x| (x - center).powi(2)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
}
pub(crate) fn local(comparison: &mut Comparison, h: usize) -> Result<(), DrugzError> {
    let n = comparison.order.len();
    let twice = h
        .checked_mul(2)
        .ok_or_else(|| DrugzError::Parameter("half_window_size overflow".into()))?;
    if n < 2
        || h >= n
            .checked_mul(2)
            .ok_or_else(|| DrugzError::Numeric("window size overflow".into()))?
    {
        return Err(DrugzError::Numeric(format!(
            "local window H={h}, guides={n}: final-tail index infeasible"
        )));
    }
    let limit = n
        .checked_add(25)
        .ok_or_else(|| DrugzError::Numeric("window index overflow".into()))?
        .saturating_sub(h);
    let ordered: Vec<_> = comparison.order.iter().map(|&i| comparison.fc[i]).collect();
    let mut scales = vec![0.0; n];
    let initial = sample_sd(&ordered[..twice.min(n)]);
    scales[..h.min(n)].fill(initial);
    for i in (h..limit).step_by(25) {
        if i > n {
            return Err(DrugzError::Numeric(format!(
                "local window H={h}, guides={n}, step={i}: previous-row index infeasible"
            )));
        }
        let end = i
            .checked_add(h)
            .ok_or_else(|| DrugzError::Numeric("window index overflow".into()))?
            .min(n);
        let start = i - h;
        if end.saturating_sub(start) < 2 {
            return Err(DrugzError::Numeric(format!(
                "local window H={h}, guides={n}, step={i}: fewer than two observations"
            )));
        }
        let sd = sample_sd(&ordered[start..end]);
        finite(
            sd,
            format!("local window H={h}, step={i}: nonfinite sample SD"),
        )?;
        let previous = scales[i - 1];
        scales[i..i.saturating_add(25).min(n)].fill(sd.max(previous));
    }
    let tail_index = if h < n { n - h - 1 } else { n - (h - n) - 1 };
    let tail_start = if h <= n {
        n - h
    } else {
        n.saturating_sub(h - n)
    };
    let tail = scales[tail_index];
    scales[tail_start..].fill(tail);
    comparison.sd = vec![0.0; n];
    for (position, &row) in comparison.order.iter().enumerate() {
        let sd = scales[position];
        if !sd.is_finite() || sd <= 0.0 {
            return Err(DrugzError::Numeric(format!(
                "local SD at ordered guide position {position}, H={h}: expected finite positive scale, got {sd}"
            )));
        }
        comparison.sd[row] = sd;
    }
    Ok(())
}
