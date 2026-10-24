use crate::{AnalysisConfig, DrugzError, design};
use std::collections::BTreeSet;
#[derive(Debug)]
pub(crate) struct Row {
    pub id: String,
    pub gene: String,
    pub position: usize,
    pub counts: Vec<f64>,
    pub integer_columns: Vec<bool>,
    pub ordering_counts: Vec<Option<u64>>,
}
pub(crate) fn parse(bytes: &[u8], c: &AnalysisConfig) -> Result<Vec<Row>, DrugzError> {
    let text =
        std::str::from_utf8(bytes).map_err(|e| DrugzError::Input(format!("invalid UTF-8: {e}")))?;
    let mut lines = text.lines();
    let header: Vec<_> = lines
        .next()
        .ok_or_else(|| DrugzError::Input("missing header".into()))?
        .split('\t')
        .collect();
    if header.len() < 3 {
        return Err(DrugzError::Input(
            "header requires guide, gene and sample columns".into(),
        ));
    }
    let mut samples = BTreeSet::new();
    for (i, name) in header.iter().enumerate().skip(2) {
        if name.is_empty() || !samples.insert(name) {
            return Err(DrugzError::Input(format!(
                "empty or duplicate sample at column {}: {name:?}",
                i + 1
            )));
        }
    }
    let mut ids = BTreeSet::new();
    let mut raw = Vec::new();
    for (position, line) in lines.enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        let row = position + 2;
        if fields.len() != header.len() {
            return Err(DrugzError::Input(format!(
                "row {row}: width {}, expected {}",
                fields.len(),
                header.len()
            )));
        }
        if fields[0].is_empty() || !ids.insert(fields[0]) {
            return Err(DrugzError::Input(format!(
                "row {row}, column 1: empty or duplicate guide {:?}",
                fields[0]
            )));
        }
        if fields[1].is_empty() {
            return Err(DrugzError::Input(format!(
                "row {row}, column 2, guide {:?}: empty gene",
                fields[0]
            )));
        }
        raw.push(fields);
    }
    let selected = design::resolve(&header, c)?;
    // Match pandas column dtype for paired raw-count serialization, including excluded rows.
    let integer_columns: Vec<bool> = selected
        .iter()
        .map(|&i| {
            raw.iter()
                .all(|fields| fields[i].trim().parse::<i64>().is_ok())
        })
        .collect();
    let mut rows = Vec::new();
    let integer_key = |s: &str| {
        s.trim().parse::<u64>().ok().or_else(|| {
            s.trim()
                .parse::<i64>()
                .ok()
                .and_then(|x| u64::try_from(x).ok())
        })
    };
    let ordering_integer_columns: Vec<_> = selected
        .iter()
        .map(|&i| raw.iter().all(|f| integer_key(f[i]).is_some()))
        .collect();
    for (position, fields) in raw.iter().enumerate() {
        let counts = selected.iter().map(|&i| {
            let context = || format!("row {}, column {}, sample {:?}, gene {:?}, guide {:?}: invalid finite nonnegative count {:?}", position + 2, i + 1, header[i], fields[1], fields[0], fields[i]);
            let value: f64 = fields[i].trim().parse().map_err(|_| DrugzError::Input(context()))?;
            if !value.is_finite() || value < 0.0 { return Err(DrugzError::Input(context())); }
            Ok(value)
        }).collect::<Result<_, _>>()?;
        if !c.exclude.iter().any(|g| g == fields[1]) {
            rows.push(Row {
                id: fields[0].into(),
                gene: fields[1].into(),
                position,
                counts,
                integer_columns: integer_columns.clone(),
                ordering_counts: selected
                    .iter()
                    .enumerate()
                    .map(|(column, &i)| {
                        if ordering_integer_columns[column] {
                            integer_key(fields[i])
                        } else {
                            None
                        }
                    })
                    .collect(),
            });
        }
    }
    if rows.is_empty() {
        return Err(DrugzError::Input(
            "no retained guides after exclusion".into(),
        ));
    }
    Ok(rows)
}
