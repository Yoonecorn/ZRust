use crate::{AnalysisConfig, DrugzError};
use std::collections::BTreeSet;
pub(crate) fn validate(c: &AnalysisConfig) -> Result<(), DrugzError> {
    if c.half_window_size == 0 || c.minobs == 0 {
        return Err(DrugzError::Parameter(
            "half_window_size and minobs must be >= 1".into(),
        ));
    }
    let mut names = BTreeSet::new();
    for (arm, list) in [("control", &c.control), ("treatment", &c.treatment)] {
        if list.is_empty() {
            return Err(DrugzError::Design(format!("empty {arm} arm")));
        }
        for name in list {
            if name.is_empty() || !names.insert(name) {
                return Err(DrugzError::Design(format!(
                    "empty, duplicate or overlapping {arm} sample {name:?}"
                )));
            }
        }
    }
    if !c.unpaired && c.control.len() != c.treatment.len() {
        return Err(DrugzError::Design(
            "paired arms must have equal nonzero sizes".into(),
        ));
    }
    Ok(())
}
pub(crate) fn resolve(header: &[&str], c: &AnalysisConfig) -> Result<Vec<usize>, DrugzError> {
    c.control
        .iter()
        .chain(&c.treatment)
        .map(|name| {
            header
                .iter()
                .enumerate()
                .skip(2)
                .find(|(_, h)| **h == name)
                .map(|(i, _)| i)
                .ok_or_else(|| DrugzError::Design(format!("missing selected sample {name:?}")))
        })
        .collect()
}
