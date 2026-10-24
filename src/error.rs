use std::fmt;

/// Failures carry a phase and contextual diagnostic. No partial analysis is returned.
#[derive(Debug)]
pub enum DrugzError {
    Parameter(String),
    Design(String),
    Input(String),
    Numeric(String),
    Output(String),
    Io(std::io::Error),
}
impl fmt::Display for DrugzError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parameter(s) => write!(f, "parameter: {s}"),
            Self::Design(s) => write!(f, "design: {s}"),
            Self::Input(s) => write!(f, "input: {s}"),
            Self::Numeric(s) => write!(f, "numerical analysis: {s}"),
            Self::Output(s) => write!(f, "output: {s}"),
            Self::Io(e) => write!(f, "I/O: {e}"),
        }
    }
}
impl std::error::Error for DrugzError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        if let Self::Io(e) = self {
            Some(e)
        } else {
            None
        }
    }
}
impl From<std::io::Error> for DrugzError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
pub(crate) fn finite(value: f64, context: impl Into<String>) -> Result<f64, DrugzError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DrugzError::Numeric(context.into()))
    }
}
