use crate::{AnalysisResult, DrugzError};
use std::{fmt::Write as _, io::Write, path::Path};

// Python's .3g chooses fixed notation for rounded exponents -4 through 2.
fn significant(value: f64) -> String {
    let scientific = format!("{value:.2e}");
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("finite scientific number");
    let exponent: i32 = exponent.parse().expect("formatted exponent");
    let trim = |s: String| {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    };
    if !(-4..3).contains(&exponent) {
        format!("{}e{exponent:+03}", trim(mantissa.into()))
    } else {
        trim(format!(
            "{value:.precision$}",
            precision = (2 - exponent) as usize
        ))
    }
}
pub(crate) fn genes(r: &AnalysisResult) -> String {
    let mut s = String::from(
        "GENE\tsumZ\tnumObs\tnormZ\tpval_synth\trank_synth\tfdr_synth\tpval_supp\trank_supp\tfdr_supp\n",
    );
    s = chunks(r.genes.len(), r.threads, 25_000, s, |range, s| {
        gene_rows(r, range, s)
    });
    s
}
fn gene_rows(r: &AnalysisResult, range: std::ops::Range<usize>, s: &mut String) {
    for g in &r.genes[range] {
        writeln!(
            s,
            "{}\t{:.2}\t{}\t{:4.2}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.name,
            g.sum,
            g.count,
            g.norm,
            significant(g.p[0]),
            g.rank[0],
            significant(g.fdr[0]),
            significant(g.p[1]),
            g.rank[1],
            significant(g.fdr[1])
        )
        .expect("string write");
    }
}
pub(crate) fn guides(r: &AnalysisResult) -> String {
    let mut s = String::new();
    if r.config.unpaired {
        s.push_str("guide_id\tGENE\traw_control_mean\traw_treatment_mean\tnormalized_control_mean\tnormalized_treatment_mean\teb_std\tzscore\tmean_fc\n");
    } else {
        s.push_str("\tGENE");
        for (rep, sample) in r.config.control.iter().enumerate() {
            write!(s, "\t{sample}\teb_std_{rep}\tzscore_fc_{rep}\tfc_{rep}").expect("string write");
        }
        s.push('\n');
    }
    s = chunks(r.rows.len(), r.threads, 50_000, s, |range, s| {
        guide_rows(r, range, s)
    });
    s
}
fn guide_rows(r: &AnalysisResult, range: std::ops::Range<usize>, s: &mut String) {
    if r.config.unpaired {
        let c = &r.comparisons[0];
        for &i in &c.order[range] {
            let row = &r.rows[i];
            writeln!(
                s,
                "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                row.id,
                row.gene,
                c.raw_control[i],
                c.raw_treatment[i],
                c.norm_control[i],
                c.norm_treatment[i],
                c.sd[i],
                c.z[i],
                c.fc[i]
            )
            .expect("string write");
        }
    } else {
        for &i in &r.comparisons[0].order[range] {
            write!(s, "{}\t{}", r.rows[i].id, r.rows[i].gene).expect("string write");
            for (rep, c) in r.comparisons.iter().enumerate() {
                write!(
                    s,
                    "\t{}\t{:.3}\t{:.3}\t{:.3}",
                    if r.rows[i].integer_columns[rep] {
                        if c.raw_control[i] == 0.0 {
                            "0".into()
                        } else {
                            c.raw_control[i].to_string()
                        }
                    } else {
                        format!("{:.3}", c.raw_control[i])
                    },
                    c.sd[i],
                    c.z[i],
                    c.fc[i]
                )
                .expect("string write");
            }
            s.push('\n');
        }
    }
}
// Independent formatting only. Handles and bytes are consumed in source order.
fn chunks(
    len: usize,
    threads: usize,
    threshold: usize,
    mut output: String,
    render: impl Fn(std::ops::Range<usize>, &mut String) + Sync,
) -> String {
    if threads == 1 || len < threshold {
        render(0..len, &mut output);
        return output;
    }
    let workers = threads.min(8).min(len.div_ceil(threshold / 8));
    let width = len.div_ceil(workers);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..len)
            .step_by(width)
            .map(|start| {
                let render = &render;
                std::thread::Builder::new()
                    .spawn_scoped(scope, move || {
                        let mut part = String::new();
                        render(start..(start + width).min(len), &mut part);
                        part
                    })
                    .map_err(|_| {
                        let mut part = String::new();
                        render(start..(start + width).min(len), &mut part);
                        part
                    })
            })
            .collect();
        let parts: Vec<_> = handles
            .into_iter()
            .map(|handle| match handle {
                Ok(handle) => match handle.join() {
                    Ok(part) => part,
                    Err(payload) => std::panic::resume_unwind(payload),
                },
                Err(part) => part,
            })
            .collect();
        output.reserve(parts.iter().map(String::len).sum());
        for part in parts {
            output.push_str(&part);
        }
        output
    })
}

fn prepare(path: &Path, bytes: &[u8]) -> Result<tempfile::NamedTempFile, DrugzError> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut file = tempfile::NamedTempFile::new_in(parent)?;
    file.write_all(bytes)?;
    file.as_file().sync_all()?;
    Ok(file)
}
fn commit(file: tempfile::NamedTempFile, path: &Path) -> Result<(), DrugzError> {
    file.persist(path)
        .map_err(|e| DrugzError::Output(format!("publishing {}: {}", path.display(), e.error)))?;
    Ok(())
}
pub(crate) fn publish(
    r: &AnalysisResult,
    gene: &Path,
    guide: Option<&Path>,
) -> Result<(), DrugzError> {
    if let Some(guide) = guide {
        // Canonical parent + filename also detects aliases through parent symlinks.
        let identity = |p: &Path| -> Result<_, DrugzError> {
            let parent = p
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            Ok(std::fs::canonicalize(parent)?.join(
                p.file_name()
                    .ok_or_else(|| DrugzError::Output("destination requires filename".into()))?,
            ))
        };
        if identity(gene)? == identity(guide)? {
            return Err(DrugzError::Output(
                "gene and guide destinations must differ".into(),
            ));
        }
    }
    let gene_file = prepare(gene, r.gene_output().as_bytes())?;
    let guide_file = guide
        .map(|p| prepare(p, r.guide_output().as_bytes()))
        .transpose()?;
    if let Some((file, path)) = guide_file.zip(guide) {
        commit(file, path)?;
    }
    commit(gene_file, gene)
}
#[cfg(test)]
mod tests {
    use super::significant;
    #[test]
    fn python_significant_format() {
        for (v, expected) in [
            (0.0, "0"),
            (1.0, "1"),
            (0.001234, "0.00123"),
            (0.00001234, "1.23e-05"),
            (0.9999, "1"),
            (1234.0, "1.23e+03"),
        ] {
            assert_eq!(significant(v), expected);
        }
    }
    #[test]
    fn formatting_chunks_preserve_source_order_at_boundaries() {
        for len in [0, 1, 24_999, 25_000, 25_003, 50_003] {
            let render = |range: std::ops::Range<usize>, out: &mut String| {
                // Slow early chunks so completion order differs from output order.
                if range.start == 0 && range.end < len {
                    std::thread::sleep(std::time::Duration::from_millis(2));
                }
                for i in range {
                    use std::fmt::Write;
                    writeln!(out, "β\t{i}").unwrap();
                }
            };
            let expected = super::chunks(len, 1, 25_000, "header\n".into(), render);
            for workers in [2, 4, 8, usize::MAX] {
                assert_eq!(
                    super::chunks(len, workers, 25_000, "header\n".into(), render),
                    expected
                );
            }
        }
    }
}
