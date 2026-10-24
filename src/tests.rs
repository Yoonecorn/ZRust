use super::*;
use std::{collections::BTreeMap, path::Path};
#[path = "tie_regression.rs"]
mod tie_regression;
fn config() -> AnalysisConfig {
    AnalysisConfig {
        control: vec!["C1".into(), "C2".into()],
        treatment: vec!["T1".into(), "T2".into()],
        half_window_size: 5,
        ..Default::default()
    }
}
fn fixture(name: &str) -> Vec<u8> {
    std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name),
    )
    .unwrap()
}
fn result() -> AnalysisResult {
    analyze(&fixture("synthetic.tsv"), &config()).unwrap()
}
fn snapshot(r: &AnalysisResult) -> BTreeMap<String, Vec<f64>> {
    let mut out = BTreeMap::new();
    for (i, row) in r.rows.iter().enumerate() {
        out.insert(format!("normalized/{}", row.id), r.normalized[i].clone());
        for (rep, c) in r.comparisons.iter().enumerate() {
            out.insert(
                format!("comparison/{rep}/{}", row.id),
                vec![
                    c.raw_control[i],
                    c.raw_treatment[i],
                    c.norm_control[i],
                    c.norm_treatment[i],
                    c.fc[i],
                    c.sd[i],
                    c.z[i],
                ],
            );
        }
    }
    for g in &r.genes {
        out.insert(
            format!("gene/{}", g.name),
            vec![
                g.sum,
                g.count as f64,
                g.pre,
                g.norm,
                g.p[0],
                g.p[1],
                g.fdr[0],
                g.fdr[1],
                g.rank[0] as f64,
                g.rank[1] as f64,
            ],
        );
    }
    out
}
#[test]
fn python_differential_all_stages() {
    let base = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cases = [
        (
            "signed_zero",
            "g2/signed_zero.tsv",
            5,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "fractional",
            "g2/fractional.tsv",
            5,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "paired",
            "synthetic.tsv",
            5,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "unpaired",
            "synthetic.tsv",
            5,
            5,
            1,
            true,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "unpaired_unequal",
            "synthetic.tsv",
            5,
            5,
            1,
            true,
            "C1,C2",
            "T1",
            "",
        ),
        (
            "single_pair",
            "synthetic.tsv",
            5,
            5,
            1,
            false,
            "C1",
            "T1",
            "",
        ),
        ("zeros", "zeros.tsv", 5, 5, 1, false, "C1,C2", "T1,T2", ""),
        (
            "pseudocount_one",
            "zeros.tsv",
            5,
            1,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "pseudocount_zero_finite",
            "synthetic.tsv",
            5,
            0,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "exclude",
            "synthetic.tsv",
            25,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "G00,G19",
        ),
        (
            "small_window",
            "synthetic.tsv",
            3,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "large_window",
            "synthetic.tsv",
            159,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "minobs_eight",
            "synthetic.tsv",
            5,
            5,
            8,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "reordered",
            "synthetic.tsv",
            5,
            5,
            1,
            false,
            "C2,C1",
            "T1,T2",
            "",
        ),
        (
            "finite_ties",
            "g2/finite_ties.tsv",
            80,
            5,
            1,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
        (
            "mixed_minobs",
            "g2/mixed.tsv",
            80,
            5,
            6,
            false,
            "C1,C2",
            "T1,T2",
            "",
        ),
    ];
    let mut report =
        String::from("fixture\tstage\tvalues\tmax_absolute_difference\tmax_tolerance_fraction\n");
    for (name, source, h, p, minobs, unpaired, control, treatment, exclude) in cases {
        let split = |s: &str| s.split(',').map(String::from).collect::<Vec<_>>();
        let c = AnalysisConfig {
            control: split(control),
            treatment: split(treatment),
            exclude: if exclude.is_empty() {
                vec![]
            } else {
                split(exclude)
            },
            pseudocount: p,
            minobs,
            half_window_size: h,
            unpaired,
        };
        let r = analyze(&fixture(source), &c).unwrap();
        let actual = snapshot(&r);
        let text =
            std::fs::read_to_string(base.join(format!("references/unrounded-probes/{name}.tsv")))
                .unwrap();
        let expected: BTreeMap<String, Vec<f64>> = text
            .lines()
            .map(|line| {
                let mut fields = line.split('\t');
                (
                    fields.next().unwrap().into(),
                    fields.map(|s| s.parse().unwrap()).collect(),
                )
            })
            .collect();
        assert_eq!(
            actual.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>(),
            "{name} row and gene membership"
        );
        let mut stages: BTreeMap<String, (usize, f64, f64)> = BTreeMap::new();
        for (key, a) in &actual {
            let b = &expected[key];
            assert_eq!(a.len(), b.len());
            for (column, (&a, &b)) in a.iter().zip(b).enumerate() {
                let stage = if key.starts_with("normalized/") {
                    "normalized"
                } else if key.starts_with("comparison/") {
                    [
                        "raw_control",
                        "raw_treatment",
                        "normalized_control",
                        "normalized_treatment",
                        "fold_change",
                        "local_sd",
                        "guide_z",
                    ][column]
                } else {
                    [
                        "sumZ",
                        "numObs",
                        "preZ",
                        "normZ",
                        "p_synth",
                        "p_supp",
                        "fdr_synth",
                        "fdr_supp",
                        "rank_synth",
                        "rank_supp",
                    ][column]
                };
                // Legacy positional tie ranks differ from accepted UTF-8 tie ordering.
                if name == "finite_ties" && stage.starts_with("rank_") {
                    continue;
                }
                assert!(
                    a.is_finite() && b.is_finite(),
                    "{name}/{key}/{column} finite classification"
                );
                assert_eq!(
                    a == 0.0,
                    b == 0.0,
                    "{name}/{key}/{column} zero classification"
                );
                let tolerance = 1e-10 + 1e-9 * a.abs().max(b.abs());
                let difference = (a - b).abs();
                if stage == "numObs" || stage.starts_with("rank_") {
                    assert_eq!(a, b, "{name}/{key}/{column}");
                }
                assert!(
                    difference <= tolerance,
                    "first divergent {stage}: {name}/{key}/{column}: Rust={a:.17e} Python={b:.17e}, tolerance={tolerance}"
                );
                if stage.starts_with("p_") && a.abs().min(b.abs()) < 1e-12 {
                    assert!(difference <= 1e-12);
                }
                let entry = stages.entry(stage.into()).or_default();
                entry.0 += 1;
                entry.1 = entry.1.max(difference);
                entry.2 = entry.2.max(difference / tolerance);
            }
        }
        for (stage, (n, maximum, fraction)) in stages {
            use std::fmt::Write;
            writeln!(
                report,
                "{name}\t{stage}\t{n}\t{maximum:.17e}\t{fraction:.17e}"
            )
            .unwrap();
        }
        for (rep, c) in r.comparisons.iter().enumerate() {
            assert_eq!(c.order.len(), r.rows.len());
            // G3R supersedes stable local ties only; historical G2 probes stay immutable.
            let orders = include_str!(
                "../tie-order-requalification/fixtures/G2-local-order-replacements.tsv"
            );
            let line = orders
                .lines()
                .find(|line| line.starts_with(&format!("{name}\t{rep}\t")))
                .unwrap();
            let expected: Vec<_> = line.split('\t').nth(2).unwrap().split(',').collect();
            assert_eq!(
                c.order
                    .iter()
                    .map(|&i| r.rows[i].id.as_str())
                    .collect::<Vec<_>>(),
                expected,
                "{name}/{rep} local order"
            );
        }
        for pair in r.genes.windows(2) {
            assert!(
                pair[0].norm < pair[1].norm
                    || (pair[0].norm == pair[1].norm && pair[0].name < pair[1].name)
            );
        }
    }
    std::fs::write(base.join("target/G2-differential-results.tsv"), report).unwrap();
}
#[test]
fn strict_tsv_and_validation_precedence() {
    let c = config();
    let mut invalid = c.clone();
    invalid.minobs = 0;
    assert!(matches!(
        analyze(b"\xff", &invalid),
        Err(DrugzError::Parameter(_))
    ));
    invalid = c.clone();
    invalid.control.push("C1".into());
    assert!(matches!(
        analyze(b"\xff", &invalid),
        Err(DrugzError::Design(_))
    ));
    assert!(matches!(analyze(b"\xff", &c), Err(DrugzError::Input(_))));
    for input in [
        "id\tgene\tC1\tC1\tT1\tT2\nx\tG\t1\t1\t2\t2\n",
        "id\tgene\tC1\tC2\tT1\tT2\nx\tG\t1\n",
        "id\tgene\tC1\tC2\tT1\tT2\n\tG\t1\t1\t1\t1\n",
        "id\tgene\tC1\tC2\tT1\tT2\nx\t\t1\t1\t1\t1\n",
        "id\tgene\tC1\tC2\tT1\tT2\nx\tG\t1\t1\t1\t1\nx\tG\t1\t1\t1\t1\n",
        "id\tgene\tC1\t\tT1\tT2\n",
    ] {
        assert!(
            matches!(analyze(input.as_bytes(), &c), Err(DrugzError::Input(_))),
            "{input}"
        );
    }
    let mut missing = c.clone();
    missing.control = vec![" MISSING".into(), "C2".into()];
    let malformed = b"id\tgene\tC1\tC2\tT1\tT2\nx\tG\tBAD\t1\t1\t1\n";
    assert!(matches!(
        analyze(malformed, &missing),
        Err(DrugzError::Design(_))
    ));
    assert!(
        analyze(malformed, &c)
            .unwrap_err()
            .to_string()
            .contains("row 2, column 3, sample \"C1\", gene \"G\", guide \"x\"")
    );
    let mut overlap = c.clone();
    overlap.treatment[0] = "C1".into();
    assert!(matches!(overlap.validate(), Err(DrugzError::Design(_))));
    overlap = c.clone();
    overlap.control[0] = "".into();
    assert!(overlap.validate().is_err());
    overlap = c.clone();
    overlap.treatment.pop();
    assert!(overlap.validate().is_err());
    overlap.unpaired = true;
    assert!(overlap.validate().is_ok());
}
#[test]
fn string_identifiers_and_selected_sample_order() {
    let c = config();
    let rows=input::parse(b"id\tgene\tT2\tC2\tT1\tC1\tignored\n001\t0001\t4\t2\t3\t1\tBAD\n002\tNA\t8\t6\t7\t5\tNaN\n",&c).unwrap();
    assert_eq!(rows[0].id, "001");
    assert_eq!(rows[0].gene, "0001");
    assert_eq!(rows[1].gene, "NA");
    assert_eq!(rows[0].counts, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(rows[1].position, 1);
}
#[test]
fn exclusion_normalization_and_pseudocount_order() {
    let c = config();
    let rows = input::parse(&fixture("synthetic.tsv"), &c).unwrap();
    let names = c.control.iter().chain(&c.treatment).collect::<Vec<_>>();
    let norm = normalization::normalize(&rows, &names).unwrap();
    let total = rows.iter().map(|r| r.counts[0]).sum::<f64>();
    assert_eq!(norm[0][0], 1e7 * rows[0].counts[0] / total);
    let fc = fold_change::calculate(&rows, &norm, &c, 0).unwrap();
    assert_eq!(fc.fc[0], ((norm[0][2] + 5.0) / (norm[0][0] + 5.0)).log2());
    let mut unpaired = c.clone();
    unpaired.unpaired = true;
    let u = fold_change::calculate(&rows, &norm, &unpaired, 0).unwrap();
    let expected =
        (((norm[0][2] + norm[0][3]) / 2.0 + 5.0) / ((norm[0][0] + norm[0][1]) / 2.0 + 5.0)).log2();
    assert_eq!(u.fc[0], expected);
    let mut excluded = c.clone();
    excluded.exclude = vec!["G00".into(), "G19".into()];
    excluded.half_window_size = 25;
    let r = analyze(&fixture("synthetic.tsv"), &excluded).unwrap();
    assert_eq!(r.rows.len(), 72);
    assert_eq!(r.rows[0].position, 4);
    assert!(!r.genes.iter().any(|g| excluded.exclude.contains(&g.name)));
    assert_ne!(r.normalized[0][0], norm[4][0]);
}
fn comparison(n: usize) -> fold_change::Comparison {
    let fc = (0..n).map(|i| (i as f64).sin() + i as f64 / 10.0).collect();
    fold_change::Comparison {
        order: (0..n).collect(),
        raw_control: vec![],
        raw_treatment: vec![],
        norm_control: vec![],
        norm_treatment: vec![],
        fc,
        sd: vec![],
        z: vec![],
    }
}
#[test]
fn exact_step25_and_final_tail_overwrite() {
    let mut c = comparison(80);
    variance::local(&mut c, 5).unwrap();
    let first = variance::sample_sd(&c.fc[..10]);
    let second = variance::sample_sd(&c.fc[25..35]).max(first);
    let third = variance::sample_sd(&c.fc[50..60]).max(second);
    assert!(c.sd[..30].iter().all(|&v| v == first));
    assert!(c.sd[30..55].iter().all(|&v| v == second));
    assert!(c.sd[55..].iter().all(|&v| v == third));
    // H=20,N=80 runs 20,45,70; final tail starts at60 and reads59,
    // overwriting the scale assigned at70 even when it is larger.
    let mut c = comparison(80);
    c.fc[75] = 1000.0;
    variance::local(&mut c, 20).unwrap();
    assert!(c.sd[60..].iter().all(|&v| v == c.sd[59]));
    assert!(variance::sample_sd(&c.fc[50..80]) > c.sd[70]);
    for (n, h) in [
        (80, 5),
        (80, 3),
        (30, 20),
        (30, 30),
        (30, 31),
        (30, 59),
        (2, 2),
    ] {
        assert!(variance::local(&mut comparison(n), h).is_ok(), "{n},{h}");
    }
    for (n, h) in [(72, 5), (80, 1), (30, 60), (2, 3), (1, 1)] {
        assert!(variance::local(&mut comparison(n), h).is_err(), "{n},{h}");
    }
    assert!((variance::sample_sd(&[1.0, 2.0, 3.0]) - 1.0).abs() < 1e-15);
}
#[test]
fn zscores_nonzero_counts_and_functional_minobs() {
    let r = result();
    for c in &r.comparisons {
        for i in 0..r.rows.len() {
            assert_eq!(c.z[i], c.fc[i] / c.sd[i]);
        }
    }
    assert!(r.genes.iter().all(|g| g.count == 8));
    let mut c = config();
    c.minobs = 8;
    assert!(analyze(&fixture("synthetic.tsv"), &c).is_ok());
    c.minobs = 9;
    assert!(matches!(
        analyze(&fixture("synthetic.tsv"), &c),
        Err(DrugzError::Numeric(_))
    ));
    let mut r = result();
    let i = r.comparisons[0].order[0];
    let name = r.rows[i].gene.clone();
    r.comparisons[0].z[i] = 0.0;
    r.comparisons[1].z[i] = -0.0;
    let genes = gene::aggregate(&r.rows, &r.comparisons, 1).unwrap();
    assert_eq!(genes.iter().find(|g| g.name == name).unwrap().count, 6);
    let mixed = AnalysisConfig {
        minobs: 6,
        half_window_size: 80,
        ..config()
    };
    let r = analyze(&fixture("g2/mixed.tsv"), &mixed).unwrap();
    assert!(!r.genes.iter().any(|g| g.name == "SMALL"));
    assert!(r.genes.iter().any(|g| g.name == "LARGE"));
}
#[test]
fn population_standardization_separate_fdr_and_ties() {
    let r = result();
    let mean = r.genes.iter().map(|g| g.norm).sum::<f64>() / r.genes.len() as f64;
    let variance = r.genes.iter().map(|g| g.norm * g.norm).sum::<f64>() / r.genes.len() as f64;
    assert!(mean.abs() < 1e-14);
    assert!((variance - 1.0).abs() < 1e-14);
    for direction in 0..2 {
        let mut order: Vec<_> = r.genes.iter().collect();
        order.sort_by_key(|g| g.rank[direction]);
        let mut minimum = f64::INFINITY;
        for (position, g) in order.iter().enumerate().rev() {
            minimum = minimum.min(g.p[direction] / ((position + 1) as f64 / order.len() as f64));
            assert_eq!(g.fdr[direction], minimum);
        }
    }
    let make = |name: &str, pre: f64| gene::Gene {
        name: name.into(),
        sum: pre,
        count: 1,
        pre,
        norm: 0.0,
        p: [0.0; 2],
        fdr: [0.0; 2],
        rank: [0; 2],
    };
    let mut genes = vec![make("é", 1.0), make("z", 1.0), make("a", -2.0)];
    statistics::rank(&mut genes).unwrap();
    assert_eq!(
        genes.iter().map(|g| g.name.as_str()).collect::<Vec<_>>(),
        vec!["a", "z", "é"]
    );
    assert_eq!(genes[1].rank, [2, 1]);
    assert_eq!(genes[2].rank, [3, 2]);
    let mut equal = vec![make("a", 1.0), make("b", 1.0)];
    assert!(statistics::rank(&mut equal).is_err());
}
#[test]
fn accurate_normal_tails() {
    for (z, p) in [
        (0.0, 0.5),
        (1.0, 0.15865525393145707),
        (8.0, 6.22096057427174e-16),
        (10.0, 7.61985302416047e-24),
        (38.0, 2.88542835e-316),
    ] {
        let a = statistics::sf(z);
        assert!((a - p).abs() <= p * 2e-9 + 1e-323, "{z} {a} {p}");
        assert!(statistics::sf(-z) >= 0.5);
    }
    assert_eq!(statistics::sf(40.0), 0.0);
    assert_eq!(statistics::sf(-40.0), 1.0);
}
#[test]
fn typed_failures_and_no_nonfinite_successes() {
    for name in [
        "nonfinite.tsv",
        "negative.tsv",
        "nonnumeric.tsv",
        "allzero.tsv",
        "missing_gene.tsv",
        "malformed.tsv",
        "single.tsv",
        "ties.tsv",
    ] {
        assert!(analyze(&fixture(name), &config()).is_err(), "{name}");
    }
    for value in ["NaN", "inf", "-inf", "-1"] {
        let text = format!("id\tgene\tC1\tC2\tT1\tT2\nx\tG\t{value}\t1\t2\t3\n");
        assert!(matches!(
            analyze(text.as_bytes(), &config()),
            Err(DrugzError::Input(_))
        ));
    }
    let mut c = config();
    c.pseudocount = 0;
    assert!(
        analyze(&fixture("zeros.tsv"), &c)
            .unwrap_err()
            .to_string()
            .contains("logarithm")
    );
    c = config();
    c.exclude = (0..20).map(|i| format!("G{i:02}")).collect();
    assert!(matches!(
        analyze(&fixture("synthetic.tsv"), &c),
        Err(DrugzError::Input(_))
    ));
    let mut rows = input::parse(&fixture("synthetic.tsv"), &config()).unwrap();
    rows[0].counts[0] = 1e308;
    let c = config();
    let names = c.control.iter().chain(&c.treatment).collect::<Vec<_>>();
    assert!(matches!(
        normalization::normalize(&rows, &names),
        Err(DrugzError::Numeric(_))
    ));
    let mut c = config();
    c.half_window_size = usize::MAX;
    assert!(analyze(&fixture("synthetic.tsv"), &c).is_err());
}
#[test]
fn output_schemas_and_legacy_rounding() {
    let r = result();
    let genes = r.gene_output();
    assert!(genes.ends_with('\n'));
    assert_eq!(genes.lines().count(), 21);
    assert_eq!(genes.lines().next().unwrap().split('\t').count(), 10);
    let expected = std::fs::read_to_string("references/legacy-oracles/paired/genes.tsv").unwrap();
    assert_eq!(genes, expected);
    assert_eq!(
        r.guide_output(),
        std::fs::read_to_string("references/legacy-oracles/paired/guides.tsv").unwrap()
    );
    let mut c = config();
    c.unpaired = true;
    let r = analyze(&fixture("synthetic.tsv"), &c).unwrap();
    let guides = r.guide_output();
    assert!(guides.starts_with("guide_id\tGENE\traw_control_mean\traw_treatment_mean\tnormalized_control_mean\tnormalized_treatment_mean\teb_std\tzscore\tmean_fc\n"));
    assert_eq!(guides.lines().count(), 81);
    assert!(guides.lines().all(|s| s.split('\t').count() == 9));
}
#[test]
fn atomic_output_success_preparation_failure_and_commit_limit() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    let dir = tempfile::tempdir_in(root).unwrap();
    let gene = dir.path().join("genes.tsv");
    let guide = dir.path().join("guides.tsv");
    let r = result();
    std::fs::write(&gene, "old gene").unwrap();
    std::fs::write(&guide, "old guide").unwrap();
    assert!(
        r.publish(&gene, Some(&dir.path().join("missing/guides.tsv")))
            .is_err()
    );
    assert_eq!(std::fs::read_to_string(&gene).unwrap(), "old gene");
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
    assert!(r.publish(&gene, Some(&gene)).is_err());
    r.publish(&gene, Some(&guide)).unwrap();
    assert_eq!(std::fs::read_to_string(&gene).unwrap(), r.gene_output());
    assert_eq!(std::fs::read_to_string(&guide).unwrap(), r.guide_output());
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
    let blocked = dir.path().join("blocked");
    std::fs::create_dir(&blocked).unwrap();
    std::fs::write(&gene, "old gene").unwrap();
    assert!(r.publish(&gene, Some(&blocked)).is_err());
    assert_eq!(std::fs::read_to_string(&gene).unwrap(), "old gene");
    std::fs::write(&guide, "old guide").unwrap();
    assert!(r.publish(&blocked, Some(&guide)).is_err());
    assert_eq!(std::fs::read_to_string(&guide).unwrap(), r.guide_output());
    assert!(blocked.is_dir());
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 3);
}
#[test]
fn signed_zero_numeric_ties_and_identical_fractional_genes() {
    let c = config();
    let rows = input::parse(
        b"id\tgene\tC1\tC2\tT1\tT2\na\tA\t-0.0\t1\t1\t1\nb\tB\t0.0\t1\t1\t1\n",
        &c,
    )
    .unwrap();
    let comp = fold_change::calculate(
        &rows,
        &[vec![0.0, 1.0, 1.0, 1.0], vec![0.0, 1.0, 1.0, 1.0]],
        &c,
        0,
    )
    .unwrap();
    assert_eq!(comp.order, vec![0, 1]);
    let mut genes: Vec<_> = ["a", "b", "c"]
        .into_iter()
        .map(|name| gene::Gene {
            name: name.into(),
            sum: 0.1,
            count: 1,
            pre: 0.1,
            norm: 0.0,
            p: [0.0; 2],
            fdr: [0.0; 2],
            rank: [0; 2],
        })
        .collect();
    assert!(
        statistics::rank(&mut genes)
            .unwrap_err()
            .to_string()
            .contains("identical")
    );
}
#[test]
fn reserved_sample_labels_preserve_gene_strings_and_all_output_columns() {
    let original = String::from_utf8(fixture("synthetic.tsv")).unwrap();
    for name in ["GENE", "fc_0", "eb_std_0", "zscore_fc_0"] {
        let text = original.replacen("C1", name, 1);
        let mut c = config();
        c.control[0] = name.into();
        let r = analyze(text.as_bytes(), &c).unwrap();
        assert_eq!(r.gene_output(), result().gene_output());
        assert_eq!(
            r.guide_output().lines().next().unwrap().split('\t').count(),
            10
        );
    }
}
#[test]
fn python_extreme_normal_tails() {
    for line in std::fs::read_to_string("references/unrounded-probes/normal_tails.tsv")
        .unwrap()
        .lines()
    {
        let (z, p) = line.split_once('\t').unwrap();
        let z: f64 = z.parse().unwrap();
        let p: f64 = p.parse().unwrap();
        let a = statistics::sf(z);
        let delta = (a - p).abs();
        // SciPy 1.10.1 underflows at z=38; the accepted accurate SF must not
        // clamp a representable probability to reproduce library underflow.
        assert!(a.is_finite() && p.is_finite());
        if a != 0.0 && p == 0.0 {
            assert_eq!(z, 38.0);
            assert!(a < 1e-300);
        }
        assert!(delta <= 1e-10 + 1e-9 * a.abs().max(p.abs()));
        if p < 1e-12 {
            assert!(delta <= 1e-12);
        }
        if p > 0.0 && p < 0.5 {
            assert!(delta <= p * 2e-9 + 1e-323);
        }
    }
}
