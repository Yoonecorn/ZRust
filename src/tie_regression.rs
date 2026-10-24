use super::*;

#[test]
fn paired_integer_keys_do_not_round() {
    for (label, start) in [
        ("integer_precision", 1u64 << 53),
        ("uint_precision", 1u64 << 63),
    ] {
        let mut input = String::from("id\tgene\tC1\tT1\n");
        for i in 0..24 {
            input.push_str(&format!("R{i}\tG{i}\t{}\t100\n", start + i % 3));
        }
        let c = AnalysisConfig {
            control: vec!["C1".into()],
            treatment: vec!["T1".into()],
            ..Default::default()
        };
        let rows = input::parse(input.as_bytes(), &c).unwrap();
        let norm = normalization::normalize(&rows, &[&c.control[0], &c.treatment[0]]).unwrap();
        let result = fold_change::calculate(&rows, &norm, &c, 0).unwrap();
        let expected = include_str!("../tie-order-requalification/fixtures/order-regressions.tsv")
            .lines()
            .find(|l| l.starts_with(&format!("{label}\t")))
            .unwrap();
        let expected: Vec<usize> = expected
            .split('\t')
            .nth(3)
            .unwrap()
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(result.order, expected);
    }
}

#[test]
fn causal_stage_regressions() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (name, h, unpaired) in [
        ("all_equal", 25, false),
        ("two_values", 25, false),
        ("zero_heavy", 50, false),
        ("signed_zero", 20, false),
        ("first_window", 500, false),
        ("step25", 35, false),
        ("final_tail", 500, false),
        ("unpaired", 50, true),
        ("ordinary_no_ties", 50, false),
        ("real_derived", 12, false),
    ] {
        let c = AnalysisConfig {
            half_window_size: h,
            unpaired,
            ..config()
        };
        let r = analyze(&fixture(&format!("g3r/{name}.tsv")), &c).unwrap();
        let base = root.join("tie-order-requalification/fixtures/stage-probes");
        let order = std::fs::read_to_string(base.join(format!("{name}.orders.tsv"))).unwrap();
        let mut expected_orders = vec![Vec::new(); r.comparisons.len()];
        for line in order.lines() {
            let (rep, guide) = line.split_once('\t').unwrap();
            expected_orders[rep.parse::<usize>().unwrap()].push(guide);
        }
        for (rep, comparison) in r.comparisons.iter().enumerate() {
            assert_eq!(
                comparison
                    .order
                    .iter()
                    .map(|&i| r.rows[i].id.as_str())
                    .collect::<Vec<_>>(),
                expected_orders[rep],
                "exact causal guide order {name}/{rep}"
            );
        }
        let actual = snapshot(&r);
        let expected = std::fs::read_to_string(base.join(format!("{name}.tsv"))).unwrap();
        let expected: BTreeMap<String, Vec<f64>> = expected
            .lines()
            .map(|line| {
                let mut fields = line.split('\t');
                (
                    fields.next().unwrap().into(),
                    fields.map(|x| x.parse().unwrap()).collect(),
                )
            })
            .collect();
        assert_eq!(
            actual.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>()
        );
        for (key, values) in actual {
            for (i, (&a, &b)) in values.iter().zip(&expected[&key]).enumerate() {
                let tolerance = 1e-10 + 1e-9 * a.abs().max(b.abs());
                assert!((a - b).abs() <= tolerance, "{name}/{key}/{i}: {a} != {b}");
                if key.starts_with("gene/") && (4..8).contains(&i) && a.abs().min(b.abs()) < 1e-12 {
                    assert!((a - b).abs() <= 1e-12);
                }
            }
        }
        assert_eq!(
            r.gene_output(),
            analyze(&fixture(&format!("g3r/{name}.tsv")), &c)
                .unwrap()
                .gene_output()
        );
        assert_eq!(
            r.guide_output(),
            analyze(&fixture(&format!("g3r/{name}.tsv")), &c)
                .unwrap()
                .guide_output()
        );
    }
}
