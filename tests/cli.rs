use std::{
    io::Write,
    path::Path,
    process::{Command, Output, Stdio},
};
fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_zrust"))
        .args(args)
        .output()
        .unwrap()
}
fn base() -> Vec<&'static str> {
    vec![
        "-i",
        "tests/fixtures/synthetic.tsv",
        "-c",
        "C1,C2",
        "-x",
        "T1,T2",
        "--half_window_size",
        "5",
    ]
}
#[test]
fn help_syntax_and_types_do_not_touch_output() {
    for help in ["-h", "--help"] {
        let r = run(&[help]);
        assert!(r.status.success());
        let text = String::from_utf8(r.stdout).unwrap();
        for flag in [
            "-i",
            "-o",
            "-f",
            "-c",
            "-x",
            "-r",
            "-p",
            "-I",
            "--minobs",
            "--half_window_size",
            "-q",
            "-unpaired",
        ] {
            assert!(text.contains(flag));
        }
    }
    let dir = tempfile::tempdir_in("target").unwrap();
    let out = dir.path().join("genes");
    std::fs::write(&out, "sentinel").unwrap();
    for extra in [
        vec!["-p", "0.5"],
        vec!["-p", "-1"],
        vec!["--minobs", "0"],
        vec!["--half_window_size", "0"],
        vec!["-I", "1"],
        vec!["--unknown"],
        vec!["-i"],
        vec!["-f", "-q"],
    ] {
        let mut args = base();
        args.extend(["-o", out.to_str().unwrap()]);
        args.extend(extra);
        let r = run(&args);
        assert_eq!(
            r.status.code(),
            Some(2),
            "{:?} {}",
            args,
            String::from_utf8_lossy(&r.stderr)
        );
        assert_eq!(std::fs::read_to_string(&out).unwrap(), "sentinel");
    }
    assert_eq!(
        run(&["-i", "-q", "-c", "C1", "-x", "T1"]).status.code(),
        Some(2)
    );
    assert_eq!(run(&[]).status.code(), Some(2));
}
#[test]
fn all_finite_cli_cases_status_and_repeated_bytes() {
    let cases = [
        ("synthetic.tsv", "C1,C2", "T1,T2", 5, false),
        ("synthetic.tsv", "C1,C2", "T1", 5, true),
        ("zeros.tsv", "C1,C2", "T1,T2", 5, false),
        ("g2/fractional.tsv", "C1,C2", "T1,T2", 5, false),
        ("g2/finite_ties.tsv", "C1,C2", "T1,T2", 80, false),
    ];
    for (file, control, treatment, h, unpaired) in cases {
        let input = Path::new("tests/fixtures").join(file);
        let hs = h.to_string();
        let mut args = vec![
            "-i",
            input.to_str().unwrap(),
            "-c",
            control,
            "-x",
            treatment,
            "--half_window_size",
            &hs,
            "-q",
        ];
        if unpaired {
            args.push("-unpaired");
        }
        let first = run(&args);
        assert_eq!(
            first.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(first.stderr.is_empty());
        for _ in 0..2 {
            let repeated = run(&args);
            assert_eq!(repeated.status.code(), Some(0));
            assert_eq!(repeated.stdout, first.stdout);
            assert!(repeated.stderr.is_empty());
        }
    }
    let normal = run(&base());
    assert!(String::from_utf8_lossy(&normal.stderr).contains("completed"));
    assert_eq!(
        normal.stdout,
        std::fs::read("references/legacy-oracles/paired/genes.tsv").unwrap()
    );
    let mut args = base();
    args.extend(["--minobs", "8"]);
    assert_eq!(run(&args).stdout, normal.stdout);
    args.pop();
    args.push("9");
    assert_eq!(run(&args).status.code(), Some(1));
}
#[test]
fn stdin_quiet_errors_and_validation_preserve_files() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_zrust"))
        .args([
            "-c",
            "C1,C2",
            "-x",
            "T1,T2",
            "--half_window_size",
            "5",
            "-q",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&std::fs::read("tests/fixtures/synthetic.tsv").unwrap())
        .unwrap();
    let r = child.wait_with_output().unwrap();
    assert!(r.status.success());
    assert!(r.stderr.is_empty());
    assert_eq!(
        r.stdout,
        std::fs::read("references/legacy-oracles/paired/genes.tsv").unwrap()
    );
    let dir = tempfile::tempdir_in("target").unwrap();
    let gene = dir.path().join("gene");
    let guide = dir.path().join("guide");
    std::fs::write(&gene, "old gene").unwrap();
    std::fs::write(&guide, "old guide").unwrap();
    for name in [
        "negative.tsv",
        "nonfinite.tsv",
        "nonnumeric.tsv",
        "allzero.tsv",
        "malformed.tsv",
        "ties.tsv",
        "single.tsv",
    ] {
        let path = Path::new("tests/fixtures").join(name);
        let args = [
            "-i",
            path.to_str().unwrap(),
            "-c",
            "C1,C2",
            "-x",
            "T1,T2",
            "--half_window_size",
            "5",
            "-o",
            gene.to_str().unwrap(),
            "-f",
            guide.to_str().unwrap(),
            "-q",
        ];
        let r = run(&args);
        assert_eq!(r.status.code(), Some(1));
        assert!(!r.stderr.is_empty());
        assert!(r.stdout.is_empty());
        assert_eq!(std::fs::read_to_string(&gene).unwrap(), "old gene");
        assert_eq!(std::fs::read_to_string(&guide).unwrap(), "old guide");
    }
    for (controls, treatments) in [
        ("C1,C1", "T1,T2"),
        ("C1,C2", "C2,T2"),
        ("C1,C2", "T1"),
        (" C1,C2", "T1,T2"),
        ("C1,", "T1,T2"),
    ] {
        let r = run(&[
            "-i",
            "tests/fixtures/synthetic.tsv",
            "-c",
            controls,
            "-x",
            treatments,
            "--half_window_size",
            "5",
            "-q",
        ]);
        assert_eq!(r.status.code(), Some(1));
        assert!(!r.stderr.is_empty());
    }
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
}
#[test]
fn cli_file_publication_and_unpaired_schema() {
    let dir = tempfile::tempdir_in("target").unwrap();
    let gene = dir.path().join("gene");
    let guide = dir.path().join("guide");
    let mut args = base();
    args.extend([
        "-o",
        gene.to_str().unwrap(),
        "-f",
        guide.to_str().unwrap(),
        "-unpaired",
        "-q",
    ]);
    let r = run(&args);
    assert!(r.status.success());
    assert!(r.stdout.is_empty());
    assert!(r.stderr.is_empty());
    let bytes = std::fs::read(&guide).unwrap();
    assert!(bytes.starts_with(b"guide_id\tGENE\traw_control_mean"));
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
}
#[test]
fn fixture_and_reference_integrity() {
    // Research checkouts retain the complete immutable-history audit. Source
    // distributions have no Git history or pinned Python environment; verify
    // their exact retained fixture bytes with the standard-library-only tool.
    let (python, verifier) = if Path::new(".git").exists()
        && Path::new("target/legacy-venv/bin/python").is_file()
        && Path::new("tools/qualification/verify.py").is_file()
    {
        (
            "target/legacy-venv/bin/python",
            "tools/qualification/verify.py",
        )
    } else {
        ("python3", "tools/release/verify_fixtures.py")
    };
    let r = Command::new(python).arg(verifier).output().unwrap();
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
}

#[test]
fn qualified_banner_and_version() {
    for flag in ["-h", "--help"] {
        let output = run(&[flag]);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        let text = String::from_utf8(output.stdout).unwrap();
        assert!(text.starts_with("ZRust finite reference profile\nQualified for source use on macOS ARM64; see README for limits\nUsage: zrust "));
        assert!(text.contains("-V, --version"));
    }
    for flag in ["-V", "--version"] {
        let output = run(&[flag]);
        assert!(output.status.success());
        assert_eq!(output.stdout, b"ZRust 0.1.0\n");
        assert!(output.stderr.is_empty());
    }
}
