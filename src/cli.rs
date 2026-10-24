use std::path::PathBuf;
use zrust::{AnalysisConfig, DrugzError};
pub(crate) struct Args {
    pub config: AnalysisConfig,
    pub input: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub guide: Option<PathBuf>,
    pub quiet: bool,
    pub threads: usize,
}
const HELP: &str = "ZRust finite reference profile\nQualified for source use on macOS ARM64; see README for limits\nUsage: zrust -c C1,C2 -x T1,T2 [options]\n  -i FILE                 UTF-8 TSV input (default stdin; - for stdin)\n  -o FILE                 gene output (default stdout; - for stdout)\n  -f FILE                 optional guide output file\n  -c LIST                 literal comma-separated controls (required)\n  -x LIST                 literal comma-separated treatments (required)\n  -r LIST                 excluded genes, before normalization\n  -p INTEGER              normalized-count pseudocount >=0 (default 5)\n  -I INTEGER              guide index; only 0 is supported\n  --minobs INTEGER        nonzero guide observations >=1 (default 1)\n  --half_window_size INT  local half window >=1 (default 500)\n  --threads INTEGER       output-format workers >=1 (default 8; cap 8)\n  -q                      suppress progress, retain errors\n  -unpaired               mean normalized counts within each arm\n  -h, --help              show help\n  -V, --version           show version\n";
pub(crate) fn parse() -> Result<Option<Args>, String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{HELP}");
        return Ok(None);
    }
    if args.iter().any(|a| a == "-V" || a == "--version") {
        println!("ZRust {}", env!("CARGO_PKG_VERSION"));
        return Ok(None);
    }
    let mut result = Args {
        config: AnalysisConfig::default(),
        input: None,
        output: None,
        guide: None,
        quiet: false,
        threads: 8,
    };
    let mut seen = std::collections::BTreeSet::new();
    let mut i = 0;
    while i < args.len() {
        let flag = &args[i];
        if !seen.insert(flag) {
            return Err(format!("duplicate option {flag}"));
        }
        if flag == "-q" {
            result.quiet = true;
            i += 1;
            continue;
        }
        if flag == "-unpaired" {
            result.config.unpaired = true;
            i += 1;
            continue;
        }
        if ![
            "-i",
            "-o",
            "-f",
            "-c",
            "-x",
            "-r",
            "-p",
            "-I",
            "--minobs",
            "--threads",
            "--half_window_size",
        ]
        .contains(&flag.as_str())
        {
            return Err(format!("unknown option {flag}"));
        }
        i += 1;
        let value = args
            .get(i)
            .ok_or_else(|| format!("{flag} requires a value"))?;
        if [
            "-i",
            "-o",
            "-f",
            "-c",
            "-x",
            "-r",
            "-p",
            "-I",
            "--minobs",
            "--threads",
            "--half_window_size",
            "-q",
            "-unpaired",
        ]
        .contains(&value.as_str())
        {
            return Err(format!("{flag} requires a value before option {value}"));
        }
        let integer = || {
            value
                .parse::<i64>()
                .map_err(|_| format!("{flag} requires an integer, got {value:?}"))
        };
        match flag.as_str() {
            "-i" => result.input = (value != "-").then(|| value.into()),
            "-o" => result.output = (value != "-").then(|| value.into()),
            "-f" => result.guide = Some(value.into()),
            "-c" => result.config.control = value.split(',').map(String::from).collect(),
            "-x" => result.config.treatment = value.split(',').map(String::from).collect(),
            "-r" => result.config.exclude = value.split(',').map(String::from).collect(),
            "-p" => {
                result.config.pseudocount =
                    u64::try_from(integer()?).map_err(|_| "-p must be >= 0")?
            }
            "-I" => {
                if integer()? != 0 {
                    return Err("-I: unsupported value; G2 accepts only 0".into());
                }
            }
            "--threads" | "--minobs" | "--half_window_size" => {
                let number = usize::try_from(integer()?)
                    .ok()
                    .filter(|&n| n >= 1)
                    .ok_or_else(|| format!("{flag} must be >= 1"))?;
                if flag == "--threads" {
                    result.threads = number;
                } else if flag == "--minobs" {
                    result.config.minobs = number;
                } else {
                    result.config.half_window_size = number;
                }
            }
            _ => unreachable!(),
        }
        i += 1;
    }
    if !seen.contains(&String::from("-c")) || !seen.contains(&String::from("-x")) {
        return Err("-c and -x are required".into());
    }
    Ok(Some(result))
}
pub(crate) fn run(args: Args) -> Result<(), DrugzError> {
    use std::io::{Read, Write};
    if !args.quiet {
        eprintln!("ZRust: calculating finite reference profile");
    }
    args.config.validate()?;
    let bytes = if let Some(path) = args.input {
        std::fs::read(path)?
    } else {
        let mut bytes = Vec::new();
        std::io::stdin().read_to_end(&mut bytes)?;
        bytes
    };
    let result = zrust::analyze_with_threads(&bytes, &args.config, args.threads)?;
    if let Some(path) = args.output {
        result.publish(&path, args.guide.as_deref())?;
    } else {
        // Prepare stdout bytes before publishing the optional guide file.
        let gene = result.gene_output();
        if let Some(path) = args.guide {
            let parent = path
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(std::path::Path::new("."));
            let mut temp = tempfile::NamedTempFile::new_in(parent)?;
            temp.write_all(result.guide_output().as_bytes())?;
            temp.as_file().sync_all()?;
            temp.persist(path)
                .map_err(|e| DrugzError::Output(e.error.to_string()))?;
        }
        std::io::stdout().lock().write_all(gene.as_bytes())?;
    }
    if !args.quiet {
        eprintln!("ZRust: completed");
    }
    Ok(())
}
