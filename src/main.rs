mod cli;
fn main() {
    match cli::parse() {
        Ok(None) => (),
        Ok(Some(args)) => {
            if let Err(error) = cli::run(args) {
                eprintln!("zrust: {error}");
                std::process::exit(1);
            }
        }
        Err(error) => {
            eprintln!("zrust: CLI: {error}");
            std::process::exit(2);
        }
    }
}
