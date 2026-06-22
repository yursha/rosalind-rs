use rosalind_rs::simulation::dynamics::AgeStructuredModel;
use std::env;
use std::fs;
use std::process;
use std::str::FromStr;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <input_file_path>", args[0]);
        process::exit(1);
    }

    let file_path = &args[1];
    let content = match fs::read_to_string(file_path) {
        Ok(c) => c,
        Err(err) => {
            eprintln!("Error reading file '{}': {}", file_path, err);
            process::exit(1);
        }
    };

    let tokens: Vec<&str> = content.split_whitespace().collect();
    if tokens.len() < 4 {
        eprintln!(
            "Error: Input file must contain exactly four integers:\n\
             <elapsed_intervals> <fecundity> <initial_cohort> <lifespan>"
        );
        process::exit(1);
    }

    let elapsed_intervals: u32 = parse_num(tokens[0], "elapsed time intervals");
    let fecundity: u64 = parse_num(tokens[1], "fecundity");
    let initial_cohort: u128 = parse_num(tokens[2], "initial cohort");
    let lifespan: usize = parse_num(tokens[3], "lifespan");

    let model = AgeStructuredModel::new(fecundity, 1.0).with_lifespan(lifespan);

    let count = model.project(initial_cohort, elapsed_intervals);

    println!("{}", count);
}

fn parse_num<T: FromStr>(token: &str, label: &str) -> T {
    match token.parse::<T>() {
        Ok(count) => count,
        Err(_) => {
            eprintln!("Error: Failed to parse {} from token '{}'", label, token);
            process::exit(1);
        }
    }
}
