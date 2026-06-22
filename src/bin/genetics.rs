use rosalind_rs::simulation::genetics::{CouplePopulation, IndividualPopulation};
use std::env;
use std::fs;
use std::process;

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

    // Extract the space-separated numbers from the input file
    let tokens: Vec<&str> = content.split_whitespace().collect();

    // Dynamically branch execution based on dataset shape
    match tokens.len() {
        3 => {
            // Path A: Unstructured Individual Mendelian Mating Pool
            let homozygous_dominant = parse_population_count(tokens[0], "homozygous dominant");
            let heterozygous = parse_population_count(tokens[1], "heterozygous");
            let homozygous_recessive = parse_population_count(tokens[2], "homozygous recessive");

            let population =
                IndividualPopulation::new(homozygous_dominant, heterozygous, homozygous_recessive);
            let probability = population.dominant_phenotype_probability();

            // Output formatted to 5 decimal places for probability
            println!("{:.5}", probability);
        }
        6 => {
            // Path B: Pre-paired Couple Expectation Values
            let dom_dom = parse_population_count(tokens[0], "AA-AA");
            let dom_het = parse_population_count(tokens[1], "AA-Aa");
            let dom_rec = parse_population_count(tokens[2], "AA-aa");
            let het_het = parse_population_count(tokens[3], "Aa-Aa");
            let het_rec = parse_population_count(tokens[4], "Aa-aa");
            let rec_rec = parse_population_count(tokens[5], "aa-aa");

            let couples =
                CouplePopulation::new(dom_dom, dom_het, dom_rec, het_het, het_rec, rec_rec);
            let expected_offspring = couples.expected_dominant_offspring();

            // Output clean raw float/integer for expected value
            println!("{}", expected_offspring);
        }
        _ => {
            // Fallthrough error boundary
            eprintln!(
                "Error: Invalid input dataset. File must contain exactly 3 tokens (for individual pools) or 6 tokens (for couple pairings). Found {} tokens.",
                tokens.len()
            );
            process::exit(1);
        }
    }
}

/// Helper function to parse input tokens safely with clear error context
fn parse_population_count(token: &str, label: &str) -> u64 {
    match token.parse::<u64>() {
        Ok(count) => count,
        Err(_) => {
            eprintln!(
                "Error: Failed to parse {} count from token '{}'",
                label, token
            );
            process::exit(1);
        }
    }
}
