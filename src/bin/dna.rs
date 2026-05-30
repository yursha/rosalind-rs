use clap::{Parser, Subcommand};
use rosalind::dna::DnaSequence;
use std::fs;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about)]
struct Cli {
    /// Path to the input DNA sequence file (available to all subcommands)
    #[arg(short, long, value_name = "FILE")]
    input_file: PathBuf,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Count distinct nucleotides in the DNA sequence
    Count,
    /// Transcribe the DNA sequence into an RNA sequence
    Transcribe,
    /// Compute the reverse complement of the DNA sequence
    Revcomp,
}

fn main() {
    let cli = Cli::parse();

    let input = fs::read_to_string(&cli.input_file).expect("Failed to read input.txt");

    let sequence: DnaSequence = input.trim().parse().expect("Failed to parse DNA sequence");

    match cli.command {
        Commands::Count => {
            let counts = sequence.count_bases();
            println!("{}", counts);
        }
        Commands::Transcribe => {
            let rna = sequence.transcribe();
            println!("{}", rna);
        }
        Commands::Revcomp => {
            let rev_comp = sequence.reverse_complement();
            println!("{}", rev_comp);
        }
    }
}
