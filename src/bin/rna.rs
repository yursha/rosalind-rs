use clap::{Parser, Subcommand};
use rosalind_rs::sequence::rna::RnaSequence;
use std::fs;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about = "RNA Sequence manipulation utility")]
struct Cli {
    /// Path to the input file
    #[arg(short, long, value_name = "FILE")]
    input_file: PathBuf,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Translate an RNA sequence into a protein string
    Translate,
}

fn main() {
    let cli = Cli::parse();

    let content = fs::read_to_string(&cli.input_file)
        .unwrap_or_else(|e| panic!("Failed to read file {:?}: {}", cli.input_file, e))
        .trim()
        .to_string();

    let rna: RnaSequence = content.parse().expect("Invalid RNA sequence in file");

    match cli.command {
        Commands::Translate => {
            match rna.translate() {
                Ok(protein) => {
                    // Convert AminoAcid enum vec to a String for output
                    let protein_str: String = protein
                        .iter()
                        .map(|aa| format!("{}", aa))
                        .collect::<Vec<String>>()
                        .join("");
                    println!("{}", protein_str);
                }
                Err(e) => eprintln!("Translation error: {:?}", e),
            }
        }
    }
}
