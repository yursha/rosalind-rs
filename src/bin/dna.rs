use clap::{Parser, Subcommand, ValueEnum};
use rosalind_rs::io::fasta::FastaReader;
use rosalind_rs::sequence::dna::DnaSequence;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about = "DNA Sequence manipulation utility")]
struct Cli {
    /// Path to the input file
    #[arg(short, long, value_name = "FILE")]
    input_file: PathBuf,

    /// Explicitly set the input format [default: auto-detect from extension]
    #[arg(short, long, value_enum)]
    format: Option<InputFormat>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum InputFormat {
    Raw,
    Fasta,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Count distinct nucleotides in the DNA sequence
    Count,
    /// Transcribe the DNA sequence into an RNA sequence
    Transcribe,
    /// Compute the reverse complement of the DNA sequence
    Revcomp,
    /// Compute the GC-content percentage of the DNA sequence
    GcContent,
}

fn main() {
    let cli = Cli::parse();

    // Determine the file format (explicit flag overrides extension auto-detection)
    let format = cli.format.unwrap_or_else(|| {
        match cli.input_file.extension().and_then(|ext| ext.to_str()) {
            Some("fasta") | Some("fa") | Some("fna") => InputFormat::Fasta,
            _ => InputFormat::Raw,
        }
    });

    let file = File::open(&cli.input_file)
        .unwrap_or_else(|e| panic!("Failed to open file {:?}: {}", cli.input_file, e));
    let mut reader = BufReader::new(file);

    match format {
        InputFormat::Raw => {
            let mut input_str = String::new();
            reader
                .read_to_string(&mut input_str)
                .expect("Failed to read raw input data");
            let sequence: DnaSequence = input_str
                .trim()
                .parse()
                .expect("Failed to parse DNA sequence");

            execute_command_on_sequence(&sequence, &cli.command, None);
        }
        InputFormat::Fasta => {
            let fasta_parser = FastaReader::new(reader);

            for result in fasta_parser {
                let record = result.expect("Error encountered while parsing FASTA stream");
                let sequence: DnaSequence = record
                    .sequence
                    .parse()
                    .expect("Malformed DNA found in FASTA record payload");

                // Print a label for clarity when handling multiple entries
                execute_command_on_sequence(&sequence, &cli.command, Some(&record.id));
            }
        }
    }
}

/// Helper to handle printing outputs cleanly with optional record identifiers
fn execute_command_on_sequence(sequence: &DnaSequence, command: &Commands, label: Option<&str>) {
    let prefix = label.map(|id| format!("[{}] ", id)).unwrap_or_default();

    match command {
        Commands::Count => {
            println!("{}{}", prefix, sequence.count_bases());
        }
        Commands::Transcribe => {
            println!("{}{}", prefix, sequence.transcribe());
        }
        Commands::Revcomp => {
            println!("{}{}", prefix, sequence.reverse_complement());
        }
        Commands::GcContent => {
            println!("{}{:.5}%", prefix, sequence.gc_content());
        }
    }
}
