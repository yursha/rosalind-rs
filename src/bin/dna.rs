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
    /// Compute the Hamming distance pairwise across sequences
    Hamming,
    /// Find all 1-based locations of a motif in the DNA sequence
    Motif,
    /// Calculate the consensus string for a set of equal-length DNA sequences
    Consensus,
}

/// A unified internal representation of a sequence entry
struct SequenceRecord {
    id: String,
    sequence: DnaSequence,
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

    let records: Vec<SequenceRecord> = match format {
        InputFormat::Raw => {
            let mut input_str = String::new();
            reader
                .read_to_string(&mut input_str)
                .expect("Failed to read raw input data");

            input_str
                .lines()
                .map(|l| l.trim())
                .filter(|l| !l.is_empty())
                .enumerate()
                .map(|(idx, line)| {
                    let sequence = line.parse().expect("Failed to parse raw DNA sequence");
                    SequenceRecord {
                        id: format!("Line {}", idx + 1),
                        sequence,
                    }
                })
                .collect()
        }
        InputFormat::Fasta => {
            let fasta_parser = FastaReader::new(reader);
            fasta_parser
                .map(|result| {
                    let record = result.expect("Error parsing FASTA stream");
                    let sequence = record
                        .sequence
                        .parse()
                        .expect("Malformed DNA found in FASTA record");
                    SequenceRecord {
                        id: record.id,
                        sequence,
                    }
                })
                .collect()
        }
    };

    if records.is_empty() {
        eprintln!("Warning: No DNA sequences found to process.");
        return;
    }

    match cli.command {
        Commands::Consensus => {
            let sequences: Vec<DnaSequence> = records.iter().map(|r| r.sequence.clone()).collect();
            match DnaSequence::consensus(&sequences) {
                Ok(result) => {
                    println!("{}", result.consensus);

                    // Flatten the profile into four vectors
                    let mut a_counts = Vec::new();
                    let mut c_counts = Vec::new();
                    let mut g_counts = Vec::new();
                    let mut t_counts = Vec::new();

                    for p in &result.profile {
                        a_counts.push(p.a.to_string());
                        c_counts.push(p.c.to_string());
                        g_counts.push(p.g.to_string());
                        t_counts.push(p.t.to_string());
                    }

                    println!("A: {}", a_counts.join(" "));
                    println!("C: {}", c_counts.join(" "));
                    println!("G: {}", g_counts.join(" "));
                    println!("T: {}", t_counts.join(" "));
                }
                Err(e) => eprintln!("Error calculating consensus: {}", e),
            }
        }
        Commands::Motif => {
            if records.len() < 2 {
                eprintln!(
                    "Error: Motif command requires at least two sequences in the input file."
                );
                std::process::exit(1);
            }
            let haystack = &records[0];
            let motif = &records[1];

            let locations = haystack.sequence.find_motif(&motif.sequence);
            let loc_str: Vec<String> = locations.iter().map(|p| p.to_string()).collect();
            println!("{}", loc_str.join(" "));
        }
        Commands::Hamming => {
            // Process sequences pairwise in sequential blocks of two
            for chunk in records.chunks(2) {
                match chunk {
                    [rec1, rec2] => match rec1.sequence.hamming_distance(&rec2.sequence) {
                        Ok(dist) => println!("Dist[{} <-> {}]: {}", rec1.id, rec2.id, dist),
                        Err(e) => {
                            eprintln!("Error comparing {} and {}: {}", rec1.id, rec2.id, e);
                            std::process::exit(1);
                        }
                    },
                    [lone_rec] => {
                        eprintln!(
                            "Warning: Sequence '{}' is unpaired and was skipped for Hamming distance calculation.",
                            lone_rec.id
                        );
                    }
                    _ => unreachable!(),
                }
            }
        }
        _ => {
            for rec in &records {
                process_sequence(&rec.sequence, &cli.command, &rec.id);
            }
        }
    }
}

/// Dispatches operations that run independently on a single DNA sequence
fn process_sequence(sequence: &DnaSequence, command: &Commands, id: &str) {
    match command {
        Commands::Count => {
            println!("[{}] {}", id, sequence.count_bases());
        }
        Commands::Transcribe => {
            println!("[{}] {}", id, sequence.transcribe());
        }
        Commands::Revcomp => {
            println!("[{}] {}", id, sequence.reverse_complement());
        }
        Commands::GcContent => {
            println!("[{}] {:.5}%", id, sequence.gc_content());
        }
        Commands::Hamming => unreachable!(),
        Commands::Motif => unreachable!(),
        Commands::Consensus => unreachable!(),
    }
}
