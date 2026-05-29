use std::fs;
use std::path::PathBuf;
use clap::Parser;
use rosalind::util::DnaSequence;

#[derive(Parser, Debug)]
#[command(version, about = "Counts distinct nucleotides in a DNA sequence file")]
struct Args {
    #[arg(short, long, value_name = "FILE")]
    input_file: PathBuf,
}

fn main() {
    let args = Args::parse();
    let input = fs::read_to_string(&args.input_file)
        .expect("Failed to read input.txt");

    let sequence: DnaSequence = input.trim().parse()
        .expect("Failed to parse DNA sequence");

    let counts = sequence.count_bases();

    println!("{}", counts);
}
