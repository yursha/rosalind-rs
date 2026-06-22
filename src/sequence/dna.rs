use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

use crate::sequence::rna::{RnaBase, RnaSequence};

#[derive(Debug, PartialEq, Eq)]
pub struct InvalidDnaSymbolError(pub char);

#[derive(Debug, PartialEq, Eq)]
pub struct LengthMismatchError {
    pub sequence_len: usize,
    pub other_len: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub struct ConsensusResult {
    pub consensus: DnaSequence,
    pub profile: Vec<DnaBaseCounts>,
}

impl std::fmt::Display for LengthMismatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Hamming distance undefined for unequal lengths: {} vs {}",
            self.sequence_len, self.other_len
        )
    }
}

impl std::error::Error for LengthMismatchError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DnaBase {
    A,
    C,
    G,
    T,
}

impl fmt::Display for DnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            DnaBase::A => 'A',
            DnaBase::C => 'C',
            DnaBase::G => 'G',
            DnaBase::T => 'T',
        };
        write!(f, "{}", c)
    }
}

impl TryFrom<char> for DnaBase {
    type Error = InvalidDnaSymbolError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(DnaBase::A),
            'C' => Ok(DnaBase::C),
            'G' => Ok(DnaBase::G),
            'T' => Ok(DnaBase::T),
            _ => Err(InvalidDnaSymbolError(c)),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaSequence(pub Vec<DnaBase>);

impl FromStr for DnaSequence {
    type Err = InvalidDnaSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bases = s
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' => Ok(DnaBase::A),
                'C' => Ok(DnaBase::C),
                'G' => Ok(DnaBase::G),
                'T' => Ok(DnaBase::T),
                _ => Err(InvalidDnaSymbolError(c)),
            })
            .collect::<Result<Vec<DnaBase>, _>>()?;
        Ok(DnaSequence(bases))
    }
}

impl fmt::Display for DnaSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in &self.0 {
            let c = match base {
                DnaBase::A => 'A',
                DnaBase::C => 'C',
                DnaBase::G => 'G',
                DnaBase::T => 'T',
            };
            write!(f, "{}", c)?;
        }
        Ok(())
    }
}

impl DnaSequence {
    /// Counts the occurrences of each distinct nucleotide in the sequence.
    /// Operates in O(n) time and O(1) auxiliary space.
    pub fn count_bases(&self) -> DnaBaseCounts {
        let mut counts = DnaBaseCounts {
            a: 0,
            c: 0,
            g: 0,
            t: 0,
        };

        for base in &self.0 {
            match base {
                DnaBase::A => counts.a += 1,
                DnaBase::C => counts.c += 1,
                DnaBase::G => counts.g += 1,
                DnaBase::T => counts.t += 1,
            }
        }

        counts
    }

    /// Transcribes the DNA sequence into an RNA sequence by swapping Thymine (T) for Uracil (U).
    /// Operates in O(n) time and returns a typed RnaSequence wrapper.
    pub fn transcribe(&self) -> RnaSequence {
        let rna_bases = self
            .0
            .iter()
            .map(|base| match base {
                DnaBase::A => RnaBase::A,
                DnaBase::C => RnaBase::C,
                DnaBase::G => RnaBase::G,
                DnaBase::T => RnaBase::U,
            })
            .collect();

        RnaSequence(rna_bases)
    }

    /// Returns the reverse complement of the DNA sequence.
    /// Reverses the order of the bases and swaps each nucleotide with its complement (A <-> T, C <-> G).
    /// Operates in O(n) time.
    pub fn reverse_complement(&self) -> DnaSequence {
        let rev_comp_bases = self
            .0
            .iter()
            .rev()
            .map(|base| match base {
                DnaBase::A => DnaBase::T,
                DnaBase::C => DnaBase::G,
                DnaBase::G => DnaBase::C,
                DnaBase::T => DnaBase::A,
            })
            .collect();

        DnaSequence(rev_comp_bases)
    }

    /// Calculates the GC-content of the DNA sequence as a percentage.
    ///
    /// Returns `0.0` if the sequence is empty to prevent a division-by-zero panic.
    pub fn gc_content(&self) -> f64 {
        if self.0.is_empty() {
            return 0.0;
        }

        // Count occurrences of C and G variants within the strongly typed vector
        let gc_count = self
            .0
            .iter()
            .filter(|&&base| matches!(base, DnaBase::C | DnaBase::G))
            .count();

        (gc_count as f64 / self.0.len() as f64) * 100.0
    }

    /// Computes the Hamming distance between this sequence and another.
    ///
    /// Returns a `LengthMismatchError` if the sequences are not of equal length.
    /// Operates in O(n) time and O(1) auxiliary space.
    pub fn hamming_distance(&self, other: &Self) -> Result<usize, LengthMismatchError> {
        if self.0.len() != other.0.len() {
            return Err(LengthMismatchError {
                sequence_len: self.0.len(),
                other_len: other.0.len(),
            });
        }

        // Zip the two iterators together and count positions where elements differ
        let distance = self
            .0
            .iter()
            .zip(other.0.iter())
            .filter(|(b1, b2)| b1 != b2)
            .count();

        Ok(distance)
    }

    /// Returns the 1-based start positions of all occurrences of a motif within the sequence.
    pub fn find_motif(&self, motif: &DnaSequence) -> Vec<usize> {
        let n = self.0.len();
        let m = motif.0.len();
        let mut positions = Vec::new();

        if m == 0 || m > n {
            return positions;
        }

        // Slide the window across the sequence
        for i in 0..=(n - m) {
            // Compare the slice of our sequence with the motif
            if &self.0[i..i + m] == &motif.0[..] {
                // Convert 0-based index to 1-based position
                positions.push(i + 1);
            }
        }

        positions
    }

    /// Calculates the consensus string and the profile matrix for a collection of equal-length DNA sequences.
    ///
    /// NOTE: This implementation is an oversimplification. It uses deterministic tie-breaking
    /// (A > C > G > T) and does not account for biological ambiguity (IUPAC codes) or
    /// base quality scores common in production-grade bioinformatics.
    pub fn consensus(sequences: &[DnaSequence]) -> Result<ConsensusResult, String> {
        if sequences.is_empty() {
            return Err("Cannot compute consensus of empty sequence set".to_string());
        }

        let len = sequences[0].0.len();
        if sequences.iter().any(|s| s.0.len() != len) {
            return Err("All sequences must have the same length".to_string());
        }

        // Initialize profile matrix
        let mut profile = vec![
            DnaBaseCounts {
                a: 0,
                c: 0,
                g: 0,
                t: 0
            };
            len
        ];

        // Fill profile matrix
        for seq in sequences {
            for (i, base) in seq.0.iter().enumerate() {
                match base {
                    DnaBase::A => profile[i].a += 1,
                    DnaBase::C => profile[i].c += 1,
                    DnaBase::G => profile[i].g += 1,
                    DnaBase::T => profile[i].t += 1,
                }
            }
        }

        // Derive consensus string
        let consensus_bases = profile
            .iter()
            .map(|p| {
                if p.a >= p.c && p.a >= p.g && p.a >= p.t {
                    DnaBase::A
                } else if p.c >= p.g && p.c >= p.t {
                    DnaBase::C
                } else if p.g >= p.t {
                    DnaBase::G
                } else {
                    DnaBase::T
                }
            })
            .collect();

        Ok(ConsensusResult {
            consensus: DnaSequence(consensus_bases),
            profile,
        })
    }

    /// Computes the directed overlap graph O_k in O(N*k + E) time.
    pub fn overlap_graph(sequences: &[DnaSequence], k: usize) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        let n = sequences.len();

        if k == 0 || n < 2 {
            return edges;
        }

        // PASS 1: Build the Prefix Index
        // Key: A zero-copy slice of the sequence (&[DnaBase])
        // Value: A list of sequence indices that start with this exact k-mer
        let mut prefix_map: HashMap<&[DnaBase], Vec<usize>> = HashMap::with_capacity(n);

        for (i, seq) in sequences.iter().enumerate() {
            let bases = &seq.0;
            if bases.len() >= k {
                let prefix = &bases[..k];
                prefix_map.entry(prefix).or_default().push(i);
            }
        }

        // PASS 2: Probe the table using Suffixes
        for (i, seq) in sequences.iter().enumerate() {
            let bases = &seq.0;
            if bases.len() < k {
                continue;
            }

            let suffix = &bases[bases.len() - k..];

            // O(1) average lookup time
            if let Some(target_indices) = prefix_map.get(suffix) {
                for &j in target_indices {
                    if i != j {
                        // Forbid self-loops
                        edges.push((i, j));
                    }
                }
            }
        }

        edges
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DnaBaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
}

impl fmt::Display for DnaBaseCounts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {} {}", self.a, self.c, self.g, self.t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_parsing_valid() {
        let input = "AGCTagct";
        let expected = vec![
            DnaBase::A,
            DnaBase::G,
            DnaBase::C,
            DnaBase::T,
            DnaBase::A,
            DnaBase::G,
            DnaBase::C,
            DnaBase::T,
        ];

        let parsed: Vec<DnaBase> = input
            .chars()
            .map(DnaBase::try_from)
            .collect::<Result<_, _>>()
            .expect("Valid DNA characters should parse cleanly");

        assert_eq!(parsed, expected);
    }

    #[test]
    fn test_dna_parsing_invalid() {
        let input = "GATXG";
        let result: Result<Vec<DnaBase>, _> = input.chars().map(DnaBase::try_from).collect();

        assert_eq!(result, Err(InvalidDnaSymbolError('X')));
    }

    #[test]
    fn test_dna_to_rna_transcription() {
        let dna = DnaSequence(vec![DnaBase::G, DnaBase::A, DnaBase::T, DnaBase::C]);
        let expected_rna = RnaSequence(vec![RnaBase::G, RnaBase::A, RnaBase::U, RnaBase::C]);

        assert_eq!(dna.transcribe(), expected_rna);
    }
}

#[cfg(test)]
mod serialization_tests {
    use super::*;

    #[test]
    fn test_dna_sequence_serialization_round_trip() {
        let original_raw = "GATGGAACTTGACTACGTAAATT";

        // 1. Deserialize (String -> Struct)
        let sequence: DnaSequence = original_raw
            .parse()
            .expect("Valid DNA string should deserialize seamlessly");

        // Verify internal structural representation
        assert_eq!(sequence.0[0], DnaBase::G);
        assert_eq!(sequence.0[2], DnaBase::T);

        // 2. Serialize (Struct -> String)
        let serialized_output = sequence.to_string();

        // Assert perfect round-trip parity
        assert_eq!(serialized_output, original_raw);
    }

    #[test]
    fn test_deserialization_failure_on_invalid_text() {
        // 'X' is invalid in DNA
        let corrupt_dna = "GATGXAACTT";
        let result: Result<DnaSequence, _> = corrupt_dna.parse();
        assert_eq!(result, Err(InvalidDnaSymbolError('X')));
    }
}

#[cfg(test)]
mod algorithm_tests {
    use super::*;

    #[test]
    fn test_count_empty_dna_sequence() {
        let empty_seq = DnaSequence(vec![]);
        let expected = DnaBaseCounts {
            a: 0,
            c: 0,
            g: 0,
            t: 0,
        };
        assert_eq!(empty_seq.count_bases(), expected);
    }

    #[test]
    fn test_count_dna_bases() {
        let sample_input = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";
        let expected_output = "20 12 17 21";

        let sequence: DnaSequence = sample_input
            .parse()
            .expect("Sample dataset must contain only valid DNA symbols");

        let counts = sequence.count_bases();

        assert_eq!(counts.a, 20);
        assert_eq!(counts.c, 12);
        assert_eq!(counts.g, 17);
        assert_eq!(counts.t, 21);

        assert_eq!(counts.to_string(), expected_output);
    }

    #[test]
    fn test_reverse_complement_empty() {
        let empty_seq = DnaSequence(vec![]);
        assert_eq!(empty_seq.reverse_complement(), DnaSequence(vec![]));
    }

    #[test]
    fn test_reverse_complement_single_bases() {
        let a = DnaSequence(vec![DnaBase::A]);
        let t = DnaSequence(vec![DnaBase::T]);
        let c = DnaSequence(vec![DnaBase::C]);
        let g = DnaSequence(vec![DnaBase::G]);

        assert_eq!(a.reverse_complement(), t);
        assert_eq!(t.reverse_complement(), a);
        assert_eq!(c.reverse_complement(), g);
        assert_eq!(g.reverse_complement(), c);
    }

    #[test]
    fn test_reverse_complement_sequence() {
        // "AAAACCCGGT" reversed is "TGGCCCAAAA", complemented is "ACCGGGTTTT"
        let input = "AAAACCCGGT";
        let expected_output = "ACCGGGTTTT";

        let sequence: DnaSequence = input.parse().expect("Valid DNA string");

        let result = sequence.reverse_complement();
        assert_eq!(result.to_string(), expected_output);
    }

    #[test]
    fn test_gc_content_empty() {
        let empty_seq = DnaSequence(vec![]);
        assert_eq!(empty_seq.gc_content(), 0.0);
    }

    #[test]
    fn test_gc_content_standard_vector() {
        // User example: "AGCTATAG" has 3 GC pairs out of 8 total bases (37.5%)
        let sequence: DnaSequence = "AGCTATAG".parse().expect("Valid test string");

        let expected = 37.5;
        assert!((sequence.gc_content() - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn test_gc_content_pure_gc() {
        let sequence: DnaSequence = "GGCCGGCCGG".parse().expect("Valid test string");

        assert_eq!(sequence.gc_content(), 100.0);
    }

    #[test]
    fn test_hamming_distance_identical() {
        let seq1: DnaSequence = "GAGCCTACTAACGGGAT".parse().expect("Valid DNA");
        let seq2: DnaSequence = "GAGCCTACTAACGGGAT".parse().expect("Valid DNA");

        assert_eq!(seq1.hamming_distance(&seq2), Ok(0));
    }

    #[test]
    fn test_hamming_distance_different() {
        let seq1: DnaSequence = "GAGCCTACTAACGGGAT".parse().expect("Valid DNA");
        let seq2: DnaSequence = "CATCGTAATGACGGCCT".parse().expect("Valid DNA");

        // Differences at positions: 0, 2, 3, 5, 8, 10, 15 (7 total mismatches)
        assert_eq!(seq1.hamming_distance(&seq2), Ok(7));
    }

    #[test]
    fn test_hamming_distance_empty() {
        let seq1 = DnaSequence(vec![]);
        let seq2 = DnaSequence(vec![]);

        assert_eq!(seq1.hamming_distance(&seq2), Ok(0));
    }

    #[test]
    fn test_hamming_distance_length_mismatch() {
        let seq1: DnaSequence = "GATTACA".parse().expect("Valid DNA");
        let seq2: DnaSequence = "GATTACAG".parse().expect("Valid DNA");

        let expected_err = Err(LengthMismatchError {
            sequence_len: 7,
            other_len: 8,
        });

        assert_eq!(seq1.hamming_distance(&seq2), expected_err);
    }

    #[test]
    fn test_find_motif() {
        let s: DnaSequence = "GATATATGCATATACTT".parse().unwrap();
        let t: DnaSequence = "ATAT".parse().unwrap();

        let expected = vec![2, 4, 10];
        assert_eq!(s.find_motif(&t), expected);
    }

    #[test]
    fn test_consensus() {
        let input = vec![
            "ATCCAGCT".parse::<DnaSequence>().unwrap(),
            "GGGCAACT".parse::<DnaSequence>().unwrap(),
            "ATGGATCT".parse::<DnaSequence>().unwrap(),
            "AAGCAACC".parse::<DnaSequence>().unwrap(),
            "TTGGAACT".parse::<DnaSequence>().unwrap(),
            "ATGCCATT".parse::<DnaSequence>().unwrap(),
            "ATGGCACT".parse::<DnaSequence>().unwrap(),
        ];

        let result = DnaSequence::consensus(&input).unwrap();
        assert_eq!(result.consensus.to_string(), "ATGCAACT");
    }

    #[test]
    fn test_overlap_graph() {
        let sequences = vec![
            "AAATAAA".parse::<DnaSequence>().unwrap(), // Index 0
            "AAATTTT".parse::<DnaSequence>().unwrap(), // Index 1
            "TTTTCCC".parse::<DnaSequence>().unwrap(), // Index 2
            "AAATCCC".parse::<DnaSequence>().unwrap(), // Index 3
            "GGGTGGG".parse::<DnaSequence>().unwrap(), // Index 4
        ];

        let edges = DnaSequence::overlap_graph(&sequences, 3);

        // For k=3, the expected directed overlaps are:
        // AAATAAA (0) --[AAA]--> AAATTTT (1)
        // AAATAAA (0) --[AAA]--> AAATCCC (3)
        // AAATTTT (1) --[TTT]--> TTTTCCC (2)
        let expected = vec![(0, 1), (0, 3), (1, 2)];

        assert_eq!(edges, expected);
    }
}
