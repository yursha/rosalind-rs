use std::fmt;
use std::str::FromStr;

use crate::rna::{RnaBase, RnaSequence};

#[derive(Debug, PartialEq, Eq)]
pub struct InvalidDnaSymbolError(pub char);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
}
