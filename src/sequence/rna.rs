use crate::sequence::protein::AminoAcid;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub struct InvalidRnaSymbolError(pub char);

#[derive(Debug, PartialEq, Eq)]
pub enum TranslationError {
    IncompleteCodon,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RnaBase {
    A,
    C,
    G,
    U,
}

impl fmt::Display for RnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            RnaBase::A => 'A',
            RnaBase::C => 'C',
            RnaBase::G => 'G',
            RnaBase::U => 'U',
        };
        write!(f, "{}", c)
    }
}

impl TryFrom<char> for RnaBase {
    type Error = InvalidRnaSymbolError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(RnaBase::A),
            'C' => Ok(RnaBase::C),
            'G' => Ok(RnaBase::G),
            'U' => Ok(RnaBase::U),
            _ => Err(InvalidRnaSymbolError(c)),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RnaSequence(pub Vec<RnaBase>);

impl FromStr for RnaSequence {
    type Err = InvalidRnaSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bases = s
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' => Ok(RnaBase::A),
                'C' => Ok(RnaBase::C),
                'G' => Ok(RnaBase::G),
                'U' => Ok(RnaBase::U),
                _ => Err(InvalidRnaSymbolError(c)),
            })
            .collect::<Result<Vec<RnaBase>, _>>()?;
        Ok(RnaSequence(bases))
    }
}

impl fmt::Display for RnaSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in &self.0 {
            let c = match base {
                RnaBase::A => 'A',
                RnaBase::C => 'C',
                RnaBase::G => 'G',
                RnaBase::U => 'U',
            };
            write!(f, "{}", c)?;
        }
        Ok(())
    }
}

impl RnaSequence {
    pub fn translate(&self) -> Result<Vec<AminoAcid>, TranslationError> {
        if self.0.len() % 3 != 0 {
            return Err(TranslationError::IncompleteCodon);
        }

        self.0
            .chunks_exact(3)
            .map(|codon| codon_to_amino_acid([codon[0], codon[1], codon[2]]))
            .collect()
    }
}

fn codon_to_amino_acid(codon: [RnaBase; 3]) -> Result<AminoAcid, TranslationError> {
    match codon {
        [RnaBase::U, RnaBase::U, RnaBase::U] | [RnaBase::U, RnaBase::U, RnaBase::C] => {
            Ok(AminoAcid::F)
        }
        [RnaBase::U, RnaBase::U, RnaBase::A]
        | [RnaBase::U, RnaBase::U, RnaBase::G]
        | [RnaBase::C, RnaBase::U, RnaBase::U]
        | [RnaBase::C, RnaBase::U, RnaBase::C]
        | [RnaBase::C, RnaBase::U, RnaBase::A]
        | [RnaBase::C, RnaBase::U, RnaBase::G] => Ok(AminoAcid::L),
        [RnaBase::A, RnaBase::U, RnaBase::U]
        | [RnaBase::A, RnaBase::U, RnaBase::C]
        | [RnaBase::A, RnaBase::U, RnaBase::A] => Ok(AminoAcid::I),
        [RnaBase::A, RnaBase::U, RnaBase::G] => Ok(AminoAcid::M),
        [RnaBase::G, RnaBase::U, RnaBase::U]
        | [RnaBase::G, RnaBase::U, RnaBase::C]
        | [RnaBase::G, RnaBase::U, RnaBase::A]
        | [RnaBase::G, RnaBase::U, RnaBase::G] => Ok(AminoAcid::V),
        [RnaBase::U, RnaBase::C, RnaBase::U]
        | [RnaBase::U, RnaBase::C, RnaBase::C]
        | [RnaBase::U, RnaBase::C, RnaBase::A]
        | [RnaBase::U, RnaBase::C, RnaBase::G]
        | [RnaBase::A, RnaBase::G, RnaBase::U]
        | [RnaBase::A, RnaBase::G, RnaBase::C] => Ok(AminoAcid::S),
        [RnaBase::C, RnaBase::C, RnaBase::U]
        | [RnaBase::C, RnaBase::C, RnaBase::C]
        | [RnaBase::C, RnaBase::C, RnaBase::A]
        | [RnaBase::C, RnaBase::C, RnaBase::G] => Ok(AminoAcid::P),
        [RnaBase::A, RnaBase::C, RnaBase::U]
        | [RnaBase::A, RnaBase::C, RnaBase::C]
        | [RnaBase::A, RnaBase::C, RnaBase::A]
        | [RnaBase::A, RnaBase::C, RnaBase::G] => Ok(AminoAcid::T),
        [RnaBase::G, RnaBase::C, RnaBase::U]
        | [RnaBase::G, RnaBase::C, RnaBase::C]
        | [RnaBase::G, RnaBase::C, RnaBase::A]
        | [RnaBase::G, RnaBase::C, RnaBase::G] => Ok(AminoAcid::A),
        [RnaBase::U, RnaBase::A, RnaBase::U] | [RnaBase::U, RnaBase::A, RnaBase::C] => {
            Ok(AminoAcid::Y)
        }
        [RnaBase::U, RnaBase::A, RnaBase::A]
        | [RnaBase::U, RnaBase::A, RnaBase::G]
        | [RnaBase::U, RnaBase::G, RnaBase::A] => Ok(AminoAcid::Stop),
        [RnaBase::C, RnaBase::A, RnaBase::U] | [RnaBase::C, RnaBase::A, RnaBase::C] => {
            Ok(AminoAcid::H)
        }
        [RnaBase::C, RnaBase::A, RnaBase::A] | [RnaBase::C, RnaBase::A, RnaBase::G] => {
            Ok(AminoAcid::Q)
        }
        [RnaBase::A, RnaBase::A, RnaBase::U] | [RnaBase::A, RnaBase::A, RnaBase::C] => {
            Ok(AminoAcid::N)
        }
        [RnaBase::A, RnaBase::A, RnaBase::A] | [RnaBase::A, RnaBase::A, RnaBase::G] => {
            Ok(AminoAcid::K)
        }
        [RnaBase::G, RnaBase::A, RnaBase::U] | [RnaBase::G, RnaBase::A, RnaBase::C] => {
            Ok(AminoAcid::D)
        }
        [RnaBase::G, RnaBase::A, RnaBase::A] | [RnaBase::G, RnaBase::A, RnaBase::G] => {
            Ok(AminoAcid::E)
        }
        [RnaBase::U, RnaBase::G, RnaBase::U] | [RnaBase::U, RnaBase::G, RnaBase::C] => {
            Ok(AminoAcid::C)
        }
        [RnaBase::U, RnaBase::G, RnaBase::G] => Ok(AminoAcid::W),
        [RnaBase::C, RnaBase::G, RnaBase::U]
        | [RnaBase::C, RnaBase::G, RnaBase::C]
        | [RnaBase::C, RnaBase::G, RnaBase::A]
        | [RnaBase::C, RnaBase::G, RnaBase::G]
        | [RnaBase::A, RnaBase::G, RnaBase::A]
        | [RnaBase::A, RnaBase::G, RnaBase::G] => Ok(AminoAcid::R),
        [RnaBase::G, RnaBase::G, RnaBase::U]
        | [RnaBase::G, RnaBase::G, RnaBase::C]
        | [RnaBase::G, RnaBase::G, RnaBase::A]
        | [RnaBase::G, RnaBase::G, RnaBase::G] => Ok(AminoAcid::G),
    }
}

#[cfg(test)]
mod serialization_tests {
    use super::*;

    #[test]
    fn test_rna_sequence_serialization_round_trip() {
        let original_raw = "GAUGGAACUUGACUACGUAAAUU";

        // 1. Deserialize
        let sequence: RnaSequence = original_raw
            .parse()
            .expect("Valid RNA string should deserialize seamlessly");

        // 2. Serialize
        let serialized_output = sequence.to_string();

        assert_eq!(serialized_output, original_raw);
    }

    #[test]
    fn test_deserialization_failure_on_invalid_text() {
        // 'T' is invalid in RNA
        let corrupt_rna = "GAUGUAACTT";
        let rna_result: Result<RnaSequence, _> = corrupt_rna.parse();
        assert_eq!(rna_result, Err(InvalidRnaSymbolError('T')));
    }
}

#[test]
fn test_rna_translation() {
    let rna: RnaSequence = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        .parse()
        .unwrap();
    let protein = rna.translate().unwrap();
    // Verify first few amino acids
    assert_eq!(protein[0], AminoAcid::M); // AUG -> Methionine
    assert_eq!(protein[1], AminoAcid::A); // GCC -> Alanine
}
