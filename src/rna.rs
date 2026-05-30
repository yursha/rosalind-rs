use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub struct InvalidRnaSymbolError(pub char);

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
