use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AminoAcid {
    A,
    R,
    N,
    D,
    C,
    Q,
    E,
    G,
    H,
    I,
    L,
    K,
    M,
    F,
    P,
    S,
    T,
    W,
    Y,
    V,
    Stop,
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let symbol = match self {
            AminoAcid::A => 'A',
            AminoAcid::R => 'R',
            AminoAcid::N => 'N',
            AminoAcid::D => 'D',
            AminoAcid::C => 'C',
            AminoAcid::Q => 'Q',
            AminoAcid::E => 'E',
            AminoAcid::G => 'G',
            AminoAcid::H => 'H',
            AminoAcid::I => 'I',
            AminoAcid::L => 'L',
            AminoAcid::K => 'K',
            AminoAcid::M => 'M',
            AminoAcid::F => 'F',
            AminoAcid::P => 'P',
            AminoAcid::S => 'S',
            AminoAcid::T => 'T',
            AminoAcid::W => 'W',
            AminoAcid::Y => 'Y',
            AminoAcid::V => 'V',
            AminoAcid::Stop => return Ok(()), // Print nothing for Stop
        };
        write!(f, "{}", symbol)
    }
}
