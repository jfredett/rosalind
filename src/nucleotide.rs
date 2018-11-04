#[derive(PartialEq, Eq, Debug)]
/// Represents a single base Nucleotide as a two-bit string.
///
/// Values are mapped such that logical NOT corresponds to nucleotide complement.
pub enum Nucleotide {
    A = 0, // !T
    C = 1, // !G
    G = 2, // !C
    T = 3  // !A
}
