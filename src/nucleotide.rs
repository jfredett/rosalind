use crate::strand::mode::StrandMode;

#[derive(PartialEq, Eq, Debug, Hash, Clone, Copy)]
/// Represents a single base Nucleotide as a two-bit string.
///
/// Values are mapped such that logical NOT corresponds to nucleotide complement.
pub enum Nucleotide {
    A = 0, // !T
    C = 1, // !G
    G = 2, // !C
    T = 3  // !A
}

impl Nucleotide {
    pub fn to_char(&self, mode : StrandMode) -> char {
        match self {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => { if mode == StrandMode::DNA { 'T' } else { 'U' } },
        }
    }
}

impl From<u64> for Nucleotide {
    fn from(v: u64) -> Self {
        match v {
            0 => Nucleotide::A,
            1 => Nucleotide::C,
            2 => Nucleotide::G,
            3 => Nucleotide::T,
            _ => panic!("Cannot convert value to Nucleotide, value `{:?}` too large", v)
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;


    mod to_char_with_rna {
        use super::*;

        #[test]
        fn to_char_a() {
            assert_eq!(Nucleotide::A.to_char(StrandMode::RNA), 'A');
        }

        #[test]
        fn to_char_c() {
            assert_eq!(Nucleotide::C.to_char(StrandMode::RNA), 'C');
        }

        #[test]
        fn to_char_g() {
            assert_eq!(Nucleotide::G.to_char(StrandMode::RNA), 'G');
        }

        #[test]
        fn to_char_t() {
            assert_eq!(Nucleotide::T.to_char(StrandMode::RNA), 'U');
        }
    }
    mod to_char_with_dna {
        use super::*;

        #[test]
        fn to_char_a() {
            assert_eq!(Nucleotide::A.to_char(StrandMode::DNA), 'A');
        }

        #[test]
        fn to_char_c() {
            assert_eq!(Nucleotide::C.to_char(StrandMode::DNA), 'C');
        }

        #[test]
        fn to_char_g() {
            assert_eq!(Nucleotide::G.to_char(StrandMode::DNA), 'G');
        }

        #[test]
        fn to_char_t() {
            assert_eq!(Nucleotide::T.to_char(StrandMode::DNA), 'T');
        }
    }

    mod from {
        use super::*;

        #[test]
        fn from_u64_a() {
            assert_eq!(Nucleotide::from(0u64), Nucleotide::A);
        }

        #[test]
        fn from_u64_t() {
            assert_eq!(Nucleotide::from(3u64), Nucleotide::T);
        }

        #[test]
        fn from_u64_c() {
            assert_eq!(Nucleotide::from(1u64), Nucleotide::C);
        }

        #[test]
        fn from_u64_g() {
            assert_eq!(Nucleotide::from(2u64), Nucleotide::G);
        }

    }
}
