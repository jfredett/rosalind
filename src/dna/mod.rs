pub use crate::nucleotide::Nucleotide;

pub mod error;
pub mod result;
pub mod iterator;

pub use crate::dna::error::*;
pub use crate::dna::result::*;
pub use crate::dna::iterator::*;

use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::collections::HashMap;

/// A Strand of DNA
#[derive(PartialEq, Eq, Debug)]
pub struct DNA {
    /// Content of 2-bit encoded nucleotides, 32 per entry, LEndian encoded (LSB = first pair)
    content: Vec<u64>,
    /// A mask indicating unused bits in the last entry of the content
    // NOTE: might want to refactor this to like, Vec<u32>, with each bit indicating if we know
    // that nucleotide, this would allow a sparse-encoding, we might have a partial DNA sequence
    // and want to represent it
    unused_mask: u64,
    /// Index of the end of the strand
    terminus_idx: usize,
}

impl DNA {
    /// Create an Empty DNA strand
    pub fn empty() -> DNA {
        DNA { content: vec![0], unused_mask: !0u64, terminus_idx: 0 }
    }

    /// Load a DNA strand from a str
    pub fn from_str(s: &str) -> DNAResult<DNA> {
        return DNA::from_string(s.to_string());
    }

    /// Load a DNA strand from a String
    pub fn from_string(s: String) -> DNAResult<DNA> {
        let mut dna = DNA::empty();
        for c in s.chars() {
            match c {
                'A' => dna.add(Nucleotide::A).ok(),
                'T' => dna.add(Nucleotide::T).ok(),
                'G' => dna.add(Nucleotide::G).ok(),
                'C' => dna.add(Nucleotide::C).ok(),
                '\n' => None,
                _ => { return Err(DNAError::InvalidNucleotide); }
            };
        }

        return Ok(dna);
    }

    pub fn nucleotide_stats(self) -> HashMap<Nucleotide,i64> {
        let mut results = HashMap::new();

        results.insert(Nucleotide::G, 0);
        results.insert(Nucleotide::C, 0);
        results.insert(Nucleotide::T, 0);
        results.insert(Nucleotide::A, 0);

        for nucleotide in self {
            *results.get_mut(&nucleotide).unwrap() += 1;
        }

        return results;
    }

    /// A string representation of the strand
    ///
    /// # Examples
    ///
    /// ```
    /// use rosalind::dna::*;
    /// let dna = DNA::from_str("GATTACA").ok().unwrap();
    /// assert_eq!(dna.full_strand(), "GATTACA");
    /// ```
    pub fn full_strand(self) -> String {
        let mut ret = String::new();
        for nucleotide in self {
            ret.push(nucleotide.to_char());
        }
        return ret;
    }

    /// Load a DNA strand from a filepath
    pub fn from_file<P: AsRef<Path>>(path: P) -> DNAResult<DNA> {
        let mut file = File::open(path).ok().unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).ok();
        return DNA::from_string(contents);
    }

    /// Add a single nucleotide to the end of the DNA strand
    ///
    /// # Psuedo-code:
    ///
    /// if end-of-strand is full
    ///      increment terminus_idx
    ///      allocate another 0 at terminus_idx
    ///      reset mask to !0
    /// end
    /// Find first available 2-bit position -- should be zeroes by default by above
    /// OR the nucleotide code << that position to the vector
    /// AND the used mask (!unused) to fix any messed up bits
    /// return Result
    ///
    pub fn add(&mut self, n: Nucleotide) -> DNAResult<()> {
        if self.is_full() {
            self.content.push(0);
            self.unused_mask = !0u64;
            self.terminus_idx += 1;
        }
        let vector = (n as u64) << self.get_shift();
        let terminal = self.content[self.terminus_idx];
        self.unused_mask <<= 2; // shift by two bits
        self.content[self.terminus_idx] = (terminal | vector) & !self.unused_mask;

        Ok(())
    }

    /// Checks if the terminal slot is full
    fn is_full(&self) -> bool {
        self.unused_mask == 0
    }

    /// find the offset of the next available slot
    fn get_shift(&self) -> usize {
        let mut u = 0;
        let mut mask = self.unused_mask;
        while mask & 1 != 1 {
            mask >>= 1;
            u += 1;
        }
        return u;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod display {
        use super::*;

        #[test]
        fn full_strand() {
            let dna = DNA::from_str("GATTACA").ok().unwrap();
            assert_eq!(dna.full_strand(), "GATTACA");
        }

    }

    mod stats {
        use super::*;

        #[test]
        fn nucleotide_stats() {
            let dna = DNA::from_str("GATTACA").ok().unwrap();
            let stats = dna.nucleotide_stats();
            assert_eq!(stats.get(&Nucleotide::A), Some(&3));
            assert_eq!(stats.get(&Nucleotide::T), Some(&2));
            assert_eq!(stats.get(&Nucleotide::C), Some(&1));
            assert_eq!(stats.get(&Nucleotide::G), Some(&1));
        }
    }

    mod creation {
        use super::*;

        #[test]
        fn empty() {
            let dna = DNA::empty();
            assert_eq!(dna.content, vec![0]);
            assert_eq!(dna.unused_mask, !0);
            assert_eq!(dna.terminus_idx, 0);
        }

        #[test]
        fn from_string() {
            let s = String::from("ACTG");
            let dna = DNA::from_string(s).ok().unwrap();
            //                           G T C A
            assert_eq!(dna.content[0], 0b10110100);
        }

        #[test]
        fn from_str() {
            let dna = DNA::from_str("ACTG").ok().unwrap();
            //                           G T C A
            assert_eq!(dna.content[0], 0b10110100);
        }
    }

    mod adding_nucleotides {
        use super::*;

        #[test]
        fn add_a() {
            let mut dna = DNA::empty();
            dna.add(Nucleotide::A).ok();

            assert_eq!(dna.content[0], 0b00);
        }

        #[test]
        fn add_c() {
            let mut dna = DNA::empty();
            dna.add(Nucleotide::C).ok();

            assert_eq!(dna.content[0], 0b01);
        }

        #[test]
        fn add_g() {
            let mut dna = DNA::empty();
            dna.add(Nucleotide::G).ok();

            assert_eq!(dna.content[0], 0b10);
        }

        #[test]
        fn add_t() {
            let mut dna = DNA::empty();
            dna.add(Nucleotide::T).ok();

            assert_eq!(dna.content[0], 0b11);
        }

        #[test]
        fn add_multiple() {
            let mut dna = DNA::empty();
            dna.add(Nucleotide::T).ok();
            dna.add(Nucleotide::A).ok();
            dna.add(Nucleotide::G).ok();
            dna.add(Nucleotide::C).ok();
            //                           C G A T
            assert_eq!(dna.content[0], 0b01100011);
        }

        #[test]
        fn add_multiple_with_allocation() {
            let mut dna = DNA::empty();

            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();
            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();
            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();
            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();

            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();
            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();
            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();
            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();

            dna.add(Nucleotide::T).ok(); dna.add(Nucleotide::A).ok(); dna.add(Nucleotide::G).ok(); dna.add(Nucleotide::C).ok();

            assert_eq!(dna.terminus_idx, 1);
            assert_eq!(dna.content[0], 0x63_63_63_63_63_63_63_63);
            assert_eq!(dna.content[1], 0x63);
        }
    }
}

#[cfg(test)]
#[allow(unused_must_use)]
mod benches {
    use super::*;

    use test::{Bencher, black_box};

    mod creation {
        use super::*;

        #[bench]
        fn empty(b: &mut Bencher) {
            b.iter(|| {
                black_box(DNA::empty());
            });
        }

        #[bench]
        fn from_string_short(b: &mut Bencher) {
            b.iter(|| {
                black_box(DNA::from_string("GATTACA".to_string()));
            });
        }

        #[bench]
        fn from_string_long(b: &mut Bencher) {
            b.iter(|| {
                black_box(DNA::from_string("GATTACAATATGGAGTATCAGCTGCATCGCGATTCGAGGATTCGAGAGACTTTGAACAGCCACCCACGTTCCTCAGAGAGAGCGCGTCA".to_string()));
            });
        }

        #[bench]
        fn from_str_short(b: &mut Bencher) {
            b.iter(|| {
                black_box(DNA::from_str("GATTACA"));
            });
        }

        #[bench]
        fn from_str_long(b: &mut Bencher) {
            b.iter(|| {
                black_box(DNA::from_str("GATTACAATATGGAGTATCAGCTGCATCGCGATTCGAGGATTCGAGAGACTTTGAACAGCCACCCACGTTCCTCAGAGAGAGCGCGTCA"));
            });
        }

        #[bench]
        fn from_file(b: &mut Bencher) {
            b.iter(|| {
                black_box(DNA::from_string("GATTACAATATGGAGTATCAGCTGCATCGCGATTCGAGGATTCGAGAGACTTTGAACAGCCACCCACGTTCCTCAGAGAGAGCGCGTCA".to_string()));
            });
        }
    }

    mod adding_nucleotides {
        use super::*;

        #[bench]
        fn add_single_nucleotide(b: &mut Bencher) {
            let mut dna_strand = DNA::empty();
            b.iter(|| {
                dna_strand.add(Nucleotide::T);
            });
        }

        #[bench]
        fn add_multiple_nucleotides_without_allocation(b: &mut Bencher) {
            let mut dna_strand = DNA::empty();
            b.iter(|| {
                dna_strand.add(Nucleotide::T);
                dna_strand.add(Nucleotide::A);
                dna_strand.add(Nucleotide::C);
                dna_strand.add(Nucleotide::G);
            });
        }

        #[bench]
        fn add_multiple_nucleotides_with_allocation(b: &mut Bencher) {
            let mut dna_strand = DNA::empty();
            b.iter(|| {
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);

                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);

                dna_strand.add(Nucleotide::T); dna_strand.add(Nucleotide::A); dna_strand.add(Nucleotide::C); dna_strand.add(Nucleotide::G);
            });
        }
    }

    mod stats {
        use super::*;

        #[bench]
        fn nucleotide_stats_short(b: &mut Bencher) {
            b.iter(|| {
                let dna = DNA::from_str("GATTACA").ok().unwrap();
                let stats = dna.nucleotide_stats();
                black_box(stats);
            });
        }

        #[bench]
        fn nucleotide_stats_long(b: &mut Bencher) {
            b.iter(|| {
                let dna = DNA::from_str("GATTACAGATTATTTATTATCACGTGACGGCATCGGCAGCAGCATACTCTGCGCAGGCGTAGCGCGATCTCTAGAGCGCGCAGATTC").ok().unwrap();
                let stats = dna.nucleotide_stats();
                black_box(stats);
            });
        }

        #[bench]
        fn nucleotide_stats_file(b: &mut Bencher) {
            b.iter(|| {
                let dna = DNA::from_file("tests/scaffolds/rosalind_dna-1.txt").ok().unwrap();
                let stats = dna.nucleotide_stats();
                black_box(stats);
            });
        }
    }

}