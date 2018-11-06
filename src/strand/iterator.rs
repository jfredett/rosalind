use crate::strand::*;

pub struct StrandIterator {
    strand: Strand,
    idx: usize
}

impl IntoIterator for Strand {
    type Item = Nucleotide;
    type IntoIter = StrandIterator;

    fn into_iter(self) -> Self::IntoIter {
        StrandIterator { strand: self, idx: 0 }
    }

}

impl Iterator for StrandIterator {
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Nucleotide> {
        // we want bitpairs
        let mask = 0b11;
        // XXX: 32 == number of nucleotides per block, probably refactor this to calculate
        // from the type?
        let block_idx = self.idx / 32;
        let shift = 2 * (self.idx % 32);

        let element = self.strand.content[block_idx];

        // We've finished all the necessary calculations, so we can increment the counter.
        self.idx += 1;

        if block_idx == self.strand.terminus_idx && shift == self.strand.get_shift() {
            // we're on the last block and past the end of strand.
            None
        } else {
            // Extract the bitpair using our mask, convert it to a nucleotide, ship it.
            //
            // Example w/ 8-bit block size (real thing is 64 bit)
            // element = 0b10_01_11_10, shift = 2.
            // (0b10_01_10_11 & (0b11 << 2)) >> 2            ==>
            // (0b10_01_10_11 & (0b00_00_00_11_00)) >> 2     ==>
            // (0b00_00_10_00) >> 2                          ==>
            // (0b00_00_00_10)
            Some(Nucleotide::from((element & (mask << shift)) >> shift))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iteration_long() {
        let strand = Strand::from_str("GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA").ok().unwrap();
        let expected_sequence = vec![
            Nucleotide::G,
            Nucleotide::A,
            Nucleotide::T,
            Nucleotide::T,
            Nucleotide::A,
            Nucleotide::C,
            Nucleotide::A
        ];
        let mut idx = 0;
        for nucleotide in strand {
            assert_eq!(expected_sequence[idx], nucleotide);
            idx += 1;
            idx %= 7; // repeated sequence makes it easy to test, if not the ideal model of real data
        }
    }

    #[test]
    fn iteration_short() {
        let strand = Strand::from_str("GATTACA").ok().unwrap();
        let expected_sequence = vec![
            Nucleotide::G,
            Nucleotide::A,
            Nucleotide::T,
            Nucleotide::T,
            Nucleotide::A,
            Nucleotide::C,
            Nucleotide::A
        ];
        let mut idx = 0;
        for nucleotide in strand {
            assert_eq!(nucleotide, expected_sequence[idx]);
            idx += 1;
        }
    }
}

#[cfg(test)]
#[allow(unused_must_use)]
mod benches {
    use super::*;

    use test::{Bencher, black_box};
    mod iter {
        use super::*;

        #[bench]
        fn iter_short(b: &mut Bencher) {
            b.iter(|| {
                let strand = black_box(Strand::from_str("GATTACA"));
                for nucleotide in strand {
                    black_box(nucleotide);
                }
            });
        }

        #[bench]
        fn iter_long(b: &mut Bencher) {
            b.iter(|| {
                let strand = black_box(Strand::from_str("GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA"));
                for nucleotide in strand {
                    black_box(nucleotide);
                }
            });
        }
    }
}
