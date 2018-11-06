use crate::dna::*;

pub struct DNAIterator {
    dna: DNA,
    idx: usize
}

impl IntoIterator for DNA {
    type Item = Nucleotide;
    type IntoIter = DNAIterator;

    fn into_iter(self) -> Self::IntoIter {
        DNAIterator { dna: self, idx: 0 }
    }

}

impl Iterator for DNAIterator {
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Nucleotide> {
        // we want bitpairs
        let mask = 0b11;
        // XXX: 32 == number of nucleotides per block, probably refactor this to calculate
        // from the type?
        let block_idx = self.idx / 32;
        let shift = 2 * (self.idx % 32);

        let element = self.dna.content[block_idx];

        // We've finished all the necessary calculations, so we can increment the counter.
        self.idx += 1;

        if block_idx == self.dna.terminus_idx && shift == self.dna.get_shift() {
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
#[allow(unused_must_use)]
mod benches {
    use super::*;

    use test::{Bencher, black_box};
    mod iter {
        use super::*;

        #[bench]
        fn iter_short(b: &mut Bencher) {
            b.iter(|| {
                let dna = black_box(DNA::from_str("GATTACA"));
                for nucleotide in dna {
                    black_box(nucleotide);
                }
            });
        }

        #[bench]
        fn iter_long(b: &mut Bencher) {
            b.iter(|| {
                let dna = black_box(DNA::from_str("GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA"));
                for nucleotide in dna {
                    black_box(nucleotide);
                }
            });
        }
    }
}
