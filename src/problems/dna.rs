extern crate rosalind;

use rosalind::dna::DNA;
use rosalind::nucleotide::Nucleotide;

use std::path::Path;
use std::collections::HashMap;
use std::env;

fn main() {
    let args : Vec<_> = env::args().collect();

    let file = Path::new(&args[1]);
    let dna = DNA::from_file(file).ok().unwrap();

    let mut results = HashMap::new();

    results.insert(Nucleotide::G, 0);
    results.insert(Nucleotide::C, 0);
    results.insert(Nucleotide::T, 0);
    results.insert(Nucleotide::A, 0);

    for nucleotide in dna {
        *results.get_mut(&nucleotide).unwrap() += 1;
    }

    for (key,val) in results {
        println!("{:?} = {}", key, val);
    }
}
