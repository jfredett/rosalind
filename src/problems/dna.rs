extern crate rosalind;

use rosalind::dna::DNA;

use std::path::Path;
use std::env;

fn main() {
    let args : Vec<_> = env::args().collect();

    let file = Path::new(&args[1]);
    let dna = DNA::from_file(file).ok().unwrap();

    for (key,val) in dna.nucleotide_stats() {
        println!("{:?} = {}", key, val);
    }
}
