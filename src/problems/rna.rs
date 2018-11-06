extern crate rosalind;

use rosalind::strand::Strand;

use std::path::Path;
use std::env;

fn main() {
    let args : Vec<_> = env::args().collect();

    let file = Path::new(&args[1]);
    let mut strand = Strand::from_file(file).ok().unwrap();

    strand.convert();

    println!("{}", strand.full_strand());
}
