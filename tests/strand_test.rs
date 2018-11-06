extern crate rosalind;

use rosalind::strand::*;
use std::path::Path;

#[test]
fn test_from_file() {
    let scaffold = Path::new("tests/scaffolds/rosalind_dna_example");
    let strand = Strand::from_file(scaffold).ok().unwrap();
    assert!(strand != Strand::empty());
    let expected = Strand::from_str("GTCAATAACTATCAGCGAAAGGATACGTTATCGCTCACACAACACGCGCAGGCCTATCCCTTCGTCTGGACGACCTATTTAGCAAGTCCGAACCACAGCCTATATGACAATTCCCTATAGCATCTTAGTGACAGCCAAAATCAGTGATACCGAAGGCACTATATACGCGGCATGGTACTACAAACACCTCTGGTATGGGCTTGCTTTGGTATCGTTTCTTAATGGTACTGTGTCCGCCTTAGGTACTGCTGAGACAGAGCATGTATTTGGTAGGGGAGAGAGCTGATTTCTGAACCTACCCCCGTGTGAGCTCTAGATTTATGGGGGTTGTTCGTCGGCATACATAGGGATTTGGTAACTCTGGTTCCCAAAGGGCAAGGCCTCCCGGATTGCTTGATGAGAGTATAATTGTTGTGGCCGCCTTGTGGAATGCTAGGAACGTCACTTTACCCCCCCACAATAGCTATCCCATGCAGAGTCTAATTGCGAAGGAGGCATGAAGCGGGACAAGAAATTCGATCTTAAAACCGCCGCCACACAAAGTCCTTAAGGGAGGTTGACGCCAGTAAAGTCTGTACCCATCACCCCTATAGACTAAGCCCTGTATTCCAAAACCAGTTGACGTAGCTCATGTTGATGAGCGCATGTTTACATGACCCGGACACGTGAAAGCGTCCGCCACACCAGAACGCTGGCACTTGCTAACTGCAGCCGTCCTCGGGAGCGCGACCCGTTGCCCGGGCAGGGTCGAATGTTTGTAAAATTTCGTTAGTCCCGTGGTCTGACGCAATCTCATCTTTACCTGGCCTGAGTACGATGATTGTACAGATTGCTTAACAACGGATACACAAA").ok().unwrap();
    assert_eq!(strand, expected);
}
