//! Rust implementation of NCBI AMRFinderPlus for identifying antimicrobial
//! resistance genes and point mutations in bacterial protein and nucleotide sequences.

pub mod alignment;
pub mod amr_report;
pub mod amrfinder;
pub mod amrfinder_index;
pub mod amrfinder_update;
pub mod columns;
pub mod common;
pub mod curl_easy;
pub mod disruption2genesymbol;
pub mod dna_mutation;
pub mod fasta2parts;
pub mod fasta_check;
pub mod fasta_extract;
pub mod gff;
pub mod gff_check;
pub mod graph;
pub mod mutate;
pub mod seq;
pub mod tsv;
