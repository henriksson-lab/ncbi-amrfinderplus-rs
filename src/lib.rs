//! Rust implementation of NCBI AMRFinderPlus for identifying antimicrobial
//! resistance genes and point mutations in bacterial protein and nucleotide sequences.

pub mod alignment;
pub mod amr_reportcli;
pub mod dna_mutation;
pub mod search;
pub mod columns;
pub mod fasta_utils;
pub mod gff;
pub mod graph;
pub mod pipeline;
pub mod report;
pub mod seq;
pub mod tsv;
pub mod update;
