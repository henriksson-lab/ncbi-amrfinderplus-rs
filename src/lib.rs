//! Rust implementation of NCBI AMRFinderPlus for identifying antimicrobial
//! resistance genes and point mutations in bacterial protein and nucleotide sequences.

pub mod alignment;
pub mod amr_reportcli;
pub mod api;
pub mod columns;
pub mod dna_mutation;
pub mod fasta_utils;
pub mod gff;
pub mod graph;
pub mod pipeline;
pub mod report;
pub mod search;
pub mod seq;
pub mod tsv;
pub mod update;

pub use api::{AmrFinder, AmrFinderBuilder, AmrFinderRun};
