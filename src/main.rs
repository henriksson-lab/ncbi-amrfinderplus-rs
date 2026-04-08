use std::io;
use std::path::PathBuf;
use std::process;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "amrfinder", version = "0.1.0")]
#[command(about = "Antimicrobial resistance gene detection")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
enum Commands {
    /// Run the AMRFinder pipeline
    Run {
        /// Input protein FASTA file
        #[arg(short = 'p', long)]
        protein: Option<PathBuf>,

        /// Input nucleotide FASTA file
        #[arg(short = 'n', long)]
        nucleotide: Option<PathBuf>,

        /// GFF file for protein locations
        #[arg(short = 'g', long)]
        gff: Option<PathBuf>,

        /// AMRFinder database directory
        #[arg(short = 'd', long)]
        database: Option<PathBuf>,

        /// Taxonomy group for mutations
        #[arg(short = 'O', long, default_value = "")]
        organism: String,

        /// Minimum identity (0..1), -1 for curated threshold
        #[arg(short = 'i', long = "ident_min", default_value = "-1")]
        ident_min: f64,

        /// Minimum coverage of reference protein (0..1)
        #[arg(short = 'c', long = "coverage_min", default_value = "0.5")]
        coverage_min: f64,

        /// Number of threads
        #[arg(long, default_value = "4")]
        threads: usize,

        /// Include plus genes
        #[arg(long)]
        plus: bool,

        /// Print hierarchy node
        #[arg(long = "print_node")]
        print_node: bool,

        /// File to report all mutations
        #[arg(long = "mutation_all")]
        mutation_all: Option<PathBuf>,

        /// GFF annotation format
        #[arg(short = 'a', long = "annotation_format", default_value = "genbank")]
        annotation_format: String,

        /// NCBI genetic code
        #[arg(short = 't', long = "translation_table", default_value = "11")]
        translation_table: u32,

        /// Name to add as first column
        #[arg(long, default_value = "")]
        name: String,

        /// BLAST binary directory
        #[arg(long = "blast_bin", default_value = "")]
        blast_bin: String,

        /// HMMer binary directory
        #[arg(long = "hmmer_bin", default_value = "")]
        hmmer_bin: String,

        /// Write output to file instead of stdout
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,

        /// Report proteins common to a taxonomy group
        #[arg(long = "report_common")]
        report_common: bool,

        /// Report all equally-scoring BLAST and HMM matches
        #[arg(long = "report_all_equal")]
        report_all_equal: bool,

        /// Print database version and exit
        #[arg(short = 'V', long = "database_version")]
        database_version: bool,

        /// List available organisms for mutation detection
        #[arg(short = 'l', long = "list_organisms")]
        list_organisms: bool,
    },

    /// Update the AMRFinder database
    Update {
        /// Force update even if already up-to-date
        #[arg(long)]
        force: bool,

        /// Database directory
        #[arg(short = 'd', long)]
        database: Option<PathBuf>,
    },

    /// Check the correctness of a FASTA file
    CheckFasta {
        /// Input FASTA file
        input: PathBuf,

        #[arg(long)]
        aa: bool,

        #[arg(long)]
        hyphen: bool,

        #[arg(long)]
        ambig: bool,

        #[arg(long = "ambig_max", default_value = "0")]
        ambig_max: usize,

        #[arg(long = "stop_codon")]
        stop_codon: bool,

        #[arg(long)]
        len: Option<PathBuf>,

        #[arg(long)]
        out: Option<PathBuf>,
    },

    /// Extract sequences from a FASTA file
    Extract {
        fasta: PathBuf,
        target: PathBuf,
        #[arg(long)]
        aa: bool,
    },

    /// Split FASTA file into parts
    SplitFasta {
        input: PathBuf,
        parts_max: usize,
        dir: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Run {
            protein,
            nucleotide,
            gff,
            database,
            organism,
            ident_min,
            coverage_min,
            threads,
            plus,
            print_node,
            mutation_all,
            annotation_format,
            translation_table,
            name,
            blast_bin,
            hmmer_bin,
            output,
            report_common: _,
            report_all_equal: _,
            database_version,
            list_organisms,
        } => {
            let db = database.unwrap_or_else(|| {
                std::env::var("AMRFINDER_DB")
                    .map(PathBuf::from)
                    .unwrap_or_default()
            });

            // Handle --database_version
            if database_version {
                let version_file = db.join("version.txt");
                if version_file.exists() {
                    let version = std::fs::read_to_string(&version_file).unwrap_or_default();
                    println!("{}", version.trim());
                } else {
                    eprintln!("Database version file not found");
                }
                return;
            }

            // Handle --list_organisms
            if list_organisms {
                let taxgroup_file = db.join("taxgroup.tsv");
                if taxgroup_file.exists() {
                    if let Ok(content) = std::fs::read_to_string(&taxgroup_file) {
                        for line in content.lines() {
                            if line.starts_with('#') {
                                continue;
                            }
                            if let Some(org) = line.split('\t').next() {
                                println!("{}", org);
                            }
                        }
                    }
                } else {
                    eprintln!("taxgroup.tsv not found in database directory");
                }
                return;
            }

            let config = amrfinder::pipeline::PipelineConfig {
                protein,
                nucleotide,
                gff,
                database: db,
                organism,
                ident_min,
                coverage_min,
                threads,
                plus,
                print_node,
                mutation_all,
                annotation_format,
                translation_table,
                name,
                blast_bin,
                hmmer_bin,
                output,
            };

            match amrfinder::pipeline::run_pipeline(&config) {
                Ok(output_text) => {
                    if config.output.is_none() {
                        print!("{}", output_text);
                    }
                    Ok(())
                }
                Err(e) => Err(e),
            }
        }
        Commands::Update { force: _, database } => {
            let db = database.unwrap_or_else(|| {
                std::env::var("AMRFINDER_DB")
                    .map(PathBuf::from)
                    .unwrap_or_else(|_| PathBuf::from("/usr/local/share/amrfinder/data"))
            });
            eprintln!("Database update not yet implemented. Database dir: {}", db.display());
            Ok(())
        }
        Commands::CheckFasta {
            input, aa, hyphen, ambig, ambig_max, stop_codon, len, out,
        } => {
            let opts = amrfinder::fasta_utils::FastaCheckOpts {
                fasta_path: &input,
                aa, hyphen, ambig, ambig_max, stop_codon,
                len_path: len.as_deref(),
                out_path: out.as_deref(),
            };
            match amrfinder::fasta_utils::fasta_check(&opts) {
                Ok((num_seqs, max_len, total_len)) => {
                    amrfinder::fasta_utils::fasta_check_print(num_seqs, max_len, total_len);
                    Ok(())
                }
                Err(e) => Err(e),
            }
        }
        Commands::Extract { fasta, target, aa } => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            amrfinder::fasta_utils::fasta_extract(&fasta, &target, aa, &mut out)
        }
        Commands::SplitFasta { input, parts_max, dir } => {
            amrfinder::fasta_utils::fasta2parts(&input, parts_max, &dir)
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
