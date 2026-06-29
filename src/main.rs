use std::ffi::OsString;
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

        /// Common directory of protein, nucleotide, and GFF input files
        #[arg(long = "dir")]
        dir: Option<PathBuf>,

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

        /// Input files were created by NCBI PGAP
        #[arg(long = "pgap")]
        pgap: bool,

        /// Interpret organism as an NCBI internal GPipe organism name
        #[arg(long = "gpipe_org")]
        gpipe_org: bool,

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

        /// Original testing parameter string
        #[arg(long, default_value = "", allow_hyphen_values = true)]
        parm: String,

        /// Write output to file instead of stdout
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,

        /// Output protein FASTA file of reported proteins
        #[arg(long = "protein_output")]
        protein_output: Option<PathBuf>,

        /// Output nucleotide FASTA file of reported nucleotide sequences
        #[arg(long = "nucleotide_output")]
        nucleotide_output: Option<PathBuf>,

        /// Output nucleotide FASTA file of reported nucleotide sequences with 5' flanking sequences
        #[arg(long = "nucleotide_flank5_output")]
        nucleotide_flank5_output: Option<PathBuf>,

        /// 5' flanking sequence size for NUC_FLANK5_FASTA_OUT
        #[arg(long = "nucleotide_flank5_size", default_value = "0")]
        nucleotide_flank5_size: usize,

        /// Report proteins common to a taxonomy group
        #[arg(long = "report_common")]
        report_common: bool,

        /// Report all equally-scoring BLAST and HMM matches
        #[arg(long = "report_all_equal")]
        report_all_equal: bool,

        /// Suppress messages to STDERR
        #[arg(short = 'q', long)]
        quiet: bool,

        /// Integrity checks
        #[arg(long)]
        debug: bool,

        /// Integrity checks (quality control)
        #[arg(long)]
        qc: bool,

        /// Level of verbosity
        #[arg(long, default_value = "0")]
        verbose: u32,

        /// Error log file, appended
        #[arg(long = "log")]
        log: Option<PathBuf>,

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
        #[arg(short = 'U', long, alias = "force_update")]
        force: bool,

        /// Database directory
        #[arg(short = 'd', long)]
        database: Option<PathBuf>,

        /// BLAST binary directory for amrfinder_index
        #[arg(long = "blast_bin")]
        blast_bin: Option<PathBuf>,

        /// HMMer binary directory for amrfinder_index
        #[arg(long = "hmmer_bin")]
        hmmer_bin: Option<PathBuf>,

        /// Quiet index generation
        #[arg(short = 'q', long)]
        quiet: bool,

        /// Pass debug flag to amrfinder_index
        #[arg(long)]
        debug: bool,

        #[arg(long, hide = true)]
        qc: bool,

        #[arg(short = 'p', long, hide = true)]
        protein: Option<PathBuf>,

        #[arg(short = 'n', long, hide = true)]
        nucleotide: Option<PathBuf>,

        #[arg(long = "__top_level_update", hide = true)]
        top_level_update: bool,
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

        #[arg(long)]
        qc: bool,

        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Extract sequences from a FASTA file
    Extract {
        fasta: PathBuf,
        target: PathBuf,
        #[arg(long)]
        aa: bool,
        #[arg(long, default_value_t = 0)]
        verbose: u32,
        #[arg(long)]
        qc: bool,
        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Split FASTA file into parts
    SplitFasta {
        input: PathBuf,
        parts_max: usize,
        dir: PathBuf,
        #[arg(long)]
        qc: bool,
        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Mutate a FASTA file
    Mutate {
        input: PathBuf,
        mutations: PathBuf,
        #[arg(long)]
        aa: bool,
        #[arg(long)]
        orig: bool,
        #[arg(long)]
        qc: bool,
        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Check a GFF file
    GffCheck {
        gff: PathBuf,
        #[arg(long = "gfftype", default_value = "genbank")]
        gfftype: String,
        #[arg(long = "prot")]
        prot: Option<PathBuf>,
        #[arg(long = "dna")]
        dna: Option<PathBuf>,
        #[arg(long = "lcl")]
        lcl: bool,
        #[arg(long = "gff_prot_match")]
        gff_prot_match: Option<PathBuf>,
        #[arg(long = "gff_dna_match")]
        gff_dna_match: Option<PathBuf>,
        #[arg(long = "verbose", default_value = "0")]
        verbose: u32,
        #[arg(long)]
        qc: bool,
        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Find point mutations from nucleotide alignments
    DnaMutation {
        blastn: PathBuf,
        mutation: PathBuf,
        organism: String,
        #[arg(long = "mutation_all")]
        mutation_all: Option<PathBuf>,
        #[arg(long = "name", default_value = "")]
        name: String,
        #[arg(long = "print_node")]
        print_node: bool,
        #[arg(long)]
        qc: bool,
        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Convert disruption gene symbols to standard gene symbols
    Disruption2GeneSymbol {
        nucl: PathBuf,
        prot: PathBuf,
        tab: PathBuf,
        #[arg(long = "gencode", default_value = "11")]
        gencode: String,
        #[arg(long = "prot_id_pos", default_value = "0")]
        prot_id_pos: String,
        #[arg(long)]
        qc: bool,
        #[arg(long = "log")]
        log: Option<PathBuf>,
    },

    /// Update the AMRFinder database using the translated standalone updater
    AmrFinderUpdate {
        #[arg(short = 'd', long = "database")]
        database: Option<PathBuf>,
        #[arg(long = "blast_bin")]
        blast_bin: Option<PathBuf>,
        #[arg(long = "hmmer_bin")]
        hmmer_bin: Option<PathBuf>,
        #[arg(long = "force_update")]
        force_update: bool,
        #[arg(short = 'q', long = "quiet")]
        quiet: bool,
        #[arg(long)]
        debug: bool,
    },

    /// Index the AMRFinder database using the translated standalone indexer
    AmrFinderIndex {
        database: PathBuf,
        #[arg(long = "blast_bin")]
        blast_bin: Option<PathBuf>,
        #[arg(long = "hmmer_bin")]
        hmmer_bin: Option<PathBuf>,
        #[arg(short = 'q', long = "quiet")]
        quiet: bool,
        #[arg(long)]
        debug: bool,
        #[arg(long)]
        log: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse_from(args_with_default_run_subcommand());

    let result = match cli.command {
        Commands::Run {
            protein,
            nucleotide,
            gff,
            database,
            dir,
            organism,
            ident_min,
            coverage_min,
            threads,
            plus,
            print_node,
            pgap,
            gpipe_org,
            mutation_all,
            annotation_format,
            translation_table,
            name,
            blast_bin,
            hmmer_bin,
            parm,
            output,
            protein_output,
            nucleotide_output,
            nucleotide_flank5_output,
            nucleotide_flank5_size,
            report_common,
            report_all_equal,
            quiet,
            debug,
            qc,
            verbose,
            log,
            database_version,
            list_organisms,
        } => {
            let mut argv = vec![std::env::current_exe()
                .map(OsString::from)
                .unwrap_or_else(|_| OsString::from("amrfinder"))];
            if let Some(protein) = protein {
                argv.push(OsString::from("-protein"));
                argv.push(protein.into());
            }
            if let Some(nucleotide) = nucleotide {
                argv.push(OsString::from("-nucleotide"));
                argv.push(nucleotide.into());
            }
            if let Some(gff) = gff {
                argv.push(OsString::from("-gff"));
                argv.push(gff.into());
            }
            if let Some(database) =
                database.or_else(|| std::env::var_os("AMRFINDER_DB").map(PathBuf::from))
            {
                argv.push(OsString::from("-database"));
                argv.push(database.into());
            }
            if let Some(dir) = dir {
                argv.push(OsString::from("-dir"));
                argv.push(dir.into());
            }
            if !organism.is_empty() {
                argv.push(OsString::from("-organism"));
                argv.push(OsString::from(organism));
            }
            argv.push(OsString::from("-ident_min"));
            argv.push(OsString::from(ident_min.to_string()));
            argv.push(OsString::from("-coverage_min"));
            argv.push(OsString::from(coverage_min.to_string()));
            argv.push(OsString::from("-threads"));
            argv.push(OsString::from(threads.to_string()));
            if plus {
                argv.push(OsString::from("-plus"));
            }
            if report_common {
                argv.push(OsString::from("-report_common"));
            }
            if report_all_equal {
                argv.push(OsString::from("-report_all_equal"));
            }
            if quiet {
                argv.push(OsString::from("-q"));
            }
            if debug {
                argv.push(OsString::from("-debug"));
            }
            if qc {
                argv.push(OsString::from("-qc"));
            }
            if verbose != 0 {
                argv.push(OsString::from("-verbose"));
                argv.push(OsString::from(verbose.to_string()));
            }
            if let Some(log) = log {
                argv.push(OsString::from("-log"));
                argv.push(log.into());
            }
            if print_node {
                argv.push(OsString::from("-print_node"));
            }
            if pgap {
                argv.push(OsString::from("-pgap"));
            }
            if gpipe_org {
                argv.push(OsString::from("-gpipe_org"));
            }
            if let Some(mutation_all) = mutation_all {
                argv.push(OsString::from("-mutation_all"));
                argv.push(mutation_all.into());
            }
            argv.push(OsString::from("-annotation_format"));
            argv.push(OsString::from(annotation_format));
            argv.push(OsString::from("-translation_table"));
            argv.push(OsString::from(translation_table.to_string()));
            if !name.is_empty() {
                argv.push(OsString::from("-name"));
                argv.push(OsString::from(name));
            }
            if !blast_bin.is_empty() {
                argv.push(OsString::from("-blast_bin"));
                argv.push(OsString::from(blast_bin));
            }
            if !hmmer_bin.is_empty() {
                argv.push(OsString::from("-hmmer_bin"));
                argv.push(OsString::from(hmmer_bin));
            }
            if !parm.is_empty() {
                argv.push(OsString::from("-parm"));
                argv.push(OsString::from(parm));
            }
            if let Some(output) = output {
                argv.push(OsString::from("-output"));
                argv.push(output.into());
            }
            if let Some(protein_output) = protein_output {
                argv.push(OsString::from("-protein_output"));
                argv.push(protein_output.into());
            }
            if let Some(nucleotide_output) = nucleotide_output {
                argv.push(OsString::from("-nucleotide_output"));
                argv.push(nucleotide_output.into());
            }
            if let Some(nucleotide_flank5_output) = nucleotide_flank5_output {
                argv.push(OsString::from("-nucleotide_flank5_output"));
                argv.push(nucleotide_flank5_output.into());
            }
            if nucleotide_flank5_size != 0 {
                argv.push(OsString::from("-nucleotide_flank5_size"));
                argv.push(OsString::from(nucleotide_flank5_size.to_string()));
            }
            if database_version {
                argv.push(OsString::from("-database_version"));
            }
            if list_organisms {
                argv.push(OsString::from("-list_organisms"));
            }
            let stdout = io::stdout();
            let mut out = stdout.lock();
            amrfinder::amrfinder::main(argv, &mut out).map(|_| ())
        }
        Commands::Update {
            force,
            database,
            blast_bin,
            hmmer_bin,
            quiet,
            debug,
            qc,
            protein,
            nucleotide,
            top_level_update,
        } => {
            if top_level_update && (protein.is_some() || nucleotide.is_some()) {
                return_with_error(
                    "AMRFinder -u/--update option cannot be run with -n/--nucleotide or -p/--protein options",
                );
            }
            if top_level_update && database.is_some() {
                eprintln!(
                    "Error: AMRFinder update option (-u/--update) only operates on the default database directory. The -d/--database option is not permitted"
                );
                process::exit(1);
            }
            let database = if let Some(database) = database {
                Some(database)
            } else if top_level_update {
                if std::env::var_os("AMRFINDER_DB").is_some() {
                    eprintln!(
                        "AMRFINDER_DB is set, but AMRFinder auto-update only downloads to the default database directory"
                    );
                }
                let mut default_db = std::env::current_exe()
                    .ok()
                    .and_then(|exe| exe.parent().map(|dir| dir.to_path_buf()))
                    .unwrap_or_default();
                default_db.push("data");
                default_db.push("latest");
                Some(default_db.parent().map(PathBuf::from).unwrap_or(default_db))
            } else {
                None
            };

            let mut argv = vec![std::env::current_exe()
                .map(OsString::from)
                .unwrap_or_else(|_| OsString::from("amrfinder_update"))];
            if let Some(database) = database {
                argv.push(OsString::from("-database"));
                argv.push(database.into());
            }
            if let Some(blast_bin) = blast_bin {
                argv.push(OsString::from("-blast_bin"));
                argv.push(blast_bin.into());
            }
            if let Some(hmmer_bin) = hmmer_bin {
                argv.push(OsString::from("-hmmer_bin"));
                argv.push(hmmer_bin.into());
            }
            if force {
                argv.push(OsString::from("--force_update"));
            }
            if quiet {
                argv.push(OsString::from("-q"));
            }
            if debug {
                argv.push(OsString::from("--debug"));
            } else if qc {
                argv.push(OsString::from("--qc"));
            }
            amrfinder::amrfinder_update::main(argv).map(|_| ())
        }
        Commands::CheckFasta {
            input,
            aa,
            hyphen,
            ambig,
            ambig_max,
            stop_codon,
            len,
            out,
            qc: _,
            log,
        } => {
            match amrfinder::fasta_check::body(
                &input,
                aa,
                hyphen,
                ambig,
                ambig_max,
                stop_codon,
                len.as_deref(),
                out.as_deref(),
            ) {
                Ok((num_seqs, max_len, total_len)) => {
                    println!("{}", num_seqs);
                    println!("{}", max_len);
                    println!("{}", total_len);
                    if let Some(log) = log {
                        let _ = std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(log);
                    }
                    Ok(())
                }
                Err(e) => Err(e),
            }
        }
        Commands::Extract {
            fasta,
            target,
            aa,
            verbose,
            qc: _,
            log,
        } => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            let result = amrfinder::fasta_extract::body(&fasta, &target, aa, verbose > 0, &mut out);
            if result.is_ok() {
                if let Some(log) = log {
                    let _ = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(log);
                }
            }
            result
        }
        Commands::SplitFasta {
            input,
            parts_max,
            dir,
            qc: _,
            log,
        } => {
            let result = amrfinder::fasta2parts::body(&input, parts_max, &dir);
            if result.is_ok() {
                if let Some(log) = log {
                    let _ = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(log);
                }
            }
            result
        }
        Commands::Mutate {
            input,
            mutations,
            aa,
            orig,
            qc,
            log,
        } => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            let mut argv = vec![OsString::from("mutate"), input.into(), mutations.into()];
            if aa {
                argv.push(OsString::from("-aa"));
            }
            if orig {
                argv.push(OsString::from("-orig"));
            }
            if qc {
                argv.push(OsString::from("-qc"));
            }
            if let Some(log) = log {
                argv.push(OsString::from("-log"));
                argv.push(log.into());
            }
            amrfinder::mutate::main(argv, &mut out).map(|_| ())
        }
        Commands::GffCheck {
            gff,
            gfftype,
            prot,
            dna,
            lcl,
            gff_prot_match,
            gff_dna_match,
            verbose,
            qc,
            log,
        } => {
            let mut argv = vec![OsString::from("gff_check"), gff.into()];
            argv.push(OsString::from("-gfftype"));
            argv.push(OsString::from(gfftype));
            if let Some(prot) = prot {
                argv.push(OsString::from("-prot"));
                argv.push(prot.into());
            }
            if let Some(dna) = dna {
                argv.push(OsString::from("-dna"));
                argv.push(dna.into());
            }
            if lcl {
                argv.push(OsString::from("-lcl"));
            }
            if let Some(gff_prot_match) = gff_prot_match {
                argv.push(OsString::from("-gff_prot_match"));
                argv.push(gff_prot_match.into());
            }
            if let Some(gff_dna_match) = gff_dna_match {
                argv.push(OsString::from("-gff_dna_match"));
                argv.push(gff_dna_match.into());
            }
            if verbose != 0 {
                argv.push(OsString::from("-verbose"));
                argv.push(OsString::from(verbose.to_string()));
            }
            if qc {
                argv.push(OsString::from("-qc"));
            }
            if let Some(log) = log {
                argv.push(OsString::from("-log"));
                argv.push(log.into());
            }
            amrfinder::gff_check::main(argv).map(|_| ())
        }
        Commands::DnaMutation {
            blastn,
            mutation,
            organism,
            mutation_all,
            name,
            print_node,
            qc,
            log,
        } => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            let mut argv = vec![
                OsString::from("dna_mutation"),
                blastn.into(),
                mutation.into(),
                organism.into(),
            ];
            if let Some(mutation_all) = mutation_all {
                argv.push(OsString::from("-mutation_all"));
                argv.push(mutation_all.into());
            }
            if !name.is_empty() {
                argv.push(OsString::from("-name"));
                argv.push(OsString::from(name));
            }
            if print_node {
                argv.push(OsString::from("-print_node"));
            }
            if qc {
                argv.push(OsString::from("-qc"));
            }
            if let Some(log) = log {
                argv.push(OsString::from("-log"));
                argv.push(log.into());
            }
            amrfinder::dna_mutation::main(argv, &mut out).map(|_| ())
        }
        Commands::Disruption2GeneSymbol {
            nucl,
            prot,
            tab,
            gencode,
            prot_id_pos,
            qc,
            log,
        } => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            let mut argv = vec![
                OsString::from("disruption2genesymbol"),
                nucl.into(),
                prot.into(),
                tab.into(),
                OsString::from("-gencode"),
                OsString::from(gencode),
                OsString::from("-prot_id_pos"),
                OsString::from(prot_id_pos),
            ];
            if qc {
                argv.push(OsString::from("-qc"));
            }
            if let Some(log) = log {
                argv.push(OsString::from("-log"));
                argv.push(log.into());
            }
            amrfinder::disruption2genesymbol::main(argv, &mut out).map(|_| ())
        }
        Commands::AmrFinderUpdate {
            database,
            blast_bin,
            hmmer_bin,
            force_update,
            quiet,
            debug,
        } => {
            let mut argv = vec![OsString::from("amrfinder_update")];
            if let Some(database) = database {
                argv.push(OsString::from("-database"));
                argv.push(database.into());
            }
            if let Some(blast_bin) = blast_bin {
                argv.push(OsString::from("-blast_bin"));
                argv.push(blast_bin.into());
            }
            if let Some(hmmer_bin) = hmmer_bin {
                argv.push(OsString::from("-hmmer_bin"));
                argv.push(hmmer_bin.into());
            }
            if force_update {
                argv.push(OsString::from("-force_update"));
            }
            if quiet {
                argv.push(OsString::from("-q"));
            }
            if debug {
                argv.push(OsString::from("--debug"));
            }
            amrfinder::amrfinder_update::main(argv).map(|_| ())
        }
        Commands::AmrFinderIndex {
            database,
            blast_bin,
            hmmer_bin,
            quiet,
            debug,
            log,
        } => {
            let mut argv = vec![OsString::from("amrfinder_index"), database.into()];
            if let Some(blast_bin) = blast_bin {
                argv.push(OsString::from("-blast_bin"));
                argv.push(blast_bin.into());
            }
            if let Some(hmmer_bin) = hmmer_bin {
                argv.push(OsString::from("-hmmer_bin"));
                argv.push(hmmer_bin.into());
            }
            if quiet {
                argv.push(OsString::from("-q"));
            }
            if debug {
                argv.push(OsString::from("-debug"));
            }
            if let Some(log) = log {
                argv.push(OsString::from("-log"));
                argv.push(log.into());
            }
            amrfinder::amrfinder_index::main(argv).map(|_| ())
        }
    };

    if let Err(e) = result {
        return_with_error(&e.to_string());
    }
}

fn return_with_error(message: &str) -> ! {
    eprintln!("Error: {}", message);
    process::exit(1);
}

fn args_with_default_run_subcommand() -> Vec<OsString> {
    let mut args: Vec<OsString> = std::env::args_os().collect();
    if let Some(command) = args
        .first()
        .and_then(|arg| std::path::Path::new(arg).file_stem())
        .and_then(|arg| arg.to_str())
        .and_then(|program| match program {
            "fasta_check" => Some("check-fasta"),
            "fasta_extract" => Some("extract"),
            "fasta2parts" => Some("split-fasta"),
            "mutate" => Some("mutate"),
            "gff_check" => Some("gff-check"),
            "dna_mutation" => Some("dna-mutation"),
            "disruption2genesymbol" => Some("disruption2-gene-symbol"),
            "amrfinder_update" => Some("amr-finder-update"),
            "amrfinder_index" => Some("amr-finder-index"),
            _ => None,
        })
    {
        args.insert(1, OsString::from(command));
        normalize_upstream_single_dash_flags(&mut args);
        return args;
    }
    let Some(first_arg) = args.get(1) else {
        return args;
    };
    let first_arg = first_arg.to_str();
    let known_subcommand = first_arg.is_some_and(|arg| {
        matches!(
            arg,
            "run"
                | "update"
                | "check-fasta"
                | "extract"
                | "split-fasta"
                | "mutate"
                | "gff-check"
                | "dna-mutation"
                | "disruption2-gene-symbol"
                | "amr-finder-update"
                | "amr-finder-index"
        )
    });
    if !known_subcommand {
        let mut update_positions = Vec::new();
        let mut force_update = false;
        let mut skip_value = false;
        for (pos, arg) in args.iter().enumerate().skip(1) {
            let Some(value) = arg.to_str() else {
                skip_value = false;
                continue;
            };
            if skip_value {
                skip_value = false;
                continue;
            }
            match value {
                "-u" | "--update" | "-update" => update_positions.push(pos),
                "-U" | "--force_update" | "-force_update" => {
                    force_update = true;
                    update_positions.push(pos);
                }
                "-p"
                | "--protein"
                | "-protein"
                | "-n"
                | "--nucleotide"
                | "-nucleotide"
                | "-g"
                | "--gff"
                | "-gff"
                | "-d"
                | "--database"
                | "-database"
                | "--dir"
                | "-dir"
                | "-O"
                | "--organism"
                | "-organism"
                | "-i"
                | "--ident_min"
                | "-ident_min"
                | "-c"
                | "--coverage_min"
                | "-coverage_min"
                | "--threads"
                | "-threads"
                | "--mutation_all"
                | "-mutation_all"
                | "-a"
                | "--annotation_format"
                | "-annotation_format"
                | "-t"
                | "--translation_table"
                | "-translation_table"
                | "--name"
                | "-name"
                | "--blast_bin"
                | "-blast_bin"
                | "--hmmer_bin"
                | "-hmmer_bin"
                | "--parm"
                | "-parm"
                | "-o"
                | "--output"
                | "-output"
                | "--protein_output"
                | "-protein_output"
                | "--nucleotide_output"
                | "-nucleotide_output"
                | "--nucleotide_flank5_output"
                | "-nucleotide_flank5_output"
                | "--nucleotide_flank5_size"
                | "-nucleotide_flank5_size"
                | "--verbose"
                | "-verbose"
                | "--log"
                | "-log" => skip_value = true,
                _ => {}
            }
        }
        if !update_positions.is_empty() {
            for pos in update_positions.into_iter().rev() {
                args.remove(pos);
            }
            args.insert(1, OsString::from("update"));
            let mut next = 2;
            if force_update {
                args.insert(next, OsString::from("--force_update"));
                next += 1;
            }
            args.insert(next, OsString::from("--__top_level_update"));
            normalize_upstream_single_dash_flags(&mut args);
            return args;
        }
    }
    match first_arg {
        Some("-v") => {
            args[1] = OsString::from("--version");
        }
        Some(arg)
            if arg.starts_with('-') && arg != "-h" && arg != "--help" && arg != "--version" =>
        {
            args.insert(1, OsString::from("run"));
        }
        _ => {}
    }
    normalize_upstream_single_dash_flags(&mut args);
    args
}

fn normalize_upstream_single_dash_flags(args: &mut [OsString]) {
    let Some(command) = args.get(1).and_then(|arg| arg.to_str()).map(str::to_owned) else {
        return;
    };
    for arg in args.iter_mut().skip(2) {
        let Some(value) = arg.to_str() else {
            continue;
        };
        let replacement = match command.as_str() {
            "run" => match value {
                "-protein" => Some("--protein".to_string()),
                "-nucleotide" => Some("--nucleotide".to_string()),
                "-gff" => Some("--gff".to_string()),
                "-database" => Some("--database".to_string()),
                "-dir" => Some("--dir".to_string()),
                "-organism" => Some("--organism".to_string()),
                "-ident_min" => Some("--ident_min".to_string()),
                "-coverage_min" => Some("--coverage_min".to_string()),
                "-threads" => Some("--threads".to_string()),
                "-plus" => Some("--plus".to_string()),
                "-report_common" => Some("--report_common".to_string()),
                "-report_all_equal" => Some("--report_all_equal".to_string()),
                "-quiet" => Some("--quiet".to_string()),
                "-debug" => Some("--debug".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-verbose" => Some("--verbose".to_string()),
                "-log" => Some("--log".to_string()),
                "-print_node" => Some("--print_node".to_string()),
                "-pgap" => Some("--pgap".to_string()),
                "-gpipe_org" => Some("--gpipe_org".to_string()),
                "-mutation_all" => Some("--mutation_all".to_string()),
                "-annotation_format" => Some("--annotation_format".to_string()),
                "-translation_table" => Some("--translation_table".to_string()),
                "-name" => Some("--name".to_string()),
                "-blast_bin" => Some("--blast_bin".to_string()),
                "-hmmer_bin" => Some("--hmmer_bin".to_string()),
                "-parm" => Some("--parm".to_string()),
                "-output" => Some("--output".to_string()),
                "-protein_output" => Some("--protein_output".to_string()),
                "-nucleotide_output" => Some("--nucleotide_output".to_string()),
                "-nucleotide_flank5_output" => Some("--nucleotide_flank5_output".to_string()),
                "-nucleotide_flank5_size" => Some("--nucleotide_flank5_size".to_string()),
                "-database_version" => Some("--database_version".to_string()),
                "-list_organisms" => Some("--list_organisms".to_string()),
                _ if value.starts_with("-protein=") => {
                    Some(value.replacen("-protein=", "--protein=", 1))
                }
                _ if value.starts_with("-nucleotide=") => {
                    Some(value.replacen("-nucleotide=", "--nucleotide=", 1))
                }
                _ if value.starts_with("-gff=") => Some(value.replacen("-gff=", "--gff=", 1)),
                _ if value.starts_with("-database=") => {
                    Some(value.replacen("-database=", "--database=", 1))
                }
                _ if value.starts_with("-dir=") => Some(value.replacen("-dir=", "--dir=", 1)),
                _ if value.starts_with("-organism=") => {
                    Some(value.replacen("-organism=", "--organism=", 1))
                }
                _ if value.starts_with("-ident_min=") => {
                    Some(value.replacen("-ident_min=", "--ident_min=", 1))
                }
                _ if value.starts_with("-coverage_min=") => {
                    Some(value.replacen("-coverage_min=", "--coverage_min=", 1))
                }
                _ if value.starts_with("-threads=") => {
                    Some(value.replacen("-threads=", "--threads=", 1))
                }
                _ if value.starts_with("-verbose=") => {
                    Some(value.replacen("-verbose=", "--verbose=", 1))
                }
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ if value.starts_with("-mutation_all=") => {
                    Some(value.replacen("-mutation_all=", "--mutation_all=", 1))
                }
                _ if value.starts_with("-annotation_format=") => {
                    Some(value.replacen("-annotation_format=", "--annotation_format=", 1))
                }
                _ if value.starts_with("-translation_table=") => {
                    Some(value.replacen("-translation_table=", "--translation_table=", 1))
                }
                _ if value.starts_with("-name=") => Some(value.replacen("-name=", "--name=", 1)),
                _ if value.starts_with("-blast_bin=") => {
                    Some(value.replacen("-blast_bin=", "--blast_bin=", 1))
                }
                _ if value.starts_with("-hmmer_bin=") => {
                    Some(value.replacen("-hmmer_bin=", "--hmmer_bin=", 1))
                }
                _ if value.starts_with("-parm=") => Some(value.replacen("-parm=", "--parm=", 1)),
                _ if value.starts_with("-output=") => {
                    Some(value.replacen("-output=", "--output=", 1))
                }
                _ if value.starts_with("-protein_output=") => {
                    Some(value.replacen("-protein_output=", "--protein_output=", 1))
                }
                _ if value.starts_with("-nucleotide_output=") => {
                    Some(value.replacen("-nucleotide_output=", "--nucleotide_output=", 1))
                }
                _ if value.starts_with("-nucleotide_flank5_output=") => Some(value.replacen(
                    "-nucleotide_flank5_output=",
                    "--nucleotide_flank5_output=",
                    1,
                )),
                _ if value.starts_with("-nucleotide_flank5_size=") => {
                    Some(value.replacen("-nucleotide_flank5_size=", "--nucleotide_flank5_size=", 1))
                }
                _ => None,
            },
            "check-fasta" => match value {
                "-aa" => Some("--aa".to_string()),
                "-hyphen" => Some("--hyphen".to_string()),
                "-ambig" => Some("--ambig".to_string()),
                "-ambig_max" => Some("--ambig_max".to_string()),
                "-stop_codon" => Some("--stop_codon".to_string()),
                "-len" => Some("--len".to_string()),
                "-out" => Some("--out".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-ambig_max=") => {
                    Some(value.replacen("-ambig_max=", "--ambig_max=", 1))
                }
                _ if value.starts_with("-len=") => Some(value.replacen("-len=", "--len=", 1)),
                _ if value.starts_with("-out=") => Some(value.replacen("-out=", "--out=", 1)),
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "extract" => match value {
                "-aa" => Some("--aa".to_string()),
                "-verbose" => Some("--verbose".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-verbose=") => {
                    Some(value.replacen("-verbose=", "--verbose=", 1))
                }
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "split-fasta" => match value {
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "mutate" => match value {
                "-aa" => Some("--aa".to_string()),
                "-orig" => Some("--orig".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "gff-check" => match value {
                "-gfftype" => Some("--gfftype".to_string()),
                "-prot" => Some("--prot".to_string()),
                "-dna" => Some("--dna".to_string()),
                "-lcl" => Some("--lcl".to_string()),
                "-gff_prot_match" => Some("--gff_prot_match".to_string()),
                "-gff_dna_match" => Some("--gff_dna_match".to_string()),
                "-verbose" => Some("--verbose".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-gfftype=") => {
                    Some(value.replacen("-gfftype=", "--gfftype=", 1))
                }
                _ if value.starts_with("-prot=") => Some(value.replacen("-prot=", "--prot=", 1)),
                _ if value.starts_with("-dna=") => Some(value.replacen("-dna=", "--dna=", 1)),
                _ if value.starts_with("-gff_prot_match=") => {
                    Some(value.replacen("-gff_prot_match=", "--gff_prot_match=", 1))
                }
                _ if value.starts_with("-gff_dna_match=") => {
                    Some(value.replacen("-gff_dna_match=", "--gff_dna_match=", 1))
                }
                _ if value.starts_with("-verbose=") => {
                    Some(value.replacen("-verbose=", "--verbose=", 1))
                }
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "dna-mutation" => match value {
                "-mutation_all" => Some("--mutation_all".to_string()),
                "-name" => Some("--name".to_string()),
                "-print_node" => Some("--print_node".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-mutation_all=") => {
                    Some(value.replacen("-mutation_all=", "--mutation_all=", 1))
                }
                _ if value.starts_with("-name=") => Some(value.replacen("-name=", "--name=", 1)),
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "disruption2-gene-symbol" => match value {
                "-gencode" => Some("--gencode".to_string()),
                "-prot_id_pos" => Some("--prot_id_pos".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-gencode=") => {
                    Some(value.replacen("-gencode=", "--gencode=", 1))
                }
                _ if value.starts_with("-prot_id_pos=") => {
                    Some(value.replacen("-prot_id_pos=", "--prot_id_pos=", 1))
                }
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            "amr-finder-update" => match value {
                "-database" => Some("--database".to_string()),
                "-d" => Some("-d".to_string()),
                "-blast_bin" => Some("--blast_bin".to_string()),
                "-hmmer_bin" => Some("--hmmer_bin".to_string()),
                "-force_update" => Some("--force_update".to_string()),
                "-debug" => Some("--debug".to_string()),
                "-q" => Some("-q".to_string()),
                "-quiet" => Some("--quiet".to_string()),
                _ if value.starts_with("-database=") => {
                    Some(value.replacen("-database=", "--database=", 1))
                }
                _ if value.starts_with("-blast_bin=") => {
                    Some(value.replacen("-blast_bin=", "--blast_bin=", 1))
                }
                _ if value.starts_with("-hmmer_bin=") => {
                    Some(value.replacen("-hmmer_bin=", "--hmmer_bin=", 1))
                }
                _ => None,
            },
            "update" => match value {
                "-protein" => Some("--protein".to_string()),
                "-nucleotide" => Some("--nucleotide".to_string()),
                "-database" => Some("--database".to_string()),
                "-blast_bin" => Some("--blast_bin".to_string()),
                "-hmmer_bin" => Some("--hmmer_bin".to_string()),
                "-force_update" => Some("--force_update".to_string()),
                "-debug" => Some("--debug".to_string()),
                "-qc" => Some("--qc".to_string()),
                "-quiet" => Some("--quiet".to_string()),
                _ if value.starts_with("-protein=") => {
                    Some(value.replacen("-protein=", "--protein=", 1))
                }
                _ if value.starts_with("-nucleotide=") => {
                    Some(value.replacen("-nucleotide=", "--nucleotide=", 1))
                }
                _ if value.starts_with("-database=") => {
                    Some(value.replacen("-database=", "--database=", 1))
                }
                _ if value.starts_with("-blast_bin=") => {
                    Some(value.replacen("-blast_bin=", "--blast_bin=", 1))
                }
                _ if value.starts_with("-hmmer_bin=") => {
                    Some(value.replacen("-hmmer_bin=", "--hmmer_bin=", 1))
                }
                _ => None,
            },
            "amr-finder-index" => match value {
                "-blast_bin" => Some("--blast_bin".to_string()),
                "-hmmer_bin" => Some("--hmmer_bin".to_string()),
                "-q" => Some("-q".to_string()),
                "-quiet" => Some("--quiet".to_string()),
                "-debug" => Some("--debug".to_string()),
                "-log" => Some("--log".to_string()),
                _ if value.starts_with("-blast_bin=") => {
                    Some(value.replacen("-blast_bin=", "--blast_bin=", 1))
                }
                _ if value.starts_with("-hmmer_bin=") => {
                    Some(value.replacen("-hmmer_bin=", "--hmmer_bin=", 1))
                }
                _ if value.starts_with("-log=") => Some(value.replacen("-log=", "--log=", 1)),
                _ => None,
            },
            _ => None,
        };
        if let Some(replacement) = replacement {
            *arg = OsString::from(replacement);
        }
    }
}
