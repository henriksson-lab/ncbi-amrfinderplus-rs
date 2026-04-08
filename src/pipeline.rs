// Main AMRFinder pipeline orchestrator — port of amrfinder.cpp
// Coordinates BLAST, HMM, and report generation

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Result};

use crate::fasta_utils::{self, FastaCheckOpts};

/// Pipeline configuration
#[allow(clippy::struct_excessive_bools)]
pub struct PipelineConfig {
    pub protein: Option<PathBuf>,
    pub nucleotide: Option<PathBuf>,
    pub gff: Option<PathBuf>,
    pub database: PathBuf,
    pub organism: String,
    pub ident_min: f64,
    pub coverage_min: f64,
    pub threads: usize,
    pub plus: bool,
    pub print_node: bool,
    pub mutation_all: Option<PathBuf>,
    pub annotation_format: String,
    pub translation_table: u32,
    pub name: String,
    pub blast_bin: String,
    pub hmmer_bin: String,
    pub output: Option<PathBuf>,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        PipelineConfig {
            protein: None,
            nucleotide: None,
            gff: None,
            database: PathBuf::new(),
            organism: String::new(),
            ident_min: -1.0,
            coverage_min: 0.5,
            threads: 4,
            plus: false,
            print_node: false,
            mutation_all: None,
            annotation_format: "genbank".to_string(),
            translation_table: 11,
            name: String::new(),
            blast_bin: String::new(),
            hmmer_bin: String::new(),
            output: None,
        }
    }
}

impl PipelineConfig {
    fn find_prog(&self, name: &str) -> PathBuf {
        let dir = match name {
            "blastp" | "blastn" | "blastx" | "tblastn" | "makeblastdb" => &self.blast_bin,
            "hmmsearch" | "hmmpress" => &self.hmmer_bin,
            _ => "",
        };
        if dir.is_empty() {
            PathBuf::from(name)
        } else {
            Path::new(dir).join(name)
        }
    }
}

/// Run the AMRFinder pipeline
/// This initially shells out to external tools (BLAST, HMM, amr_report)
/// matching the C++ behavior exactly
pub fn run_pipeline(config: &PipelineConfig) -> Result<String> {
    // Validate inputs
    if config.protein.is_none() && config.nucleotide.is_none() {
        bail!("Either --protein or --nucleotide must be provided");
    }

    if config.database.as_os_str().is_empty() {
        bail!("Database directory must be specified");
    }

    if !config.database.exists() {
        bail!("Database directory does not exist: {}", config.database.display());
    }

    // Create temp directory
    let tmp_dir = tempfile::tempdir()?;
    let tmp = tmp_dir.path();

    let db = config.database.canonicalize()?.to_str().unwrap().to_string();

    let mut amr_report_blastp = String::new();
    let mut amr_report_blastx = String::new();

    // Process protein input
    if let Some(ref prot_path) = config.protein {
        if !prot_path.exists() || fs::metadata(prot_path)?.len() == 0 {
            bail!("Protein FASTA file is empty or missing: {}", prot_path.display());
        }

        // fasta_check
        let (n_prot, _prot_len_max, prot_len_total) = fasta_utils::fasta_check(
            &FastaCheckOpts {
                fasta_path: prot_path,
                aa: true,
                hyphen: false,
                ambig: false,
                ambig_max: 20,
                stop_codon: true,
                len_path: None,
                out_path: None,
            },
        )?;

        // Run blastp
        let blastp_out = tmp.join("blastp");
        let blastp_prog = config.find_prog("blastp");
        let blastp_threads = std::cmp::min(
            config.threads,
            std::cmp::min(n_prot, prot_len_total / 10000 + 1),
        );

        // Run blastp and hmmsearch in parallel
        let hmmsearch_out = tmp.join("hmmsearch");
        let dom_out = tmp.join("dom");
        let hmmsearch_prog = config.find_prog("hmmsearch");
        let hmm_lib = format!("{}/AMR.LIB", db);
        let hmm_cpu = config.threads.saturating_sub(1).to_string();
        let prot_str = prot_path.to_str().unwrap().to_string();
        let hmmsearch_out_str = hmmsearch_out.to_str().unwrap().to_string();
        let dom_out_str = dom_out.to_str().unwrap().to_string();

        let hmm_handle = std::thread::spawn(move || {
            Command::new(&hmmsearch_prog)
                .args([
                    "--tblout", &hmmsearch_out_str,
                    "--noali",
                    "--domtblout", &dom_out_str,
                    "--cut_tc",
                    "-Z", "10000",
                    "--cpu", &hmm_cpu,
                    &hmm_lib,
                    &prot_str,
                ])
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .status()
        });

        let blastp_output = Command::new(&blastp_prog)
            .args([
                "-query", prot_path.to_str().unwrap(),
                "-db", &format!("{}/AMRProt.fa", db),
                "-comp_based_stats", "0",
                "-seg", "no",
                "-max_target_seqs", "10000",
                "-dbsize", "10000",
                "-evalue", "1e-10",
                "-word_size", "5",
                "-task", "blastp-fast",
                "-num_threads", &blastp_threads.to_string(),
                "-outfmt", "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
                "-out", blastp_out.to_str().unwrap(),
            ])
            .output()?;

        if !blastp_output.status.success() {
            let stderr = String::from_utf8_lossy(&blastp_output.stderr);
            bail!("blastp failed: {}", stderr);
        }

        // Wait for hmmsearch
        let hmm_result = hmm_handle.join().map_err(|_| anyhow::anyhow!("hmmsearch thread panicked"))??;
        if !hmm_result.success() {
            bail!("hmmsearch failed");
        }

        amr_report_blastp = format!(
            "-blastp {} -hmmsearch {} -hmmdom {}",
            blastp_out.to_str().unwrap(),
            hmmsearch_out.to_str().unwrap(),
            dom_out.to_str().unwrap(),
        );

        if let Some(ref gff_path) = config.gff {
            amr_report_blastp += &format!(
                " -gff {} -gfftype {}",
                gff_path.to_str().unwrap(),
                config.annotation_format,
            );
        }
    } else {
        // Create empty files
        File::create(tmp.join("blastp"))?;
        File::create(tmp.join("hmmsearch"))?;
        File::create(tmp.join("dom"))?;
    }

    // Process nucleotide input
    if let Some(ref dna_path) = config.nucleotide {
        if !dna_path.exists() || fs::metadata(dna_path)?.len() == 0 {
            bail!("Nucleotide FASTA file is empty or missing: {}", dna_path.display());
        }

        // fasta_check for DNA
        let (n_dna, dna_len_max, dna_len_total) = fasta_utils::fasta_check(
            &FastaCheckOpts {
                fasta_path: dna_path,
                aa: false,
                hyphen: true,
                ambig: true,
                ambig_max: 0,
                stop_codon: false,
                len_path: Some(&tmp.join("len")),
                out_path: None,
            },
        )?;

        // Choose blastx or tblastn based on sequence length
        let use_tblastn = dna_len_max > 100_000;
        let blast_prog_name = if use_tblastn { "tblastn" } else { "blastx" };
        let blast_prog = config.find_prog(blast_prog_name);

        let blastx_out = tmp.join("blastx");

        if use_tblastn {
            // tblastn: protein query vs translated nucleotide DB
            let blast_status = Command::new(&blast_prog)
                .args([
                    "-query", &format!("{}/AMRProt.fa", db),
                    "-subject", dna_path.to_str().unwrap(),
                    "-comp_based_stats", "0",
                    "-seg", "no",
                    "-max_target_seqs", "10000",
                    "-dbsize", "10000",
                    "-evalue", "1e-10",
                    "-word_size", "5",
                    "-task", "tblastn-fast",
                    "-window_size", "15",
                    "-threshold", "100",
                    "-db_gencode", &config.translation_table.to_string(),
                    "-num_threads", &config.threads.to_string(),
                    "-outfmt", "6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq",
                    "-out", blastx_out.to_str().unwrap(),
                ])
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .status()?;

            if !blast_status.success() {
                bail!("{} failed", blast_prog_name);
            }
        } else {
            // blastx: nucleotide query vs protein DB
            let blast_threads = std::cmp::min(
                config.threads,
                std::cmp::min(n_dna, dna_len_total / 10002 + 1),
            );
            let blast_status = Command::new(&blast_prog)
                .args([
                    "-query", dna_path.to_str().unwrap(),
                    "-db", &format!("{}/AMRProt.fa", db),
                    "-comp_based_stats", "0",
                    "-seg", "no",
                    "-max_target_seqs", "10000",
                    "-dbsize", "10000",
                    "-evalue", "1e-10",
                    "-word_size", "5",
                    "-query_gencode", &config.translation_table.to_string(),
                    "-num_threads", &blast_threads.to_string(),
                    "-outfmt", "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
                    "-out", blastx_out.to_str().unwrap(),
                ])
                .stderr(std::process::Stdio::null())
                .stdout(std::process::Stdio::null())
                .status()?;

            if !blast_status.success() {
                bail!("{} failed", blast_prog_name);
            }
        }

        amr_report_blastx = format!("-blastx {}", blastx_out.to_str().unwrap());

        // DNA mutation search (blastn)
        if !config.organism.is_empty() {
            let dna_db = format!("{}/AMR_DNA-{}.fa", db, config.organism);
            if Path::new(&dna_db).exists() {
                let blastn_out = tmp.join("blastn");
                let blastn_prog = config.find_prog("blastn");

                let blastn_status = Command::new(&blastn_prog)
                    .args([
                        "-query", dna_path.to_str().unwrap(),
                        "-db", &dna_db,
                        "-evalue", "1e-20",
                        "-dust", "no",
                        "-max_target_seqs", "10000",
                        "-num_threads", &config.threads.to_string(),
                        "-outfmt", "6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq",
                        "-out", blastn_out.to_str().unwrap(),
                    ])
                    .stderr(std::process::Stdio::null())
                    .stdout(std::process::Stdio::null())
                    .status()?;

                if !blastn_status.success() {
                    bail!("blastn failed");
                }
            }
        }
    }

    // Try C++ amr_report first; if not available, use Rust implementation
    let amr_report_bin = which_amr_report();

    if amr_report_bin.is_err() {
        // C++ amr_report not available — use Rust implementation
        let rust_result = run_rust_amr_report(config, tmp, &db, &amr_report_blastp, &amr_report_blastx)?;
        let result = sort_tsv_output(&rust_result, config)?;
        if let Some(ref output_path) = config.output {
            fs::write(output_path, &result)?;
        }
        return Ok(result);
    }

    let amr_report_bin = amr_report_bin.unwrap();

    let mut amr_report_args: Vec<String> = Vec::new();
    amr_report_args.push("-fam".to_string());
    amr_report_args.push(format!("{}/fam.tsv", db));

    // Parse the accumulated args
    for part in amr_report_blastp.split_whitespace() {
        amr_report_args.push(part.to_string());
    }
    if !amr_report_blastx.is_empty() {
        for part in amr_report_blastx.split_whitespace() {
            amr_report_args.push(part.to_string());
        }
    }

    // Add DNA length file if we have nucleotide input
    if config.nucleotide.is_some() {
        let len_file = tmp.join("len");
        if len_file.exists() {
            amr_report_args.extend([
                "-dna_len".to_string(),
                len_file.to_str().unwrap().to_string(),
            ]);
        }
    }

    amr_report_args.extend([
        "-organism".to_string(), config.organism.clone(),
        "-mutation".to_string(), format!("{}/AMRProt-mutation.tsv", db),
        "-susceptible".to_string(), format!("{}/AMRProt-susceptible.tsv", db),
    ]);

    if !config.plus {
        amr_report_args.push("-core".to_string());
    }
    if config.print_node {
        amr_report_args.push("-print_node".to_string());
    }
    if config.ident_min >= 0.0 {
        amr_report_args.extend([
            "-ident_min".to_string(),
            config.ident_min.to_string(),
        ]);
    }
    amr_report_args.extend([
        "-coverage_min".to_string(),
        config.coverage_min.to_string(),
    ]);
    if !config.name.is_empty() {
        amr_report_args.extend([
            "-name".to_string(),
            config.name.clone(),
        ]);
    }

    // Suppress common proteins for organisms
    let suppress_prot_file = tmp.join("suppress_prot");
    if !config.organism.is_empty() {
        let suppress_tsv = format!("{}/AMRProt-suppress.tsv", db);
        if Path::new(&suppress_tsv).exists() {
            let mut suppress_out = File::create(&suppress_prot_file)?;
            let file = File::open(&suppress_tsv)?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line?;
                if line.starts_with('#') {
                    continue;
                }
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() >= 2 && fields[0] == config.organism {
                    writeln!(suppress_out, "{}", fields[1])?;
                }
            }
            drop(suppress_out);
            amr_report_args.extend([
                "-suppress_prot".to_string(),
                suppress_prot_file.to_str().unwrap().to_string(),
            ]);
        }
    }

    // Force CDS report for DNA + organism combination
    if config.nucleotide.is_some() && !config.organism.is_empty() {
        amr_report_args.push("-force_cds_report".to_string());
    }

    let amr_report_output = Command::new(&amr_report_bin)
        .args(&amr_report_args)
        .output()?;

    if !amr_report_output.status.success() {
        let stderr = String::from_utf8_lossy(&amr_report_output.stderr);
        bail!("amr_report failed: {}", stderr);
    }

    let mut raw_result = String::from_utf8(amr_report_output.stdout)?;

    // Run dna_mutation for DNA SNP detection (blastn results)
    if config.nucleotide.is_some() && !config.organism.is_empty() {
        let blastn_file = tmp.join("blastn");
        let dna_tsv = format!("{}/AMR_DNA-{}.tsv", db, config.organism);
        if blastn_file.exists() && Path::new(&dna_tsv).exists() {
            let dna_mutation_bin = which_binary("dna_mutation")?;
            let mut dna_mut_args = vec![
                blastn_file.to_str().unwrap().to_string(),
                dna_tsv,
                config.organism.clone(),
            ];

            if config.mutation_all.is_some() {
                dna_mut_args.extend([
                    "-mutation_all".to_string(),
                    tmp.join("mutation_all.dna").to_str().unwrap().to_string(),
                ]);
            }

            if !config.name.is_empty() {
                dna_mut_args.extend(["-name".to_string(), config.name.clone()]);
            }
            if config.print_node {
                dna_mut_args.push("-print_node".to_string());
            }

            let dna_mut_output = Command::new(&dna_mutation_bin)
                .args(&dna_mut_args)
                .output()?;

            if dna_mut_output.status.success() {
                let snp_output = String::from_utf8(dna_mut_output.stdout)?;
                // Append SNP results (skip header line)
                for (i, line) in snp_output.lines().enumerate() {
                    if i > 0 && !line.is_empty() {
                        raw_result.push_str(line);
                        raw_result.push('\n');
                    }
                }
            }
        }
    }

    // Post-processing: sort by sort columns (matching C++ amrfinder.cpp behavior)
    let result = sort_tsv_output(&raw_result, config)?;

    // Handle output
    if let Some(ref output_path) = config.output {
        fs::write(output_path, &result)?;
    }

    Ok(result)
}

/// Run the Rust amr_report implementation
fn run_rust_amr_report(
    config: &PipelineConfig,
    tmp: &Path,
    db: &str,
    amr_report_blastp: &str,
    amr_report_blastx: &str,
) -> Result<String> {
    use crate::amr_reportcli::{AmrReportConfig, run_amr_report};

    // Parse blastp/hmmsearch/dom file paths from the accumulated args
    let mut blastp_path = None;
    let mut hmmsearch_path = None;
    let mut hmmdom_path = None;
    let mut gff_path = None;
    let mut gff_type = "genbank";

    let parts: Vec<&str> = amr_report_blastp.split_whitespace().collect();
    let mut i = 0;
    while i < parts.len() {
        match parts[i] {
            "-blastp" => { i += 1; if i < parts.len() { blastp_path = Some(PathBuf::from(parts[i])); } }
            "-hmmsearch" => { i += 1; if i < parts.len() { hmmsearch_path = Some(PathBuf::from(parts[i])); } }
            "-hmmdom" => { i += 1; if i < parts.len() { hmmdom_path = Some(PathBuf::from(parts[i])); } }
            "-gff" => { i += 1; if i < parts.len() { gff_path = Some(PathBuf::from(parts[i])); } }
            "-gfftype" => { i += 1; if i < parts.len() { gff_type = parts[i]; } }
            _ => {}
        }
        i += 1;
    }

    let mut blastx_path = None;
    let bx_parts: Vec<&str> = amr_report_blastx.split_whitespace().collect();
    i = 0;
    while i < bx_parts.len() {
        if bx_parts[i] == "-blastx" { i += 1; if i < bx_parts.len() { blastx_path = Some(PathBuf::from(bx_parts[i])); } }
        i += 1;
    }

    let cds_exist = gff_path.is_some() || blastx_path.is_some() || (config.nucleotide.is_some() && !config.organism.is_empty());

    let fam_file = PathBuf::from(format!("{}/fam.tsv", db));
    let mutation_file = PathBuf::from(format!("{}/AMRProt-mutation.tsv", db));
    let susceptible_file = PathBuf::from(format!("{}/AMRProt-susceptible.tsv", db));

    let report_config = AmrReportConfig {
        fam_file: &fam_file,
        blastp_file: blastp_path.as_deref(),
        blastx_file: blastx_path.as_deref(),
        hmmsearch_file: hmmsearch_path.as_deref(),
        hmmdom_file: hmmdom_path.as_deref(),
        gff_file: gff_path.as_deref(),
        gff_type,
        organism: &config.organism,
        mutation_file: Some(&mutation_file),
        susceptible_file: Some(&susceptible_file),
        coverage_min: config.coverage_min,
        ident_min: config.ident_min,
        print_node: config.print_node,
        report_core_only: !config.plus,
        cds_exist,
    };

    let mut output = Vec::new();
    run_amr_report(&report_config, &mut output)?;

    // Add dna_mutation results if applicable
    let mut raw_result = String::from_utf8(output)?;

    // DNA mutation detection using Rust implementation
    if config.nucleotide.is_some() && !config.organism.is_empty() {
        let blastn_file = tmp.join("blastn");
        let dna_tsv_path = PathBuf::from(format!("{}/AMR_DNA-{}.tsv", db, config.organism));
        if blastn_file.exists() && dna_tsv_path.exists() {
            let mut dna_output = Vec::new();
            if crate::dna_mutation::run_dna_mutation(
                &blastn_file,
                &dna_tsv_path,
                &config.organism,
                config.print_node,
                &config.name,
                &mut dna_output,
                None,
            ).is_ok() {
                let snp_output = String::from_utf8_lossy(&dna_output);
                for (i, line) in snp_output.lines().enumerate() {
                    if i > 0 && !line.is_empty() {
                        raw_result.push_str(line);
                        raw_result.push('\n');
                    }
                }
            }
        }
    }

    Ok(raw_result)
}

/// Find a C++ binary (amr_report, dna_mutation, etc.)
fn which_binary(name: &str) -> Result<PathBuf> {
    // Look in the amr/ source directory
    let amr_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr").join(name);
    if amr_dir.exists() {
        return Ok(amr_dir);
    }

    // Look in the same directory as current executable
    if let Ok(exe) = std::env::current_exe() {
        let dir = exe.parent().unwrap_or(Path::new("."));
        let candidate = dir.join(name);
        if candidate.exists() {
            return Ok(candidate);
        }
    }

    // Try PATH
    let output = Command::new("which")
        .arg(name)
        .output()?;
    if output.status.success() {
        let path = String::from_utf8(output.stdout)?.trim().to_string();
        return Ok(PathBuf::from(path));
    }

    bail!("{} binary not found", name);
}

/// Find the amr_report binary
fn which_amr_report() -> Result<PathBuf> {
    // Look in the same directory as the current executable
    if let Ok(exe) = std::env::current_exe() {
        let dir = exe.parent().unwrap_or(Path::new("."));
        let candidate = dir.join("amr_report");
        if candidate.exists() {
            return Ok(candidate);
        }
    }

    // Look in the source amr/ directory
    let amr_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amr/amr_report");
    if amr_dir.exists() {
        return Ok(amr_dir);
    }

    // Try PATH
    let output = Command::new("which")
        .arg("amr_report")
        .output()?;
    if output.status.success() {
        let path = String::from_utf8(output.stdout)?.trim().to_string();
        return Ok(PathBuf::from(path));
    }

    bail!("amr_report binary not found");
}

/// Sort the TSV output to match C++ amrfinder's post-processing
/// Sort columns: Contig id, Start, Stop, Strand, Protein id, Element symbol
fn sort_tsv_output(tsv: &str, config: &PipelineConfig) -> Result<String> {
    let mut lines: Vec<&str> = tsv.lines().collect();
    if lines.len() <= 1 {
        return Ok(tsv.to_string());
    }

    // First line is header
    let header = lines[0];
    let header_fields: Vec<&str> = header.split('\t').collect();

    // Find column indices for sort keys
    let find_col = |name: &str| -> Option<usize> {
        // Skip '#' prefix in header
        let h0 = header_fields[0].strip_prefix('#').unwrap_or(header_fields[0]);
        if h0 == name {
            return Some(0);
        }
        header_fields.iter().position(|&f| f == name)
    };

    let contig_col = find_col(crate::columns::CONTIG_COL_NAME);
    let start_col = find_col(crate::columns::START_COL_NAME);
    let stop_col = find_col(crate::columns::STOP_COL_NAME);
    let strand_col = find_col(crate::columns::STRAND_COL_NAME);
    let prot_col = find_col(crate::columns::PROT_COL_NAME);
    let gene_col = find_col(crate::columns::GENESYMBOL_COL_NAME);

    let has_cds = config.nucleotide.is_some() || config.gff.is_some();

    let data_lines = &mut lines[1..];

    // Remove duplicates first
    data_lines.sort();
    let mut unique: Vec<&str> = Vec::with_capacity(data_lines.len());
    for &line in data_lines.iter() {
        if unique.last().is_none_or(|&last| last != line) {
            unique.push(line);
        }
    }

    // Sort by the appropriate columns
    unique.sort_by(|a, b| {
        let a_fields: Vec<&str> = a.split('\t').collect();
        let b_fields: Vec<&str> = b.split('\t').collect();

        let get_field = |fields: &Vec<&str>, col: Option<usize>| -> String {
            col.and_then(|c| fields.get(c).copied())
                .unwrap_or("")
                .to_string()
        };

        let mut ord = std::cmp::Ordering::Equal;

        if has_cds {
            ord = ord.then(get_field(&a_fields, contig_col).cmp(&get_field(&b_fields, contig_col)));
            if let Some(c) = start_col {
                let a_val: i64 = a_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                let b_val: i64 = b_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                ord = ord.then(a_val.cmp(&b_val));
            }
            if let Some(c) = stop_col {
                let a_val: i64 = a_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                let b_val: i64 = b_fields.get(c).and_then(|s| s.parse().ok()).unwrap_or(0);
                ord = ord.then(a_val.cmp(&b_val));
            }
            ord = ord.then(get_field(&a_fields, strand_col).cmp(&get_field(&b_fields, strand_col)));
        }

        ord = ord.then(get_field(&a_fields, prot_col).cmp(&get_field(&b_fields, prot_col)));
        ord = ord.then(get_field(&a_fields, gene_col).cmp(&get_field(&b_fields, gene_col)));

        ord
    });

    let mut result = String::new();
    result.push_str(header);
    result.push('\n');
    for line in unique {
        if !line.is_empty() {
            result.push_str(line);
            result.push('\n');
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_config_default() {
        let config = PipelineConfig::default();
        assert_eq!(config.threads, 4);
        assert_eq!(config.coverage_min, 0.5);
        assert!(!config.plus);
    }
}
