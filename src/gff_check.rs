use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Result};

use crate::common::{Application, Key};
use crate::gff::{Annot, GffType};

pub const LOCUS_TAG_S: &str = "[locus_tag=";
pub const PRODIGAL_ID: &str = " ID=";
pub const NO_FILE: &str = "emptystring";

pub struct ThisApplication {
    application: Application,
}

impl ThisApplication {
    pub fn new() -> Self {
        Self {
            application: Application {
                description: "Check the correctness of a GFF file. Exit with an error if it is incorrect.",
                usage: "Usage: gff_check <gff> [-gfftype <type>] [-prot <file>] [-dna <file>] [-lcl] [-gff_prot_match <file>] [-gff_dna_match <file>]",
                positionals: 1,
                needs_arg: false,
                threads_used: false,
                key_args: vec![
                    Key {
                        name: "gfftype",
                        flag: false,
                        default_value: "genbank",
                        acronym: None,
                    },
                    Key {
                        name: "prot",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "dna",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "lcl",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "gff_prot_match",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "gff_dna_match",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                ],
            },
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn body(
    gff_name: &str,
    gff_type: GffType,
    prot_fname: &Path,
    dna_fname: &Path,
    prot_match_fname: &Path,
    dna_match_fname: &Path,
    lcl: bool,
    verbose: bool,
    out: &mut dyn Write,
) -> Result<()> {
    if lcl && gff_type != GffType::Pgap {
        bail!("-lcl requires type pgap");
    }

    if gff_name.ends_with(NO_FILE) {
        return Ok(());
    }

    let annot = Annot::from_gff(
        gff_name,
        gff_type,
        !prot_match_fname.as_os_str().is_empty(),
        lcl,
    )?;

    if !prot_fname.as_os_str().is_empty() {
        let mut gff_ids = Vec::with_capacity(10000);
        {
            let mut out = if prot_match_fname.as_os_str().is_empty() {
                None
            } else {
                Some(File::create(prot_match_fname)?)
            };
            let mut fasta_ids = Vec::with_capacity(gff_ids.capacity());
            let reader = BufReader::new(File::open(prot_fname)?);
            for line in reader.lines() {
                let line = line?;
                let line = line.trim_end();
                if line.is_empty() || !line.starts_with('>') {
                    continue;
                }
                let line_orig = line.to_string();
                let fasta_id = line[1..].split_whitespace().next().unwrap_or("");
                if fasta_id.is_empty() {
                    bail!("{}: No protein identifier in: {}", file!(), line_orig);
                }
                assert!(!fasta_id.contains(' '));
                fasta_ids.push(fasta_id.to_string());

                let mut gff_id = fasta_id.to_string();
                if out.is_some() {
                    match gff_type {
                        GffType::Genbank => {
                            let Some(pos) = line.find(LOCUS_TAG_S) else {
                                bail!(
                                    "{}: {:?} is not found in: {}",
                                    file!(),
                                    LOCUS_TAG_S,
                                    line_orig
                                );
                            };
                            gff_id = line[pos + LOCUS_TAG_S.len()..].to_string();
                            let Some(end) = gff_id.find(']') else {
                                bail!(
                                    "{}: ']' is not found after {:?} in: {}",
                                    file!(),
                                    LOCUS_TAG_S,
                                    line_orig
                                );
                            };
                            gff_id.truncate(end);
                        }
                        GffType::Microscope => {
                            let mut parts = gff_id.split('|');
                            let _acc = parts.next();
                            let Some(id) = parts.next() else {
                                bail!("{}: 'ID:' is not found in: {}", file!(), line_orig);
                            };
                            let Some(id) = id.strip_prefix("ID:") else {
                                bail!("{}: 'ID:' is not found in: {}", file!(), line_orig);
                            };
                            gff_id = id.to_string();
                        }
                        GffType::Prodigal => {
                            let Some(pos) = line.find(PRODIGAL_ID) else {
                                bail!(
                                    "{}: {:?} is not found in: {}",
                                    file!(),
                                    PRODIGAL_ID,
                                    line_orig
                                );
                            };
                            gff_id = line[pos + PRODIGAL_ID.len()..].to_string();
                            let Some(end) = gff_id.find(';') else {
                                bail!(
                                    "{}: ';' is not found after {:?} in: {}",
                                    file!(),
                                    PRODIGAL_ID,
                                    line_orig
                                );
                            };
                            gff_id.truncate(end);
                        }
                        _ => {}
                    }
                }
                if gff_id.contains(' ') {
                    bail!("{}: {:?} contains space", file!(), gff_id);
                }
                if gff_id.is_empty() {
                    bail!("{}: No protein identifier in: {}", file!(), line_orig);
                }
                gff_ids.push(gff_id.clone());
                if let Some(out) = &mut out {
                    writeln!(out, "{}\t{}", fasta_id, gff_id)?;
                }
            }

            let n = fasta_ids.len();
            fasta_ids.sort();
            fasta_ids.dedup();
            if fasta_ids.len() != n {
                bail!("{}: Duplicate FASTA ids", file!());
            }
            gff_ids.sort();
            for pair in gff_ids.windows(2) {
                if pair[0] == pair[1] {
                    bail!("{}: GFF identifier {:?} is not unique", file!(), pair[0]);
                }
            }
        }
        if verbose {
            writeln!(out, "# Proteins in GFF: {}", annot.prot2loci.len())?;
        }
        for seqid in &gff_ids {
            if !annot.prot2loci.contains_key(seqid) {
                bail!(
                    "{}: Protein FASTA id {:?} is not in the GFF file",
                    file!(),
                    seqid
                );
            }
        }
    }

    if !dna_fname.as_os_str().is_empty() {
        let mut contig_ids = Vec::with_capacity(10000);
        let mut gff_ids = Vec::with_capacity(10000);
        {
            let mut out = if dna_match_fname.as_os_str().is_empty() {
                None
            } else {
                Some(File::create(dna_match_fname)?)
            };
            let reader = BufReader::new(File::open(dna_fname)?);
            for line in reader.lines() {
                let line = line?;
                let line = line.trim_end();
                if line.is_empty() || !line.starts_with('>') {
                    continue;
                }
                let contig_id = line[1..].split_whitespace().next().unwrap_or("");
                assert!(!contig_id.contains(' '));

                let mut gff_id = contig_id.to_string();
                if out.is_some() && gff_type == GffType::Pseudomonasdb {
                    if let Some(pos) = gff_id.rfind('|') {
                        gff_id = gff_id[pos + 1..].to_string();
                    }
                }
                if gff_id.is_empty() {
                    bail!("{}: No contig identifier in:\n{}", file!(), line);
                }
                if lcl && !gff_id.starts_with("lcl|") {
                    bail!(
                        "{}: Contig identifier does not start with {:?}:\n{}",
                        file!(),
                        "lcl|",
                        line
                    );
                }
                gff_ids.push(gff_id.clone());
                contig_ids.push(contig_id.to_string());
                if let Some(out) = &mut out {
                    writeln!(out, "{}\t{}", contig_id, gff_id)?;
                }
            }
        }

        gff_ids.sort();
        for pair in gff_ids.windows(2) {
            if pair[0] == pair[1] {
                bail!(
                    "{}: DNA GFF identifier {:?} is not unique",
                    file!(),
                    pair[0]
                );
            }
        }
        contig_ids.sort();
        for pair in contig_ids.windows(2) {
            if pair[0] == pair[1] {
                bail!(
                    "{}: DNA contig identifier {:?} is not unique",
                    file!(),
                    pair[0]
                );
            }
        }

        let gff_id_set: BTreeSet<&str> = gff_ids.iter().map(String::as_str).collect();
        for loci in annot.prot2loci.values() {
            for cds in loci {
                if !gff_id_set.contains(cds.contig.as_str()) {
                    bail!(
                        "{}: GFF contig id {:?} is not in the DNA FASTA file",
                        file!(),
                        cds.contig
                    );
                }
            }
        }
    }

    Ok(())
}

pub fn main(argv: Vec<OsString>) -> Result<i32> {
    let mut app = ThisApplication::new();
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    let gff_type = GffType::name2type(&run.key_args["gfftype"])?;
    let prot_fname = PathBuf::from(&run.key_args["prot"]);
    let dna_fname = PathBuf::from(&run.key_args["dna"]);
    let prot_match_fname = PathBuf::from(&run.key_args["gff_prot_match"]);
    let dna_match_fname = PathBuf::from(&run.key_args["gff_dna_match"]);
    let verbose = run.key_args["verbose"].parse::<usize>().unwrap_or(0) > 0;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    body(
        &run.positional_args[0].to_string_lossy(),
        gff_type,
        &prot_fname,
        &dna_fname,
        &prot_match_fname,
        &dna_match_fname,
        run.key_args["lcl"] == "true",
        verbose,
        &mut out,
    )?;
    Ok(0)
}
