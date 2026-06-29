use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::alignment::AmrMutation;
use crate::common::{Application, Key};
use crate::seq::{EXT_DNA_ALPHABET, EXT_TERM_PEPTIDE_ALPHABET};

pub struct ThisApplication {
    application: Application,
}

impl ThisApplication {
    pub fn new() -> Self {
        Self {
            application: Application {
                description: "Mutate a FASTA file",
                usage: "Usage: mutate <in> <mut> [-aa] [-orig]",
                positionals: 2,
                needs_arg: false,
                threads_used: false,
                key_args: vec![
                    Key {
                        name: "aa",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "orig",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                ],
            },
        }
    }
}

pub fn body(
    in_fname: &Path,
    mut_fname: &Path,
    aa: bool,
    orig: bool,
    out: &mut dyn Write,
) -> Result<()> {
    let mut id2mutation: BTreeMap<String, Vec<AmrMutation>> = BTreeMap::new();
    {
        let reader = BufReader::new(File::open(mut_fname)?);
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 4 {
                bail!("Malformed mutation row: {line}");
            }
            let seq_id = fields[0].to_string();
            let pos: usize = fields[1].parse()?;
            let mutation_std = fields[2];
            let mutation_report = fields[3];
            assert!(!mutation_report.is_empty());
            let mutation = AmrMutation::new(pos, mutation_std, mutation_report, "X", "X", "X");
            mutation.qc();
            id2mutation.entry(seq_id).or_default().push(mutation);
        }
    }

    let reader = BufReader::new(File::open(in_fname)?);
    let mut lines = reader.lines();
    let mut line_num = 0usize;
    let mut current = match lines.next() {
        Some(line) => {
            line_num += 1;
            Some(line?)
        }
        None => None,
    };
    while let Some(header) = current.take() {
        if header.is_empty() || !header.starts_with('>') {
            bail!("Error in Multifasta, line {}", line_num + 1);
        }

        let seq_name = header[1..].replace('\t', " ");
        if seq_name.trim().is_empty() {
            bail!("Blank sequence name");
        }

        let mut seq0 = String::with_capacity(if aa { 1000 } else { 100000 });
        loop {
            current = match lines.next() {
                Some(line) => {
                    line_num += 1;
                    Some(line?)
                }
                None => None,
            };
            let Some(line) = current.as_ref() else {
                break;
            };
            if line.trim().is_empty() || line.starts_with('>') {
                break;
            }
            seq0.push_str(line.trim_end());
        }

        while current.as_ref().is_some_and(|line| line.trim().is_empty()) {
            current = match lines.next() {
                Some(line) => {
                    line_num += 1;
                    Some(line?)
                }
                None => None,
            };
        }

        if aa {
            seq0 = seq0.to_ascii_uppercase();
        } else {
            seq0 = seq0.to_ascii_lowercase();
        }
        seq0.retain(|c| c != '-');

        let alphabet = if aa {
            EXT_TERM_PEPTIDE_ALPHABET
        } else {
            EXT_DNA_ALPHABET
        };
        if let Some((i, c)) = seq0
            .chars()
            .enumerate()
            .find(|(_, c)| !alphabet.contains(*c))
        {
            let seq_id = seq_name.split_whitespace().next().unwrap_or("");
            bail!(
                "{seq_name}\nBad sequence character in {:?}: ASCII={} pos={}",
                seq_id,
                c as u32,
                i + 1
            );
        }
        if orig {
            writeln!(out, ">{seq_name}")?;
            for chunk in seq0.as_bytes().chunks(80) {
                writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
            }
        }
        let seq_id = seq_name.split_whitespace().next().unwrap_or("");
        if let Some(mutations) = id2mutation.get(seq_id) {
            for mutation in mutations {
                let mut seq1 = seq0.clone();
                mutation
                    .apply(&mut seq1)
                    .with_context(|| seq_name.clone())?;
                if !aa {
                    seq1 = seq1.to_ascii_lowercase();
                }
                if let Some((i, c)) = seq1
                    .chars()
                    .enumerate()
                    .find(|(_, c)| !alphabet.contains(*c))
                {
                    let seq_id = seq_name.split_whitespace().next().unwrap_or("");
                    bail!(
                        "{seq_name}\nBad sequence character in {:?}: ASCII={} pos={}",
                        seq_id,
                        c as u32,
                        i + 1
                    );
                }
                writeln!(
                    out,
                    ">{}:{}:{}",
                    seq_name,
                    mutation.pos_real + 1,
                    mutation.gene_mutation
                )?;
                for chunk in seq1.as_bytes().chunks(80) {
                    writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
                }
            }
        }
    }

    Ok(())
}

pub fn main(argv: Vec<OsString>, out: &mut dyn Write) -> Result<i32> {
    let mut app = ThisApplication::new();
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    let in_fname = PathBuf::from(&run.positional_args[0]);
    let mut_fname = PathBuf::from(&run.positional_args[1]);

    body(
        &in_fname,
        &mut_fname,
        run.key_args["aa"] == "true",
        run.key_args["orig"] == "true",
        out,
    )?;
    Ok(0)
}
