// DNA-level point mutation detection — port of dna_mutation.cpp

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Result};

use crate::alignment::{Alignment, AmrMutation};
use crate::columns;
use crate::common::{Application, Key};
use crate::seq::Hsp;
use crate::tsv::TsvOut;

const FLANKING_LEN: usize = 200;

/// BlastnAlignment — DNA-level alignment with mutation detection
struct BlastnAlignment {
    alignment: Alignment,
    organism: String,
    ref_accession_frag: String,
    product: String,
    gene: String,
    ref_mutations: Vec<AmrMutation>,
}

struct Batch {
    blast_als: Vec<BlastnAlignment>,
    accession2mutations: HashMap<String, Vec<AmrMutation>>,
}

pub struct ThisApplication {
    application: Application,
}

impl ThisApplication {
    pub fn new() -> Self {
        Self {
            application: Application {
                description: "Find mutations at DNA level and report in the format of amr_report.cpp",
                usage: "Usage: dna_mutation <blastn> <mutation> <organism> [-mutation_all <file>] [-name <text>] [-print_node]",
                positionals: 3,
                needs_arg: false,
                threads_used: false,
                key_args: vec![
                    Key {
                        name: "mutation_all",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "name",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "print_node",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                ],
            },
        }
    }
}

impl BlastnAlignment {
    fn parse(
        line: &str,
        organism: &str,
        accession2mutations: &HashMap<String, Vec<AmrMutation>>,
    ) -> Result<Self> {
        let mut hsp = Hsp::from_blast_line(line, false, false, false, false, false)?;
        hsp.qseq.make_ascii_uppercase();
        hsp.sseq.make_ascii_uppercase();

        let organism = organism.replace('_', " ");

        // qseqid = accession@gene_name@gene_symbol@offset:start-stop
        let qseqid = &hsp.qseqid;
        let (ref_accession_frag, product, gene) = if !qseqid.is_empty() {
            let mut rest = qseqid.as_str();
            let Some(pos) = rest.find('@') else {
                bail!("{line}\nMissing @ in qseqid");
            };
            let accession = &rest[..pos];
            rest = &rest[pos + 1..];
            let Some(pos) = rest.find('@') else {
                bail!("{line}\nMissing second @ in qseqid");
            };
            let product = rest[..pos].replace('_', " ");
            rest = &rest[pos + 1..];
            let Some(pos) = rest.find(':') else {
                bail!("{line}\nMissing : in qseqid");
            };
            let gene = rest[..pos].to_string();
            rest = &rest[pos + 1..];

            (format!("{accession}:{rest}"), product, gene)
        } else {
            bail!("{line}\nEmpty qseqid");
        };

        let ref_mutations = accession2mutations.get(qseqid).cloned().unwrap_or_default();
        let mut alignment = Alignment {
            hsp,
            ref_mutation: AmrMutation::default(),
            seq_changes: Vec::new(),
        };
        if !ref_mutations.is_empty() {
            alignment.set_seq_changes(&ref_mutations, FLANKING_LEN)?;
        }

        Ok(BlastnAlignment {
            alignment,
            organism,
            ref_accession_frag,
            product,
            gene,
            ref_mutations,
        })
    }

    fn good(&self) -> bool {
        let min_len = std::cmp::min(self.alignment.hsp.qlen, 2 * FLANKING_LEN + 1);
        self.alignment.hsp.sseq.len() >= min_len
    }

    fn report(
        &self,
        td: &mut TsvOut,
        mutation_all: bool,
        print_node: bool,
        name: &str,
    ) -> Result<()> {
        for change in &self.alignment.seq_changes {
            let mutations: Vec<Option<&AmrMutation>> = if change.mutations.is_empty() {
                vec![None]
            } else {
                change
                    .mutations
                    .iter()
                    .map(|&idx| Some(&self.ref_mutations[idx]))
                    .collect()
            };
            for mutation in mutations {
                if !mutation_all
                    && (change.empty() || mutation.is_none() || change.replacement.is_some())
                {
                    continue;
                }
                if change.empty() && mutation.is_none() {
                    continue;
                }

                let gene_symbol = if let Some(mutation) = mutation {
                    if change.empty() {
                        mutation.wildtype()
                    } else {
                        mutation.gene_mutation.clone()
                    }
                } else {
                    format!("{}_{}", self.gene, change.get_mutation_str())
                };

                let elem_name = if let Some(mutation) = mutation {
                    if change.empty() {
                        format!("{} {} [WILDTYPE]", self.organism, self.product)
                    } else {
                        mutation.name.clone()
                    }
                } else {
                    format!("{} {} [UNKNOWN]", self.organism, self.product)
                };

                if !name.is_empty() {
                    td.write_field(&name)?;
                }
                td.write_field(&columns::NA)?; // Protein id
                td.write_field(&self.alignment.hsp.sseqid)?; // Contig
                td.write_field(&(self.alignment.hsp.s_int.start + 1))?; // Start
                td.write_field(&self.alignment.hsp.s_int.stop)?; // Stop
                td.write_field(
                    &crate::seq::strand2char(self.alignment.hsp.s_int.strand).to_string(),
                )?; // Strand
                td.write_field(&gene_symbol)?; // Element symbol
                td.write_field(&elem_name)?; // Element name
                td.write_field(&"core")?; // Scope
                td.write_field(&"AMR")?; // Type
                td.write_field(&"POINT")?; // Subtype
                if let Some(mutation) = mutation {
                    td.write_field(&if mutation.class.is_empty() {
                        columns::NA
                    } else {
                        &mutation.class
                    })?;
                    td.write_field(&if mutation.subclass.is_empty() {
                        columns::NA
                    } else {
                        &mutation.subclass
                    })?;
                } else {
                    td.write_field(&columns::NA)?;
                    td.write_field(&columns::NA)?;
                }
                td.write_field(&"POINTN")?; // Method
                td.write_field(&self.alignment.hsp.s_int.len())?; // Target length
                td.write_field(&self.alignment.hsp.qlen)?; // Reference length
                td.write_field(&format!(
                    "{:.2}",
                    self.alignment.hsp.q_rel_coverage() * 100.0
                ))?;
                td.write_field(&format!("{:.2}", self.alignment.hsp.rel_identity() * 100.0))?;
                td.write_field(&self.alignment.hsp.sseq.len())?; // Alignment length
                td.write_field(&self.ref_accession_frag)?;
                td.write_field(&self.product)?;
                td.write_field(&columns::NA)?; // HMM accession
                td.write_field(&columns::NA)?; // HMM description
                if print_node {
                    td.write_field(&columns::NA)?;
                }
                td.new_ln()?;
            }
        }
        Ok(())
    }
}

impl Batch {
    fn new(mutation_table: &Path) -> Result<Self> {
        let mut accession2mutations: HashMap<String, Vec<AmrMutation>> = HashMap::new();
        let file = File::open(mutation_table)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 7 {
                bail!("Malformed mutation row: {line}");
            }
            let accession = fields[0].to_string();
            let pos: usize = fields[1].parse()?;
            if pos == 0 {
                bail!("Mutation position must be positive: {line}");
            }
            let gene_mutation_std = fields[2];
            let gene_mutation_report = fields[3];
            let class = fields[4];
            let subclass = fields[5];
            let name = fields[6];

            let mutation = AmrMutation::new(
                pos,
                gene_mutation_std,
                gene_mutation_report,
                class,
                subclass,
                name,
            );

            accession2mutations
                .entry(accession)
                .or_default()
                .push(mutation);
        }

        for (accession, mutations) in &mut accession2mutations {
            mutations.sort();
            if mutations
                .windows(2)
                .any(|pair| pair[0].gene_mutation_std == pair[1].gene_mutation_std)
            {
                bail!("Duplicate reference mutations for {accession}");
            }
        }

        Ok(Self {
            blast_als: Vec::new(),
            accession2mutations,
        })
    }

    fn report(
        &self,
        td: &mut TsvOut,
        mutation_all: bool,
        print_node: bool,
        input_name: &str,
    ) -> Result<()> {
        if !input_name.is_empty() {
            td.write_field(&"Name")?;
        }
        td.write_field(&columns::PROT_COL_NAME)?;
        td.write_field(&columns::CONTIG_COL_NAME)?;
        td.write_field(&columns::START_COL_NAME)?;
        td.write_field(&columns::STOP_COL_NAME)?;
        td.write_field(&columns::STRAND_COL_NAME)?;
        td.write_field(&columns::GENESYMBOL_COL_NAME)?;
        td.write_field(&columns::ELEM_NAME_COL_NAME)?;
        td.write_field(&columns::SCOPE_COL_NAME)?;
        td.write_field(&columns::TYPE_COL_NAME)?;
        td.write_field(&columns::SUBTYPE_COL_NAME)?;
        td.write_field(&columns::CLASS_COL_NAME)?;
        td.write_field(&columns::SUBCLASS_COL_NAME)?;
        td.write_field(&columns::METHOD_COL_NAME)?;
        td.write_field(&columns::TARGET_LEN_COL_NAME)?;
        td.write_field(&columns::REF_LEN_COL_NAME)?;
        td.write_field(&columns::REF_COV_COL_NAME)?;
        td.write_field(&columns::REF_IDENT_COL_NAME)?;
        td.write_field(&columns::ALIGN_LEN_COL_NAME)?;
        td.write_field(&columns::CLOSEST_REF_ACCESSION_COL_NAME)?;
        td.write_field(&columns::CLOSEST_REF_NAME_COL_NAME)?;
        td.write_field(&columns::HMM_ACCESSION_COL_NAME)?;
        td.write_field(&columns::HMM_DESCR_COL_NAME)?;
        if print_node {
            td.write_field(&columns::HIERARCHY_NODE_COL_NAME)?;
        }
        td.new_ln()?;

        for al in &self.blast_als {
            al.report(td, mutation_all, print_node, input_name)?;
        }
        Ok(())
    }
}

pub fn body(
    blastn_file: &Path,
    mutation_table: &Path,
    organism: &str,
    mutation_all_file: Option<&Path>,
    input_name: &str,
    print_node: bool,
    out: &mut dyn Write,
) -> Result<()> {
    let mut batch = Batch::new(mutation_table)?;

    let file = File::open(blastn_file)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let al = BlastnAlignment::parse(&line, organism, &batch.accession2mutations)?;
        if al.good() {
            batch.blast_als.push(al);
        }
    }

    for i in 0..batch.blast_als.len() {
        for change_i in 0..batch.blast_als[i].alignment.seq_changes.len() {
            let start_target = batch.blast_als[i].alignment.seq_changes[change_i].start_target;
            let rel_identity = batch.blast_als[i].alignment.hsp.rel_identity();
            for j in 0..batch.blast_als.len() {
                if i == j
                    || batch.blast_als[j].alignment.hsp.sseqid
                        != batch.blast_als[i].alignment.hsp.sseqid
                    || batch.blast_als[j].alignment.hsp.s_int.strand
                        != batch.blast_als[i].alignment.hsp.s_int.strand
                {
                    continue;
                }
                for change_j in 0..batch.blast_als[j].alignment.seq_changes.len() {
                    if start_target
                        != batch.blast_als[j].alignment.seq_changes[change_j].start_target
                    {
                        continue;
                    }
                    let better = batch.blast_als[i].alignment.seq_changes[change_i].better(
                        &batch.blast_als[j].alignment.seq_changes[change_j],
                        rel_identity,
                        batch.blast_als[j].alignment.hsp.rel_identity(),
                    );
                    if better {
                        batch.blast_als[j].alignment.seq_changes[change_j].replacement =
                            Some(change_i);
                    }
                }
            }
        }
    }

    let mut td = TsvOut::new(Some(out));
    td.use_pound = false;
    batch.report(&mut td, false, print_node, input_name)?;

    if let Some(mutation_all_file) = mutation_all_file {
        let mut file = File::create(mutation_all_file)?;
        let mut td = TsvOut::new(Some(&mut file));
        td.use_pound = false;
        batch.report(&mut td, true, print_node, input_name)?;
    }

    Ok(())
}

pub fn main(argv: Vec<OsString>, out: &mut dyn Write) -> Result<i32> {
    let mut app = ThisApplication::new();
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    let blastn_file = PathBuf::from(&run.positional_args[0]);
    let mutation_table = PathBuf::from(&run.positional_args[1]);
    let organism = run.positional_args[2].to_string_lossy();
    let mutation_all = run.key_args["mutation_all"].is_empty().then_some(None);
    let mutation_all =
        mutation_all.unwrap_or_else(|| Some(PathBuf::from(&run.key_args["mutation_all"])));

    body(
        &blastn_file,
        &mutation_table,
        &organism,
        mutation_all.as_deref(),
        &run.key_args["name"],
        run.key_args["print_node"] == "true",
        out,
    )?;
    Ok(0)
}
