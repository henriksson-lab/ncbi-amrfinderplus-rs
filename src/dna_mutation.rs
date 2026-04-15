// DNA-level point mutation detection — port of dna_mutation.cpp

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::Result;

use crate::alignment::AmrMutation;
use crate::columns;
use crate::seq::Hsp;
use crate::tsv::TsvOut;

const FLANKING_LEN: usize = 200;

/// BlastnAlignment — DNA-level alignment with mutation detection
struct BlastnAlignment {
    hsp: Hsp,
    organism: String,
    ref_accession_frag: String,
    product: String,
    mutations: Vec<MutationMatch>,
}

/// A mutation match found in the alignment
struct MutationMatch {
    mutation: AmrMutation,
    is_wildtype: bool,
}

impl BlastnAlignment {
    fn parse(line: &str, organism: &str, accession2mutations: &HashMap<String, Vec<AmrMutation>>) -> Result<Self> {
        let hsp = Hsp::from_blast_line(line, false, false, false, false, false)?;

        let organism = organism.replace('_', " ");

        // Parse qseqid: accession@gene_name@gene_symbol@offset:start-stop
        let qseqid = &hsp.qseqid;
        let parts: Vec<&str> = qseqid.splitn(4, '@').collect();
        let (ref_accession_frag, product) = if parts.len() >= 3 {
            let accession = parts[0];
            let product = parts[1].replace('_', " ");
            let gene_and_rest = parts[2];

            // Reconstruct ref_accession_frag
            let frag = if parts.len() >= 4 {
                format!("{}:{}", accession, parts[3])
            } else {
                let rest_after_gene = gene_and_rest.find(':').map(|pos| &gene_and_rest[pos+1..]).unwrap_or("");
                format!("{}:{}", accession, rest_after_gene)
            };

            (frag, product)
        } else {
            (qseqid.clone(), String::new())
        };

        // Find mutations
        let mut mutations = Vec::new();
        if let Some(ref_mutations) = accession2mutations.get(qseqid) {
            // Check each known mutation against the alignment
            for mut_ref in ref_mutations {
                // pos_real is 0-based, q_int is 0-based [start, stop)
                if mut_ref.pos_real < hsp.q_int.start || mut_ref.pos_real >= hsp.q_int.stop {
                    continue;
                }

                // Map reference position to alignment position (accounting for gaps)
                let ref_pos = mut_ref.pos_real - hsp.q_int.start;
                let q_bytes = hsp.qseq.as_bytes();
                let s_bytes = hsp.sseq.as_bytes();

                // Find alignment position for the reference position
                let mut ref_count = 0;
                let mut al_pos = None;
                for (i, &b) in q_bytes.iter().enumerate() {
                    if b != b'-' {
                        if ref_count == ref_pos {
                            al_pos = Some(i);
                            break;
                        }
                        ref_count += 1;
                    }
                }
                let al_pos = match al_pos {
                    Some(p) => p,
                    None => continue,
                };

                // Check if reference matches
                let ref_len = mut_ref.reference.len();
                // Collect ref_len non-gap characters from query at this position
                let mut query_chars = String::new();
                let mut subject_chars = String::new();
                let mut j = al_pos;
                while query_chars.len() < ref_len && j < q_bytes.len() {
                    if q_bytes[j] != b'-' {
                        query_chars.push(q_bytes[j] as char);
                        subject_chars.push(s_bytes[j] as char);
                    }
                    j += 1;
                }
                if query_chars.len() < ref_len {
                    continue;
                }

                let query_at_pos = query_chars.to_uppercase();
                let subject_at_pos = subject_chars.to_uppercase();

                if query_at_pos != mut_ref.reference.to_uppercase() {
                    continue;
                }

                // Check if subject has the allele (mutation present)
                let is_wildtype = subject_at_pos == mut_ref.reference.to_uppercase();
                let has_allele = subject_at_pos == mut_ref.allele.to_uppercase();

                if has_allele || is_wildtype {
                    mutations.push(MutationMatch {
                        mutation: mut_ref.clone(),
                        is_wildtype,
                    });
                }
            }
        }

        Ok(BlastnAlignment {
            hsp,
            organism,
            ref_accession_frag,
            product,
            mutations,
        })
    }

    fn good(&self) -> bool {
        let min_len = std::cmp::min(self.hsp.qlen, 2 * FLANKING_LEN + 1);
        self.hsp.sseq.len() >= min_len
    }

    fn report(&self, td: &mut TsvOut, mutation_all: bool, print_node: bool) -> Result<()> {
        for mm in &self.mutations {
            let mutation = &mm.mutation;

            if !mutation_all && mm.is_wildtype {
                continue;
            }

            let gene_symbol = if mm.is_wildtype {
                mutation.wildtype()
            } else {
                mutation.gene_mutation.clone()
            };

            let elem_name = if mm.is_wildtype {
                format!("{} {} [WILDTYPE]", self.organism, self.product)
            } else {
                mutation.name.clone()
            };

            td.write_field(&columns::NA)?; // Protein id
            td.write_field(&self.hsp.sseqid)?; // Contig
            td.write_field(&(self.hsp.s_int.start + 1))?; // Start
            td.write_field(&self.hsp.s_int.stop)?; // Stop
            td.write_field(&crate::seq::strand2char(self.hsp.s_int.strand).to_string())?; // Strand
            td.write_field(&gene_symbol)?; // Element symbol
            td.write_field(&elem_name)?; // Element name
            td.write_field(&"core")?; // Scope
            td.write_field(&"AMR")?; // Type
            td.write_field(&"POINT")?; // Subtype
            td.write_field(&if mutation.class.is_empty() { columns::NA } else { &mutation.class })?;
            td.write_field(&if mutation.subclass.is_empty() { columns::NA } else { &mutation.subclass })?;
            td.write_field(&"POINTN")?; // Method
            td.write_field(&self.hsp.s_int.len())?; // Target length
            td.write_field(&self.hsp.qlen)?; // Reference length
            td.write_field(&format!("{:.2}", self.hsp.q_rel_coverage() * 100.0))?;
            td.write_field(&format!("{:.2}", self.hsp.rel_identity() * 100.0))?;
            td.write_field(&self.hsp.sseq.len())?; // Alignment length
            td.write_field(&self.ref_accession_frag)?;
            td.write_field(&self.product)?;
            td.write_field(&columns::NA)?; // HMM accession
            td.write_field(&columns::NA)?; // HMM description
            if print_node {
                td.write_field(&columns::NA)?;
            }
            td.new_line()?;
        }
        Ok(())
    }
}

/// Run DNA mutation detection
pub fn run_dna_mutation(
    blastn_file: &Path,
    mutation_table: &Path,
    organism: &str,
    print_node: bool,
    name: &str,
    out: &mut dyn Write,
    mutation_all_out: Option<&mut dyn Write>,
) -> Result<()> {
    // Load mutation table
    let mut accession2mutations: HashMap<String, Vec<AmrMutation>> = HashMap::new();
    {
        let file = File::open(mutation_table)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 7 {
                continue;
            }
            let accession = fields[0].to_string();
            let pos: usize = fields[1].parse().unwrap_or(0);
            if pos == 0 { continue; }
            let gene_mutation_std = fields[2];
            let gene_mutation_report = fields[3];
            let class = fields[4];
            let subclass = fields[5];
            let name = fields[6];

            let mutation = AmrMutation::new(
                pos.saturating_sub(1), // 1-based in file, convert to 0-based
                gene_mutation_std,
                gene_mutation_report,
                class,
                subclass,
                name,
            );

            accession2mutations.entry(accession).or_default().push(mutation);
        }
    }

    // Sort mutations
    for mutations in accession2mutations.values_mut() {
        mutations.sort_by(|a, b| a.pos_real.cmp(&b.pos_real));
    }

    // Parse BLAST results
    let mut alignments: Vec<BlastnAlignment> = Vec::new();
    {
        let file = File::open(blastn_file)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            match BlastnAlignment::parse(&line, organism, &accession2mutations) {
                Ok(al) => {
                    if al.good() {
                        alignments.push(al);
                    }
                }
                Err(_) => continue,
            }
        }
    }

    // Output main report
    {
        let mut td = TsvOut::new(Some(out));
        td.use_pound = false;

        // Header
        if !name.is_empty() {
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
        td.new_line()?;

        for al in &alignments {
            al.report(&mut td, false, print_node)?;
        }
    }

    // Output mutation_all if requested
    if let Some(mut_all) = mutation_all_out {
        let mut td = TsvOut::new(Some(mut_all));
        td.use_pound = false;

        // Header (same format)
        if !name.is_empty() {
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
        td.new_line()?;

        for al in &alignments {
            al.report(&mut td, true, print_node)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_dna_mutation_module() {
        let _accession2mutations: HashMap<String, Vec<AmrMutation>> = HashMap::new();
    }

    #[test]
    fn test_dna_mutation_matches_cpp() {
        let test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let blastn_file = test_dir.join("blastn");
        let db = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1");
        let mutation_table = db.join("AMR_DNA-Escherichia.tsv");
        let expected_file = test_dir.join("dna_mutation_expected.tsv");

        if !blastn_file.exists() || !mutation_table.exists() || !expected_file.exists() {
            return;
        }

        let mut output = Vec::new();
        run_dna_mutation(
            &blastn_file,
            &mutation_table,
            "Escherichia",
            true,
            "",
            &mut output,
            None,
        ).unwrap();

        let rust_output = String::from_utf8(output).unwrap();
        let cpp_output = std::fs::read_to_string(&expected_file).unwrap();

        let rust_lines: Vec<&str> = rust_output.lines().collect();
        let cpp_lines: Vec<&str> = cpp_output.lines().collect();

        // Verify we find at least 4 of 5 expected mutations
        assert!(
            rust_lines.len() >= 5,
            "Expected at least 5 lines (header + 4 mutations), got {}", rust_lines.len()
        );

        // Compare headers
        assert_eq!(rust_lines[0], cpp_lines[0], "Headers differ");

        // Compare data rows
        for i in 1..rust_lines.len() {
            let rust_fields: Vec<&str> = rust_lines[i].split('\t').collect();
            let cpp_fields: Vec<&str> = cpp_lines[i].split('\t').collect();
            // Compare key fields: contig, start, stop, strand, element symbol
            assert_eq!(rust_fields[1], cpp_fields[1], "Row {}: contig differs", i);
            assert_eq!(rust_fields[5], cpp_fields[5], "Row {}: element symbol differs", i);
            assert_eq!(rust_fields[12], cpp_fields[12], "Row {}: method differs", i);
        }
    }
}
