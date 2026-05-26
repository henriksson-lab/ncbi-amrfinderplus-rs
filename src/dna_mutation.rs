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
    gene: String,
    seq_changes: Vec<DnaSeqChange>,
}

#[derive(Clone)]
struct DnaSeqChange {
    start: usize,
    len: usize,
    reference: String,
    allele: String,
    start_ref: usize,
    stop_ref: usize,
    start_target: usize,
    neighborhood_mismatch: f64,
    mutations: Vec<AmrMutation>,
    is_wildtype: bool,
    replacement: bool,
}

impl BlastnAlignment {
    fn parse(
        line: &str,
        organism: &str,
        accession2mutations: &HashMap<String, Vec<AmrMutation>>,
    ) -> Result<Self> {
        let hsp = Hsp::from_blast_line(line, false, false, false, false, false)?;

        let organism = organism.replace('_', " ");

        // Parse qseqid: accession@gene_name@gene_symbol@offset:start-stop
        let qseqid = &hsp.qseqid;
        let parts: Vec<&str> = qseqid.splitn(4, '@').collect();
        let (ref_accession_frag, product, gene) = if parts.len() >= 3 {
            let accession = parts[0];
            let product = parts[1].replace('_', " ");
            let gene_and_rest = parts[2];
            let gene = gene_and_rest
                .split(':')
                .next()
                .unwrap_or(gene_and_rest)
                .to_string();

            // Reconstruct ref_accession_frag
            let frag = if parts.len() >= 4 {
                format!("{}:{}", accession, parts[3])
            } else {
                let rest_after_gene = gene_and_rest
                    .find(':')
                    .map(|pos| &gene_and_rest[pos + 1..])
                    .unwrap_or("");
                format!("{}:{}", accession, rest_after_gene)
            };

            (frag, product, gene)
        } else {
            (qseqid.clone(), String::new(), String::new())
        };

        let mut seq_changes = accession2mutations
            .get(qseqid)
            .map(|ref_mutations| build_seq_changes(&hsp, ref_mutations))
            .unwrap_or_default();
        seq_changes.sort_by_key(|change| change.start_target);

        Ok(BlastnAlignment {
            hsp,
            organism,
            ref_accession_frag,
            product,
            gene,
            seq_changes,
        })
    }

    fn good(&self) -> bool {
        let min_len = std::cmp::min(self.hsp.qlen, 2 * FLANKING_LEN + 1);
        self.hsp.sseq.len() >= min_len
    }

    fn report(
        &self,
        td: &mut TsvOut,
        mutation_all: bool,
        print_node: bool,
        name: &str,
    ) -> Result<()> {
        for change in &self.seq_changes {
            let mutations: Vec<Option<&AmrMutation>> = if change.mutations.is_empty() {
                vec![None]
            } else {
                change.mutations.iter().map(Some).collect()
            };
            for mutation in mutations {
                if !mutation_all && (change.is_wildtype || mutation.is_none() || change.replacement)
                {
                    continue;
                }
                self.report_change(td, change, mutation, print_node, name)?;
            }
        }
        Ok(())
    }

    fn report_change(
        &self,
        td: &mut TsvOut,
        change: &DnaSeqChange,
        mutation: Option<&AmrMutation>,
        print_node: bool,
        name: &str,
    ) -> Result<()> {
        if change.is_wildtype && mutation.is_none() {
            return Ok(());
        }

        let gene_symbol = if let Some(mutation) = mutation {
            if change.is_wildtype {
                mutation.wildtype()
            } else {
                mutation.gene_mutation.clone()
            }
        } else {
            format!("{}_{}", self.gene, change.mutation_string())
        };

        let elem_name = if let Some(mutation) = mutation {
            if change.is_wildtype {
                format!("{} {} [WILDTYPE]", self.organism, self.product)
            } else {
                mutation.name.replace('_', " ")
            }
        } else {
            format!("{} {} [UNKNOWN]", self.organism, self.product)
        };

        if !name.is_empty() {
            td.write_field(&name)?;
        }
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
        Ok(())
    }
}

impl DnaSeqChange {
    fn unknown(hsp: &Hsp, mut start: usize, mut len: usize, flanking_len: usize) -> Option<Self> {
        if hsp.qseq.as_bytes().get(start) == Some(&b'-') {
            if start == 0 {
                return None;
            }
            start -= 1;
            len += 1;
        }
        let mut change = DnaSeqChange {
            start,
            len,
            reference: String::new(),
            allele: String::new(),
            start_ref: 0,
            stop_ref: 0,
            start_target: 0,
            neighborhood_mismatch: 0.0,
            mutations: Vec::new(),
            is_wildtype: false,
            replacement: false,
        };
        change.set_seq(hsp);
        if change.reference.is_empty() {
            return None;
        }
        change.set_start_stop_ref(hsp);
        change.set_start_target(hsp);
        change.set_neighborhood_mismatch(hsp, flanking_len);
        (change.neighborhood_mismatch <= 0.04).then_some(change)
    }

    fn wildtype(mutation: AmrMutation) -> Self {
        DnaSeqChange {
            start: 0,
            len: 0,
            reference: String::new(),
            allele: String::new(),
            start_ref: mutation.pos_real,
            stop_ref: mutation.get_stop(),
            start_target: 0,
            neighborhood_mismatch: 0.0,
            mutations: vec![mutation],
            is_wildtype: true,
            replacement: false,
        }
    }

    fn set_seq(&mut self, hsp: &Hsp) {
        self.reference = hsp.qseq[self.start..self.start + self.len]
            .replace('-', "")
            .to_uppercase();
        self.allele = hsp.sseq[self.start..self.start + self.len]
            .replace('-', "")
            .to_uppercase();
    }

    fn set_start_stop_ref(&mut self, hsp: &Hsp) {
        self.start_ref = hsp.q_int.start
            + hsp.qseq.as_bytes()[..self.start]
                .iter()
                .filter(|&&b| b != b'-')
                .count();
        self.stop_ref = self.start_ref
            + hsp.qseq.as_bytes()[self.start..self.start + self.len]
                .iter()
                .filter(|&&b| b != b'-')
                .count();
    }

    fn set_start_target(&mut self, hsp: &Hsp) {
        self.start_target = hsp.s_int.start;
        if hsp.s_int.strand == 1 {
            self.start_target += hsp.sseq.as_bytes()[..self.start]
                .iter()
                .filter(|&&b| b != b'-')
                .count();
        } else {
            self.start_target += hsp.sseq.as_bytes()[self.start + self.len..]
                .iter()
                .filter(|&&b| b != b'-')
                .count();
        }
    }

    fn set_neighborhood_mismatch(&mut self, hsp: &Hsp, flanking_len: usize) {
        self.neighborhood_mismatch = 0.0;
        if flanking_len == 0 {
            return;
        }

        let q_bytes = hsp.qseq.as_bytes();
        let s_bytes = hsp.sseq.as_bytes();
        let mut span = 0usize;
        let mut mismatches = 0usize;

        let mut j_left = self.start;
        if self.start > 0 {
            let mut j = self.start - 1;
            while self.start - j <= flanking_len {
                span += 1;
                if s_bytes[j] != q_bytes[j] {
                    mismatches += 1;
                }
                if j == 0 {
                    break;
                }
                j -= 1;
            }
            j_left = j;
        }
        let covered_left = self.start.saturating_sub(j_left);
        let missed_left = s_tail(hsp, true)
            .min(hsp.q_int.start)
            .min(flanking_len + 1 - covered_left.min(flanking_len + 1));
        span += missed_left;
        mismatches += missed_left;

        let stop = self.start + self.len;
        let mut j = stop + 1;
        while j < s_bytes.len() && j - stop <= flanking_len {
            span += 1;
            if s_bytes[j] != q_bytes[j] {
                mismatches += 1;
            }
            j += 1;
        }
        let covered_right = j.saturating_sub(stop);
        let missed_right = s_tail(hsp, false)
            .min(hsp.qlen - hsp.q_int.stop)
            .min(flanking_len + 1 - covered_right.min(flanking_len + 1));
        span += missed_right;
        mismatches += missed_right;

        if span > 0 {
            self.neighborhood_mismatch = mismatches as f64 / span as f64;
        }
    }

    fn matches_mutation(&self, mutation: &AmrMutation) -> bool {
        if self.len == 0
            || mutation.pos_real < self.start_ref
            || mutation.get_stop() > self.stop_ref
        {
            return false;
        }
        let head = mutation.pos_real - self.start_ref;
        let tail = self.stop_ref - mutation.get_stop();
        if head + mutation.reference.len() + tail != self.reference.len() {
            return false;
        }
        if self.reference[head..head + mutation.reference.len()] != mutation.reference {
            return false;
        }
        if head + mutation.allele.len() + tail != self.allele.len() {
            return false;
        }
        self.allele[head..head + mutation.allele.len()] == mutation.allele
    }

    fn mutation_string(&self) -> String {
        let allele = if self.allele.is_empty() {
            "DEL".to_string()
        } else if self.allele == "*" {
            "STOP".to_string()
        } else {
            self.allele.clone()
        };
        format!("{}{}{}", self.reference, self.start_ref + 1, allele)
    }
}

fn build_seq_changes(hsp: &Hsp, ref_mutations: &[AmrMutation]) -> Vec<DnaSeqChange> {
    let q_bytes = hsp.qseq.as_bytes();
    let s_bytes = hsp.sseq.as_bytes();
    let mut changes = Vec::new();
    let mut in_change = false;
    let mut start = 0usize;
    let mut len = 0usize;

    for i in 0..q_bytes.len() {
        if in_change {
            if s_bytes[i] == q_bytes[i] {
                if let Some(change) = DnaSeqChange::unknown(hsp, start, len, FLANKING_LEN) {
                    changes.push(change);
                }
                in_change = false;
                len = 0;
            } else {
                len += 1;
            }
        } else if s_bytes[i] != q_bytes[i] {
            start = i;
            len = 1;
            in_change = true;
        }
    }
    if in_change {
        if let Some(change) = DnaSeqChange::unknown(hsp, start, len, FLANKING_LEN) {
            changes.push(change);
        }
    }

    for change in &mut changes {
        for mutation in ref_mutations {
            if change.matches_mutation(mutation) {
                change.mutations.push(mutation.clone());
            }
        }
    }

    add_insertion_mutation_fallbacks(hsp, ref_mutations, &mut changes);

    add_wildtype_changes(hsp, ref_mutations, &mut changes);
    changes
}

fn add_insertion_mutation_fallbacks(
    hsp: &Hsp,
    ref_mutations: &[AmrMutation],
    changes: &mut Vec<DnaSeqChange>,
) {
    let subject_ungapped: String = hsp.sseq.chars().filter(|&c| c != '-').collect();
    for mutation in ref_mutations {
        if mutation.allele.len() <= mutation.reference.len()
            || mutation.pos_real < hsp.q_int.start
            || mutation.pos_real >= hsp.q_int.stop
            || changes.iter().any(|change| {
                !change.is_wildtype
                    && change
                        .mutations
                        .iter()
                        .any(|existing| existing.gene_mutation == mutation.gene_mutation)
            })
        {
            continue;
        }
        let ref_pos = mutation.pos_real - hsp.q_int.start;
        let Some(observed) = subject_ungapped.get(ref_pos..ref_pos + mutation.allele.len()) else {
            continue;
        };
        if observed.eq_ignore_ascii_case(&mutation.allele) {
            suppress_insertion_unknown(changes, mutation);
            changes.push(DnaSeqChange {
                start: ref_pos,
                len: mutation.allele.len(),
                reference: mutation.reference.clone(),
                allele: mutation.allele.clone(),
                start_ref: mutation.pos_real,
                stop_ref: mutation.get_stop(),
                start_target: insertion_sort_target(hsp, ref_pos, mutation),
                neighborhood_mismatch: 0.0,
                mutations: vec![mutation.clone()],
                is_wildtype: false,
                replacement: false,
            });
        }
    }
}

fn suppress_insertion_unknown(changes: &mut Vec<DnaSeqChange>, mutation: &AmrMutation) {
    let inserted = mutation
        .allele
        .strip_prefix(&mutation.reference)
        .unwrap_or(&mutation.allele);
    changes.retain(|change| {
        let anchored_insertion_unknown = !change.is_wildtype
            && change.mutations.is_empty()
            && change.start_ref <= mutation.pos_real
            && change.stop_ref <= mutation.get_stop()
            && mutation.pos_real.saturating_sub(change.stop_ref) <= 1
            && !inserted.is_empty()
            && change.allele.contains(inserted);
        !anchored_insertion_unknown
    });
}

fn insertion_sort_target(hsp: &Hsp, ref_pos: usize, mutation: &AmrMutation) -> usize {
    let offset = ref_pos + mutation.allele.len() + mutation.reference.len();
    if hsp.s_int.strand == 1 {
        hsp.s_int.start + offset
    } else {
        hsp.s_int.stop.saturating_sub(offset)
    }
}

fn add_wildtype_changes(hsp: &Hsp, ref_mutations: &[AmrMutation], changes: &mut Vec<DnaSeqChange>) {
    let q_bytes = hsp.qseq.as_bytes();
    let mut ref_pos = hsp.q_int.start;
    let mut i = 0usize;

    for mutation in ref_mutations {
        while ref_pos < mutation.pos_real {
            if i >= q_bytes.len() {
                return;
            }
            i += 1;
            while i < q_bytes.len() && q_bytes[i] == b'-' {
                i += 1;
            }
            if i >= q_bytes.len() {
                return;
            }
            ref_pos += 1;
        }
        if i >= q_bytes.len() {
            return;
        }
        if ref_pos == mutation.pos_real
            && starts_with_ignore_case(&hsp.qseq[i..], &mutation.reference)
            && starts_with_ignore_case(&hsp.sseq[i..], &mutation.reference)
        {
            changes.push(DnaSeqChange::wildtype(mutation.clone()));
        }
    }
}

fn starts_with_ignore_case(s: &str, prefix: &str) -> bool {
    s.get(..prefix.len())
        .is_some_and(|head| head.eq_ignore_ascii_case(prefix))
}

fn s_tail(hsp: &Hsp, left: bool) -> usize {
    if hsp.s_int.empty() {
        return 0;
    }
    if (left && hsp.s_int.strand == 1) || (!left && hsp.s_int.strand == -1) {
        hsp.s_int.start
    } else {
        hsp.slen.saturating_sub(hsp.s_int.stop)
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
            if pos == 0 {
                continue;
            }
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

            accession2mutations
                .entry(accession)
                .or_default()
                .push(mutation);
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
            let al = BlastnAlignment::parse(&line, organism, &accession2mutations)?;
            if al.good() {
                alignments.push(al);
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
            al.report(&mut td, false, print_node, name)?;
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
            al.report(&mut td, true, print_node, name)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::process::Command;

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
        )
        .unwrap();

        let rust_output = String::from_utf8(output).unwrap();
        let cpp_output = std::fs::read_to_string(&expected_file).unwrap();

        assert_eq!(
            rust_output, cpp_output,
            "DNA mutation output must be byte-identical to the C++ golden output"
        );
    }

    #[test]
    fn test_dna_mutation_all_matches_cpp_byte_for_byte() {
        let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let cpp = manifest.join("amr/dna_mutation");
        let blastn_file = manifest.join("tests/golden/blastn");
        let mutation_table = manifest.join("amrfinder_db/2026-03-24.1/AMR_DNA-Escherichia.tsv");

        if !cpp.exists() || !blastn_file.exists() || !mutation_table.exists() {
            return;
        }

        let cpp_mut_all = tempfile::NamedTempFile::new().unwrap();
        let cpp_out = Command::new(cpp)
            .arg(&blastn_file)
            .arg(&mutation_table)
            .arg("Escherichia")
            .arg("-mutation_all")
            .arg(cpp_mut_all.path())
            .arg("-print_node")
            .output()
            .expect("C++ dna_mutation failed");
        assert!(
            cpp_out.status.success(),
            "C++ dna_mutation failed: {}",
            String::from_utf8_lossy(&cpp_out.stderr)
        );

        let mut rust_output = Vec::new();
        let mut rust_mut_all = Vec::new();
        run_dna_mutation(
            &blastn_file,
            &mutation_table,
            "Escherichia",
            true,
            "",
            &mut rust_output,
            Some(&mut rust_mut_all),
        )
        .unwrap();

        let cpp_mut_all = std::fs::read_to_string(cpp_mut_all.path()).unwrap();
        let rust_mut_all = String::from_utf8(rust_mut_all).unwrap();

        assert_eq!(
            rust_mut_all, cpp_mut_all,
            "DNA mutation_all output must be byte-identical to C++"
        );
    }

    #[test]
    fn ampc_insertion_replaces_generic_unknown_in_mutation_all() {
        let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let blastn_file = manifest.join("tests/golden/blastn");
        let mutation_table = manifest.join("amrfinder_db/2026-03-24.1/AMR_DNA-Escherichia.tsv");

        if !blastn_file.exists() || !mutation_table.exists() {
            return;
        }

        let blastn = std::fs::read_to_string(&blastn_file).unwrap();
        let amp_c_line = blastn
            .lines()
            .find(|line| line.contains("\tcontig17\t"))
            .expect("golden BLASTN should contain the ampC insertion case");
        let focused_blastn = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(focused_blastn.path(), format!("{amp_c_line}\n")).unwrap();

        let mut report = Vec::new();
        let mut mutation_all = Vec::new();
        run_dna_mutation(
            focused_blastn.path(),
            &mutation_table,
            "Escherichia",
            true,
            "",
            &mut report,
            Some(&mut mutation_all),
        )
        .unwrap();

        let mutation_all = String::from_utf8(mutation_all).unwrap();
        assert!(
            mutation_all
                .contains("\tampC_T-14TGT\tEscherichia cephalosporin resistant ampC promoter\t"),
            "known ampC promoter insertion should be reported"
        );
        assert!(
            !mutation_all.contains("\tampC_C38CGT\t"),
            "known insertion should suppress the generic anchored unknown insertion artifact"
        );

        let lines: Vec<&str> = mutation_all.lines().collect();
        let c11 = lines
            .iter()
            .position(|line| line.contains("\tampC_C-11C\t"))
            .expect("wildtype ampC_C-11C row should be present");
        let insertion = lines
            .iter()
            .position(|line| line.contains("\tampC_T-14TGT\t"))
            .expect("ampC insertion row should be present");
        assert!(
            c11 < insertion,
            "C++ sorts empty wildtype SeqChanges before the real insertion change"
        );
    }

    #[test]
    fn reverse_reference_blastn_hit_reports_wildtype_mutation_all_rows() {
        let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let input = manifest.join("tests/golden/test_dna.fa");
        let db = manifest.join("amrfinder_db/2026-03-24.1");
        let dna_db = db.join("AMR_DNA-Escherichia.fa");
        let mutation_table = db.join("AMR_DNA-Escherichia.tsv");

        if !input.exists() || !dna_db.exists() || !mutation_table.exists() {
            return;
        }

        let blastn = tempfile::NamedTempFile::new().unwrap();
        let blastn_status = Command::new("blastn")
            .args([
                "-query",
                input.to_str().unwrap(),
                "-db",
                dna_db.to_str().unwrap(),
                "-evalue",
                "1e-20",
                "-dust",
                "no",
                "-max_target_seqs",
                "10000",
                "-outfmt",
                "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
                "-out",
                blastn.path().to_str().unwrap(),
            ])
            .status();
        let Ok(blastn_status) = blastn_status else {
            return;
        };
        assert!(
            blastn_status.success(),
            "blastn failed for reverse-hit regression"
        );

        let blastn_output = std::fs::read_to_string(blastn.path()).unwrap();
        let reverse_23s_line = blastn_output
            .lines()
            .find(|line| line.contains("\tcontig05\t"))
            .expect("test DNA should generate a reverse-orientation 23S BLASTN hit");
        let focused_blastn = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(focused_blastn.path(), format!("{reverse_23s_line}\n")).unwrap();

        let mut report = Vec::new();
        let mut mutation_all = Vec::new();
        run_dna_mutation(
            focused_blastn.path(),
            &mutation_table,
            "Escherichia",
            true,
            "",
            &mut report,
            Some(&mut mutation_all),
        )
        .unwrap();

        let mutation_all = String::from_utf8(mutation_all).unwrap();
        for symbol in [
            "23S_A2058A",
            "23S_C2611C",
            "23S_G2057G",
            "23S_G2447G",
            "23S_T2609T",
        ] {
            assert!(
                mutation_all.contains(&format!("NA\tcontig05\t237\t1224\t-\t{symbol}\t")),
                "reverse reference hit should report {symbol} as a contig05 wildtype row"
            );
        }
    }
}
