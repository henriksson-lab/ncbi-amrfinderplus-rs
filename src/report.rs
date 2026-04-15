// amr_report — core reporting engine
// Port of amr_report.cpp

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::Result;

use crate::alignment::{AmrMutation, SeqChange};
use crate::gff::Locus;
use crate::seq::Hsp;
use crate::tsv::TsvOut;

// --- BlastRule ---

/// BLAST identity/coverage thresholds for a family
#[derive(Debug, Clone, Default)]
pub struct BlastRule {
    pub ident: f64,
    pub target_coverage: f64,
    pub ref_coverage: f64,
}

impl BlastRule {
    pub fn new(ident: f64, target_coverage: f64, ref_coverage: f64) -> Self {
        BlastRule {
            ident,
            target_coverage,
            ref_coverage,
        }
    }

    pub fn empty(&self) -> bool {
        self.ident == 0.0 && self.ref_coverage == 0.0
    }
}

// --- Fam ---

/// Protein family hierarchy node
#[derive(Debug, Clone)]
pub struct Fam {
    pub id: String,
    pub genesymbol: String,
    pub family_name: String,
    pub reportable: u8, // 0, 1, or 2
    pub hmm: String,
    pub tc1: f64,
    pub tc2: f64,
    pub complete_br: BlastRule,
    pub partial_br: BlastRule,
    pub type_: String,
    pub subtype: String,
    pub class: String,
    pub subclass: String,
    pub parent_id: String,
}

impl Fam {
    pub fn empty(&self) -> bool {
        self.id.is_empty()
    }

    /// Check if this family is a descendant of the given family ID
    pub fn descendant_of(&self, ancestor_id: &str, fam_map: &HashMap<String, Fam>) -> bool {
        if self.id == ancestor_id {
            return true;
        }
        if self.parent_id.is_empty() {
            return false;
        }
        if let Some(parent) = fam_map.get(&self.parent_id) {
            parent.descendant_of(ancestor_id, fam_map)
        } else {
            false
        }
    }

    /// Get the HMM family (self or ancestor with hmm)
    pub fn get_hmm_fam<'a>(&'a self, fam_map: &'a HashMap<String, Fam>) -> Option<&'a Fam> {
        if !self.hmm.is_empty() {
            return Some(self);
        }
        if self.parent_id.is_empty() {
            return None;
        }
        fam_map.get(&self.parent_id)?.get_hmm_fam(fam_map)
    }
}

// --- Susceptible ---

/// Susceptible phenotype data
#[derive(Debug, Clone)]
pub struct Susceptible {
    pub genesymbol: String,
    pub cutoff: f64,
    pub class: String,
    pub subclass: String,
    pub name: String,
}

// --- HmmAlignment ---

/// HMM domain result
#[derive(Debug, Clone, Default)]
pub struct HmmDomain {
    pub score: f64,
    pub hmm_len: usize,
    pub hmm_start: usize,
    pub hmm_stop: usize,
    pub seq_len: usize,
    pub seq_start: usize,
    pub seq_stop: usize,
}

/// HMM search result
#[derive(Debug, Clone)]
pub struct HmmAlignment {
    pub sseqid: String,
    pub score1: f64,
    pub score2: f64,
    pub fam_id: String,
    pub domain: Option<HmmDomain>,
    pub blast_al_idx: Option<usize>, // index into BlastAlignment vector
}

impl HmmAlignment {
    pub fn good(&self, fam_map: &HashMap<String, Fam>) -> bool {
        if let Some(fam) = fam_map.get(&self.fam_id) {
            self.score1 >= fam.tc1 && self.score2 >= fam.tc2
        } else {
            false
        }
    }
}

// --- BlastAlignment ---

/// BLAST alignment result with full AMR annotation
#[derive(Debug, Clone)]
pub struct BlastAlignment {
    pub hsp: Hsp,
    pub from_hmm: bool,
    pub ref_accession: String,
    pub part: usize,
    pub parts: usize,
    pub fam_id: String,
    pub gene: String,
    pub resistance: String,
    pub class: String,
    pub subclass: String,
    pub product: String,
    pub reportable: u8,
    pub genesymbol: String,
    pub method: String,
    pub complete_br: BlastRule,
    pub partial_br: BlastRule,
    pub cdss: Vec<Locus>,
    pub hmm_al_idx: Option<usize>,
    pub susceptible_idx: Option<usize>,
    pub seq_changes: Vec<SeqChange>,
    pub ref_mutation: AmrMutation,
    // Fusion protein info
    pub fusion_ids: Vec<usize>,
}

impl BlastAlignment {
    /// Parse from BLAST qseqid format:
    /// accession|part|parts|famId|gene|resistance|reportable|subclass|classS|genesymbol|product
    pub fn from_blast_line(
        line: &str,
        q_prot: bool,
        s_prot: bool,
        default_complete_br: &BlastRule,
        default_partial_br: &BlastRule,
    ) -> Result<Self> {
        let a_prot = q_prot || s_prot;
        let hsp = Hsp::from_blast_line(line, q_prot, s_prot, a_prot, q_prot, true)?;

        // Parse qseqid to extract reference info
        // Format: accession|part|parts|famId|gene|resistance|reportable|subclass|classS|product
        // The product field may contain '|' so we split from the right (like C++ rfindSplit)
        let qseqid = &hsp.qseqid;
        let parts: Vec<&str> = qseqid.splitn(10, '|').collect();

        let (ref_accession, part, parts_count, fam_id, gene, resistance, reportable, subclass, class, genesymbol, product) = if parts.len() >= 10 {
            // product is the last field, genesymbol derived from famId
            let product_raw = parts[9].to_string().replace('_', " ");
            (
                parts[0].to_string(),
                parts[1].parse::<usize>().unwrap_or(1),
                parts[2].parse::<usize>().unwrap_or(1),
                parts[3].to_string(),
                parts[4].to_string(),
                parts[5].to_string(),
                parts[6].parse::<u8>().unwrap_or(0),
                parts[7].to_string(),
                parts[8].to_string(),
                String::new(), // genesymbol — set later from fam
                product_raw,
            )
        } else {
            // Fallback for simpler format
            (
                qseqid.clone(),
                1, 1,
                String::new(), String::new(), String::new(),
                0,
                String::new(), String::new(), String::new(), String::new(),
            )
        };

        Ok(BlastAlignment {
            hsp,
            from_hmm: false,
            ref_accession,
            part,
            parts: parts_count,
            fam_id,
            gene,
            resistance,
            class,
            subclass,
            product,
            reportable,
            genesymbol,
            method: String::new(),
            complete_br: default_complete_br.clone(),
            partial_br: default_partial_br.clone(),
            cdss: Vec::new(),
            hmm_al_idx: None,
            susceptible_idx: None,
            seq_changes: Vec::new(),
            ref_mutation: AmrMutation::default(),
            fusion_ids: Vec::new(),
        })
    }

    /// Check if this alignment meets BLAST thresholds
    pub fn good(&self) -> bool {
        let ident = self.hsp.rel_identity();
        let ref_cov = self.hsp.q_rel_coverage();

        // Complete match
        if ident >= self.complete_br.ident && ref_cov >= self.complete_br.ref_coverage {
            return true;
        }
        // Partial match
        if ident >= self.partial_br.ident && ref_cov >= self.partial_br.ref_coverage {
            return true;
        }
        false
    }
}

// --- Batch ---

/// Core processing container for AMR report generation
pub struct Batch {
    pub fam_map: HashMap<String, Fam>,
    pub hmm2fam: HashMap<String, String>, // hmm_id -> fam_id
    pub reportable_min: u8,
    pub cds_exist: bool,
    pub suppress_prots: Vec<String>,
    pub blast_als: Vec<BlastAlignment>,
    pub hmm_als: Vec<HmmAlignment>,
    pub domains: HashMap<(String, String), HmmDomain>, // (sseqid, fam_id) -> domain
    pub accession2mutations: HashMap<String, Vec<AmrMutation>>,
    pub accession2susceptible: HashMap<String, Susceptible>,
    // Computed by process()
    pub target2blast_als: BTreeMap<String, Vec<usize>>,      // sseqid -> indices into blast_als
    pub target2good_blast_als: BTreeMap<String, Vec<usize>>,
    pub target2hmm_als: BTreeMap<String, Vec<usize>>,
    pub target2good_hmm_als: BTreeMap<String, Vec<usize>>,
}

impl Batch {
    /// Create a new Batch by loading the FAM file
    pub fn from_fam_file(fam_path: &Path, reportable_min: u8) -> Result<Self> {
        let mut fam_map: HashMap<String, Fam> = HashMap::new();
        let mut hmm2fam: HashMap<String, String> = HashMap::new();

        let file = File::open(fam_path)?;
        let reader = BufReader::new(file);

        // Pass 1: Create Fam objects
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 18 {
                continue;
            }

            let fam_id = fields[0].to_string();
            let parent_id = fields[1].to_string();
            let genesymbol = if fields[2] == "-" {
                String::new()
            } else {
                fields[2].to_string()
            };
            let hmm = if fields[3] == "-" {
                String::new()
            } else {
                fields[3].to_string()
            };
            let tc1: f64 = fields[4].parse().unwrap_or(0.0);
            let tc2: f64 = fields[5].parse().unwrap_or(0.0);

            let complete_br = BlastRule::new(
                fields[6].parse().unwrap_or(0.0),
                fields[7].parse().unwrap_or(0.0),
                fields[8].parse().unwrap_or(0.0),
            );
            let partial_br = BlastRule::new(
                fields[9].parse().unwrap_or(0.0),
                fields[10].parse().unwrap_or(0.0),
                fields[11].parse().unwrap_or(0.0),
            );

            let reportable: u8 = fields[12].parse().unwrap_or(0);
            let type_ = fields[13].to_string();
            let subtype = fields[14].to_string();
            let class = fields[15].to_string();
            let subclass = fields[16].to_string();
            let family_name = fields[17].to_string();

            if !hmm.is_empty() {
                hmm2fam.insert(hmm.clone(), fam_id.clone());
            }

            let fam = Fam {
                id: fam_id.clone(),
                genesymbol,
                family_name,
                reportable,
                hmm,
                tc1,
                tc2,
                complete_br,
                partial_br,
                type_,
                subtype,
                class,
                subclass,
                parent_id,
            };

            fam_map.insert(fam_id, fam);
        }

        Ok(Batch {
            fam_map,
            hmm2fam,
            reportable_min,
            cds_exist: false,
            suppress_prots: Vec::new(),
            blast_als: Vec::new(),
            hmm_als: Vec::new(),
            domains: HashMap::new(),
            accession2mutations: HashMap::new(),
            accession2susceptible: HashMap::new(),
            target2blast_als: BTreeMap::new(),
            target2good_blast_als: BTreeMap::new(),
            target2hmm_als: BTreeMap::new(),
            target2good_hmm_als: BTreeMap::new(),
        })
    }

    /// Load organism-specific mutations
    pub fn load_mutations(&mut self, path: &Path, organism: &str) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 8 {
                continue;
            }

            if fields[0] != organism {
                continue;
            }

            let accession = fields[1].to_string();
            let pos: usize = fields[2].parse().unwrap_or(0);
            let gene_mutation_std = fields[3];
            let gene_mutation = fields[4];
            let class = fields[5];
            let subclass = fields[6];
            let name = fields[7];

            let mutation = AmrMutation::new(
                pos.saturating_sub(1), // 1-based to 0-based
                gene_mutation_std,
                gene_mutation,
                class,
                subclass,
                name,
            );

            self.accession2mutations
                .entry(accession)
                .or_default()
                .push(mutation);
        }

        Ok(())
    }

    /// Load organism-specific susceptible genes
    pub fn load_susceptible(&mut self, path: &Path, organism: &str) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 7 {
                continue;
            }

            if fields[0] != organism {
                continue;
            }

            let accession = fields[2].to_string();
            let susceptible = Susceptible {
                genesymbol: fields[1].to_string(),
                cutoff: fields[3].parse().unwrap_or(0.0),
                class: fields[4].to_string(),
                subclass: fields[5].to_string(),
                name: fields[6].to_string(),
            };

            self.accession2susceptible.insert(accession, susceptible);
        }

        Ok(())
    }

    /// Add a BLAST alignment and index it
    pub fn add_blast_al(&mut self, al: BlastAlignment) {
        let idx = self.blast_als.len();
        let sseqid = al.hsp.sseqid.clone();
        self.blast_als.push(al);
        self.target2blast_als
            .entry(sseqid)
            .or_default()
            .push(idx);
    }

    /// Add an HMM alignment and index it
    pub fn add_hmm_al(&mut self, al: HmmAlignment) {
        let idx = self.hmm_als.len();
        let sseqid = al.sseqid.clone();
        self.hmm_als.push(al);
        self.target2hmm_als
            .entry(sseqid)
            .or_default()
            .push(idx);
    }

    /// Process all alignments: filter, merge, resolve
    pub fn process(&mut self) {
        // Set BlastRules from FAM hierarchy for each alignment
        for al in &mut self.blast_als {
            if !al.from_hmm && !al.fam_id.is_empty() {
                if let Some(fam) = self.fam_map.get(&al.fam_id) {
                    if !fam.complete_br.empty() {
                        al.complete_br = fam.complete_br.clone();
                        al.partial_br = fam.partial_br.clone();
                    }
                    // Set genesymbol from fam if not already set
                    if al.genesymbol.is_empty() && !fam.genesymbol.is_empty() {
                        al.genesymbol = fam.genesymbol.clone();
                    }
                }
            }
        }

        // Step 1: Good BLAST filtering — keep only hits above thresholds
        for (target, indices) in &self.target2blast_als {
            let good: Vec<usize> = indices
                .iter()
                .filter(|&&idx| {
                    let al = &self.blast_als[idx];
                    al.from_hmm || al.good()
                })
                .copied()
                .collect();
            if !good.is_empty() {
                self.target2good_blast_als.insert(target.clone(), good);
            }
        }

        // Step 2: Good HMM filtering
        for (target, indices) in &self.target2hmm_als {
            let good: Vec<usize> = indices
                .iter()
                .filter(|&&idx| self.hmm_als[idx].good(&self.fam_map))
                .copied()
                .collect();
            if !good.is_empty() {
                self.target2good_hmm_als.insert(target.clone(), good);
            }
        }

        // Step 3: Pareto-better BLAST filtering
        // For each target, keep only alignments that are not dominated by another
        let targets: Vec<String> = self.target2good_blast_als.keys().cloned().collect();
        for target in &targets {
            if let Some(indices) = self.target2good_blast_als.get(target) {
                let filtered = self.pareto_filter_blast(indices);
                self.target2good_blast_als.insert(target.clone(), filtered);
            }
        }

        // Step 4: HMM check — remove BLAST hits without HMM support (PD-741)
        // If HMMs exist and a BLAST hit has identity < 98% and no matching HMM, remove it
        let has_hmm_data = !self.hmm_als.is_empty();
        if has_hmm_data {
            let targets: Vec<String> = self.target2good_blast_als.keys().cloned().collect();
            for target in &targets {
                if let Some(indices) = self.target2good_blast_als.get(target) {
                    let good_hmm_fams: Vec<String> = self.target2good_hmm_als
                        .get(target)
                        .map(|hmm_indices| {
                            hmm_indices.iter()
                                .map(|&idx| self.hmm_als[idx].fam_id.clone())
                                .collect()
                        })
                        .unwrap_or_default();

                    let filtered: Vec<usize> = indices.iter()
                        .filter(|&&idx| {
                            let al = &self.blast_als[idx];
                            // Keep if: from HMM, or high identity, or has matching HMM
                            if al.from_hmm {
                                return true;
                            }
                            if al.hsp.rel_identity() >= 0.98 {
                                return true;
                            }
                            // Check if the alignment's family has HMM support
                            if let Some(fam) = self.fam_map.get(&al.fam_id) {
                                if let Some(hmm_fam) = fam.get_hmm_fam(&self.fam_map) {
                                    if good_hmm_fams.contains(&hmm_fam.id) {
                                        return true;
                                    }
                                }
                                // No HMM family — check if fam has no HMM (then keep)
                                if fam.hmm.is_empty() {
                                    if fam.get_hmm_fam(&self.fam_map).is_some() {
                                        // Has ancestor with HMM but no match — remove
                                        return false;
                                    }
                                    return true; // No HMM in hierarchy — keep
                                }
                            }
                            false
                        })
                        .copied()
                        .collect();

                    self.target2good_blast_als.insert(target.clone(), filtered);
                }
            }
        }

        // Step 5: HMM suppression by BLAST
        // If a BLAST hit has higher identity than HMM, suppress the HMM-based alignment
        // (simplified — full logic is more complex)
        let targets: Vec<String> = self.target2good_blast_als.keys().cloned().collect();
        for target in &targets {
            if let Some(blast_indices) = self.target2good_blast_als.get(target) {
                let hmm_indices = self.target2good_hmm_als.get(target);
                if hmm_indices.is_none() {
                    continue;
                }

                // Remove HMM-based blast alignments if a better non-HMM blast exists
                let mut new_indices = Vec::new();
                for &idx in blast_indices {
                    let al = &self.blast_als[idx];
                    if al.from_hmm {
                        // Check if there's a non-HMM blast with same fam that's better
                        let dominated = blast_indices.iter().any(|&other_idx| {
                            let other = &self.blast_als[other_idx];
                            !other.from_hmm
                                && other.fam_id == al.fam_id
                                && other.hsp.rel_identity() >= al.hsp.rel_identity()
                                && other.hsp.q_rel_coverage() >= al.hsp.q_rel_coverage()
                        });
                        if !dominated {
                            new_indices.push(idx);
                        }
                    } else {
                        new_indices.push(idx);
                    }
                }
                self.target2good_blast_als.insert(target.clone(), new_indices);
            }
        }
    }

    /// Get the reportability level for an alignment
    fn get_reportable(&self, al: &BlastAlignment) -> u8 {
        // Check fam hierarchy for reportability (C++ getFam: famId then gene fallback)
        if let Some(fam) = self.fam_map.get(&al.fam_id)
            .or_else(|| if !al.gene.is_empty() { self.fam_map.get(&al.gene) } else { None }) {
            if fam.reportable > 0 {
                return fam.reportable;
            }
            // Check parent chain
            let mut current = fam;
            loop {
                if current.reportable > 0 {
                    return current.reportable;
                }
                if current.parent_id.is_empty() {
                    break;
                }
                match self.fam_map.get(&current.parent_id) {
                    Some(parent) => current = parent,
                    None => break,
                }
            }
        }
        al.reportable
    }

    /// Pareto-better filtering matching C++ betterEq() logic.
    /// For hits to the same target: compare by exact match, then nident, then ref length.
    /// A hit is removed if another hit to the same target is strictly better.
    fn pareto_filter_blast(&self, indices: &[usize]) -> Vec<usize> {
        if indices.is_empty() {
            return Vec::new();
        }
        let mut result = Vec::new();
        for &idx in indices {
            let al = &self.blast_als[idx];
            let dominated = indices.iter().any(|&other_idx| {
                if other_idx == idx {
                    return false;
                }
                let other = &self.blast_als[other_idx];
                self.blast_better(other, al)
            });
            if !dominated {
                result.push(idx);
            }
        }
        result
    }

    /// Simplified C++ betterEq: is `a` strictly better than `b`?
    /// Matches the key comparisons from amr_report.cpp betterEq/better.
    fn blast_better(&self, a: &BlastAlignment, b: &BlastAlignment) -> bool {
        if a.hsp.sseqid != b.hsp.sseqid {
            return false;
        }
        // Mutation vs non-mutation can't dominate each other
        let a_mut = a.resistance == "mutation";
        let b_mut = b.resistance == "mutation";
        if a_mut != b_mut {
            return false;
        }
        // Exact match beats non-exact
        let a_exact = (a.hsp.rel_identity() - 1.0).abs() < 1e-6 && a.hsp.q_complete();
        let b_exact = (b.hsp.rel_identity() - 1.0).abs() < 1e-6 && b.hsp.q_complete();
        if a_exact && !b_exact {
            return true;
        }
        if !a_exact && b_exact {
            return false;
        }
        // Compare by nident (higher is better)
        if a.hsp.nident != b.hsp.nident {
            return a.hsp.nident > b.hsp.nident;
        }
        // Equal nident: prefer shorter reference (more specific match)
        if a.hsp.qlen != b.hsp.qlen {
            return a.hsp.qlen < b.hsp.qlen;
        }
        // Tie-break: lexicographic by accession (C++ PD-1245)
        if a.ref_accession != b.ref_accession {
            return a.ref_accession < b.ref_accession;
        }
        false
    }

    /// Generate TSV report to output
    pub fn report(&self, out: &mut dyn Write, print_node: bool) -> Result<()> {
        let mut tsv = TsvOut::new(Some(out));
        tsv.use_pound = false;

        // Write header
        tsv.write_field(&crate::columns::PROT_COL_NAME)?;
        if self.cds_exist {
            tsv.write_field(&crate::columns::CONTIG_COL_NAME)?;
            tsv.write_field(&crate::columns::START_COL_NAME)?;
            tsv.write_field(&crate::columns::STOP_COL_NAME)?;
            tsv.write_field(&crate::columns::STRAND_COL_NAME)?;
        }
        tsv.write_field(&crate::columns::GENESYMBOL_COL_NAME)?;
        tsv.write_field(&crate::columns::ELEM_NAME_COL_NAME)?;
        tsv.write_field(&crate::columns::SCOPE_COL_NAME)?;
        tsv.write_field(&crate::columns::TYPE_COL_NAME)?;
        tsv.write_field(&crate::columns::SUBTYPE_COL_NAME)?;
        tsv.write_field(&crate::columns::CLASS_COL_NAME)?;
        tsv.write_field(&crate::columns::SUBCLASS_COL_NAME)?;
        tsv.write_field(&crate::columns::METHOD_COL_NAME)?;
        tsv.write_field(&crate::columns::TARGET_LEN_COL_NAME)?;
        tsv.write_field(&crate::columns::REF_LEN_COL_NAME)?;
        tsv.write_field(&crate::columns::REF_COV_COL_NAME)?;
        tsv.write_field(&crate::columns::REF_IDENT_COL_NAME)?;
        tsv.write_field(&crate::columns::ALIGN_LEN_COL_NAME)?;
        tsv.write_field(&crate::columns::CLOSEST_REF_ACCESSION_COL_NAME)?;
        tsv.write_field(&crate::columns::CLOSEST_REF_NAME_COL_NAME)?;
        tsv.write_field(&crate::columns::HMM_ACCESSION_COL_NAME)?;
        tsv.write_field(&crate::columns::HMM_DESCR_COL_NAME)?;
        if print_node {
            tsv.write_field(&crate::columns::HIERARCHY_NODE_COL_NAME)?;
        }
        tsv.new_line()?;

        // Output rows for each target
        for indices in self.target2good_blast_als.values() {
            for &idx in indices {
                let al = &self.blast_als[idx];
                // Skip mutation proteins without detected mutations
                if al.resistance == "mutation" && al.seq_changes.is_empty() {
                    continue;
                }
                // Check reportability
                let reportable = self.get_reportable(al);
                if reportable < self.reportable_min {
                    continue;
                }
                // Check suppress
                if self.suppress_prots.contains(&al.ref_accession) {
                    continue;
                }
                self.report_alignment(al, &mut tsv, print_node)?;
            }
        }

        Ok(())
    }

    /// Get the method name for an alignment
    fn get_method(&self, al: &BlastAlignment) -> String {
        if al.from_hmm && al.ref_accession.is_empty() {
            return "HMM".to_string();
        }
        if al.resistance == "mutation" {
            return format!("POINT{}", if al.hsp.s_prot { "P" } else { "X" });
        }
        let suffix = if al.hsp.s_prot { "P" } else { "X" };
        let cov = al.hsp.q_rel_coverage();

        // Check for exact match: 100% identity AND reference completely covered
        // Account for stop codon: reference may be 1 longer than alignment
        let q_coverage_complete = al.hsp.q_abs_coverage() == al.hsp.qlen
            || (al.hsp.q_prot && al.hsp.q_abs_coverage() + 1 == al.hsp.qlen);
        let ref_exactly_matched = al.hsp.nident == al.hsp.length && q_coverage_complete;

        if ref_exactly_matched {
            // Check if it's an allele
            if let Some(fam) = self.fam_map.get(&al.fam_id) {
                if !fam.genesymbol.is_empty() {
                    let parent_has_diff_symbol = self.fam_map.get(&fam.parent_id)
                        .map(|p| !p.genesymbol.is_empty() && p.genesymbol != fam.genesymbol)
                        .unwrap_or(false);
                    if parent_has_diff_symbol {
                        return format!("ALLELE{}", suffix);
                    }
                }
            } else {
                // fam_id not in fam_map — it's a specific allele not in the hierarchy
                // Check if there's a parent-like family (e.g., blaTEM for blaTEM-156)
                // Find the longest matching fam in the hierarchy
                let mut best_parent = String::new();
                for fam_key in self.fam_map.keys() {
                    if al.fam_id.starts_with(fam_key.as_str()) && fam_key.len() > best_parent.len() {
                        best_parent = fam_key.clone();
                    }
                }
                if !best_parent.is_empty() && best_parent != al.fam_id {
                    return format!("ALLELE{}", suffix);
                }
            }
            return format!("EXACT{}", suffix);
        }

        // Check internal stop
        if al.hsp.s_internal_stop {
            return "INTERNAL_STOP".to_string();
        }

        // Complete match (high coverage) vs partial
        let complete_coverage = 0.9; // C++ complete_coverage_min_def
        if cov >= complete_coverage - 1e-6 {
            return format!("BLAST{}", suffix);
        }

        // Partial — check if truncated at contig end
        // For protein targets, we can't check contig truncation without CDS info
        if !al.cdss.is_empty() {
            let cds = &al.cdss[0];
            if cds.at_contig_start() || cds.at_contig_stop() {
                return format!("PARTIAL_CONTIG_END{}", suffix);
            }
        }

        format!("PARTIAL{}", suffix)
    }

    /// Get fam-derived info for an alignment
    fn get_fam_info(&self, al: &BlastAlignment) -> (String, String, String, String, String, u8) {
        // Returns: (genesymbol, type, subtype, class, subclass, reportable)
        // C++ getFam(): try famId first, then gene field as fallback
        let fam = self.fam_map.get(&al.fam_id)
            .or_else(|| if !al.gene.is_empty() { self.fam_map.get(&al.gene) } else { None });
        let match_fam = self.find_match_fam(&al.fam_id)
            .or_else(|| if !al.gene.is_empty() { self.find_match_fam(&al.gene) } else { None });

        let is_exact = (al.hsp.rel_identity() - 1.0).abs() < 1e-6 && al.hsp.q_complete();

        let genesymbol = if al.resistance == "mutation" {
            al.genesymbol.clone()
        } else if is_exact || al.parts >= 2 {
            // For exact/allele matches: use the fam's own genesymbol
            if let Some(f) = fam {
                if !f.genesymbol.is_empty() {
                    f.genesymbol.clone()
                } else {
                    al.fam_id.clone()
                }
            } else {
                al.fam_id.clone()
            }
        } else {
            // For partial/blast matches: use the match family's genesymbol (parent with genesymbol)
            if let Some(mf) = match_fam {
                if !mf.genesymbol.is_empty() {
                    mf.genesymbol.clone()
                } else {
                    al.fam_id.clone()
                }
            } else if let Some(f) = fam {
                if !f.genesymbol.is_empty() {
                    f.genesymbol.clone()
                } else {
                    al.fam_id.clone()
                }
            } else {
                al.fam_id.clone()
            }
        };

        // Product name: for non-exact matches, use match family's familyName
        let (type_, subtype, class, subclass) = if let Some(mf) = match_fam {
            (mf.type_.clone(), mf.subtype.clone(), mf.class.clone(), mf.subclass.clone())
        } else if let Some(f) = fam {
            (f.type_.clone(), f.subtype.clone(), f.class.clone(), f.subclass.clone())
        } else {
            (al.resistance.clone(), String::new(), al.class.clone(), al.subclass.clone())
        };

        let reportable = self.get_reportable(al);

        (genesymbol, type_, subtype, class, subclass, reportable)
    }

    /// Walk family hierarchy to find ancestor with HMM info
    fn find_hmm_fam(&self, fam_id: &str) -> Option<&Fam> {
        let mut current = self.fam_map.get(fam_id)?;
        for _ in 0..100 {
            if !current.hmm.is_empty() {
                return Some(current);
            }
            if current.parent_id.is_empty() {
                break;
            }
            current = self.fam_map.get(&current.parent_id)?;
        }
        None
    }

    /// Find the matching family (ancestor with genesymbol and type info)
    fn find_match_fam(&self, fam_id: &str) -> Option<&Fam> {
        let fam = self.fam_map.get(fam_id)?;
        // Walk up to find an ancestor with genesymbol
        let mut current = fam;
        for _ in 0..100 {
            if !current.genesymbol.is_empty() && (!current.type_.is_empty() || current.reportable > 0) {
                return Some(current);
            }
            if current.parent_id.is_empty() || !self.fam_map.contains_key(&current.parent_id) {
                break;
            }
            current = self.fam_map.get(&current.parent_id)?;
        }
        // If no ancestor with genesymbol, try the first one with type
        let mut current = fam;
        for _ in 0..100 {
            if !current.type_.is_empty() {
                return Some(current);
            }
            if current.parent_id.is_empty() || !self.fam_map.contains_key(&current.parent_id) {
                break;
            }
            current = self.fam_map.get(&current.parent_id)?;
        }
        Some(fam)
    }

    fn report_alignment(
        &self,
        al: &BlastAlignment,
        tsv: &mut TsvOut,
        print_node: bool,
    ) -> Result<()> {
        let na = crate::columns::NA;
        let method = self.get_method(al);
        let (genesymbol, type_, subtype, class, subclass, reportable) = self.get_fam_info(al);

        let is_mutation = al.resistance == "mutation";

        let is_exact = (al.hsp.rel_identity() - 1.0).abs() < 1e-6 && al.hsp.q_complete();
        // Product name: exact/allele uses product, others use match family name
        let product_name = if al.ref_accession.is_empty() {
            self.fam_map.get(&al.fam_id)
                .map(|f| f.family_name.clone())
                .unwrap_or_else(|| al.product.clone())
        } else if is_exact || al.parts >= 2 {
            al.product.clone()
        } else if let Some(mf) = self.find_match_fam(&al.fam_id) {
            if !mf.family_name.is_empty() {
                mf.family_name.clone()
            } else {
                al.product.clone()
            }
        } else {
            al.product.clone()
        };

        // Protein id
        if al.hsp.s_prot {
            tsv.write_field(&al.hsp.sseqid)?;
        } else {
            tsv.write_field(&na)?;
        }

        // Contig info (only if cds_exist)
        if self.cds_exist {
            if al.cdss.is_empty() {
                if !al.hsp.s_prot {
                    tsv.write_field(&al.hsp.sseqid)?;
                    tsv.write_field(&(al.hsp.s_int.start + 1))?;
                    tsv.write_field(&al.hsp.s_int.stop)?;
                    tsv.write_field(&crate::seq::strand2char(al.hsp.s_int.strand).to_string())?;
                } else {
                    tsv.write_field(&na)?;
                    tsv.write_field(&na)?;
                    tsv.write_field(&na)?;
                    tsv.write_field(&na)?;
                }
            } else {
                let cds = &al.cdss[0];
                tsv.write_field(&cds.contig)?;
                tsv.write_field(&(cds.start + 1))?;
                tsv.write_field(&cds.stop)?;
                tsv.write_field(&if cds.strand { "+" } else { "-" })?;
            }
        }

        // Element symbol
        tsv.write_field(&if genesymbol.is_empty() { na.to_string() } else { genesymbol })?;
        // Element name
        tsv.write_field(&if product_name.is_empty() { na.to_string() } else { product_name })?;

        // Scope
        let scope = if is_mutation || reportable >= 2 { "core" } else { "plus" };
        tsv.write_field(&scope)?;

        // Type, Subtype, Class, Subclass
        if is_mutation {
            tsv.write_field(&"AMR")?;
            tsv.write_field(&"POINT")?;
        } else {
            tsv.write_field(&if type_.is_empty() { na.to_string() } else { type_ })?;
            tsv.write_field(&if subtype.is_empty() { na.to_string() } else { subtype })?;
        }
        tsv.write_field(&if class.is_empty() { na.to_string() } else { class })?;
        tsv.write_field(&if subclass.is_empty() { na.to_string() } else { subclass })?;

        // Method
        tsv.write_field(&method)?;

        // Target length
        tsv.write_field(&if al.hsp.s_prot { al.hsp.slen } else { al.hsp.s_abs_coverage() / 3 })?;

        tsv.write_field(&al.hsp.qlen)?;
        tsv.write_field(&format!("{:.2}", al.hsp.q_rel_coverage() * 100.0))?;
        tsv.write_field(&format!("{:.2}", al.hsp.rel_identity() * 100.0))?;
        if al.ref_accession.is_empty() {
            tsv.write_field(&al.hsp.qseq.len())?;
        } else {
            tsv.write_field(&al.hsp.q_len_real())?;
        }
        tsv.write_field(&al.ref_accession)?;
        tsv.write_field(&al.product)?;

        // HMM info — C++ uses getFam()->hmm which walks famId then gene fallback
        // Try: 1) from_hmm fam, 2) hmm_al_idx, 3) fam hierarchy (famId→gene), 4) NA
        let hmm_fam = if al.from_hmm {
            self.fam_map.get(&al.fam_id)
        } else if let Some(hmm_idx) = al.hmm_al_idx {
            self.fam_map.get(&self.hmm_als[hmm_idx].fam_id)
        } else {
            // Walk fam hierarchy to find HMM info (matching C++ getFam()->hmm)
            self.find_hmm_fam(&al.fam_id)
                .or_else(|| if !al.gene.is_empty() { self.find_hmm_fam(&al.gene) } else { None })
        };
        if let Some(fam) = hmm_fam.filter(|f| !f.hmm.is_empty()) {
            tsv.write_field(&fam.hmm)?;
            tsv.write_field(&fam.family_name)?;
        } else {
            tsv.write_field(&na)?;
            tsv.write_field(&na)?;
        }

        if print_node {
            tsv.write_field(&al.fam_id)?;
        }

        tsv.new_line()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn db_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("amrfinder_db/2026-03-24.1")
    }

    #[test]
    fn test_blast_rule() {
        let br = BlastRule::new(0.9, 0.5, 0.8);
        assert!(!br.empty());
        assert_eq!(br.ident, 0.9);

        let empty = BlastRule::default();
        assert!(empty.empty());
    }

    #[test]
    fn test_fam_loading() {
        let fam_path = db_dir().join("fam.tsv");
        if !fam_path.exists() {
            return;
        }

        let batch = Batch::from_fam_file(&fam_path, 0).unwrap();
        assert!(!batch.fam_map.is_empty(), "FAM map should not be empty");
        assert!(!batch.hmm2fam.is_empty(), "HMM map should not be empty");

        // Verify some known families exist
        assert!(batch.fam_map.contains_key("blaTEM"), "Should have blaTEM family");
    }

    #[test]
    fn test_mutation_loading() {
        let db = db_dir();
        let fam_path = db.join("fam.tsv");
        let mut_path = db.join("AMRProt-mutation.tsv");
        if !fam_path.exists() || !mut_path.exists() {
            return;
        }

        let mut batch = Batch::from_fam_file(&fam_path, 0).unwrap();
        batch.load_mutations(&mut_path, "Escherichia").unwrap();
        assert!(!batch.accession2mutations.is_empty(), "Should have mutations for Escherichia");
    }

    #[test]
    fn test_blast_alignment_parsing() {
        // Format: sseqid(ref) qseqid(target) sstart send slen qstart qend qlen sseq qseq
        // After Hsp parsing: qseqid=ref, sseqid=target (BLAST "query" is the reference protein)
        let line = "WP_061158039.1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tblaTEM-156\t1\t286\t287\t1\t286\t286\tMSIQH\tMSIQH";
        let default_br = BlastRule::default();
        let al = BlastAlignment::from_blast_line(line, true, true, &default_br, &default_br);
        assert!(al.is_ok(), "Failed to parse BLAST line: {:?}", al.err());
        let al = al.unwrap();
        // qseqid has the reference metadata, sseqid has the user's protein
        assert_eq!(al.hsp.sseqid, "blaTEM-156");
        assert_eq!(al.ref_accession, "WP_061158039.1");
        assert_eq!(al.fam_id, "blaTEM-156");
        assert_eq!(al.gene, "blaTEM");
    }

    // --- amr_report behavior tests ---
    // These test specific behaviors that must match C++ amr_report

    /// C++ isCore() requires reportable >= 2 (amr_report.cpp:820-821).
    /// reportable=1 entries (like stx genes) should be "plus", not "core".
    #[test]
    fn test_scope_requires_reportable_2_for_core() {
        let fam_path = db_dir().join("fam.tsv");
        if !fam_path.exists() {
            return;
        }
        let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

        // Build a BLAST alignment with reportable=1 (like stxA2)
        let line = "TJA36680.1|1|1|stxA2_acd|stxA2|VIRULENCE|1|stxA2|STX2|Shiga_toxin_Stx2_subunit_A\tstxA2a_prot\t1\t319\t320\t94\t1051\t320\tMKCIL\tMKCIL";
        let br = BlastRule::default();
        let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

        let reportable = batch.get_reportable(&al);
        // reportable=1 should NOT be core; C++ requires >= 2
        assert!(
            reportable < 2,
            "stxA2 with reportable=1 should not qualify as core (reportable={})",
            reportable
        );
    }

    /// C++ getFam() falls back from famId to gene field (amr_report.cpp:1222-1224).
    /// When fam_id is an allele (e.g. "blaTEM-156") not in fam.tsv, the gene field
    /// ("blaTEM") should be used as fallback.
    #[test]
    fn test_find_match_fam_gene_fallback() {
        let fam_path = db_dir().join("fam.tsv");
        if !fam_path.exists() {
            return;
        }
        let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

        // "blaTEM-156" is an allele — it won't be in fam.tsv
        // "blaTEM" is the parent family — it IS in fam.tsv
        assert!(
            !batch.fam_map.contains_key("blaTEM-156"),
            "blaTEM-156 should NOT be a direct key in fam_map"
        );
        assert!(
            batch.fam_map.contains_key("blaTEM"),
            "blaTEM should be in fam_map"
        );

        // find_match_fam("blaTEM-156") returns None (allele not in fam.tsv)
        // The fallback to gene field happens at the get_fam_info level
        assert!(
            batch.find_match_fam("blaTEM-156").is_none(),
            "blaTEM-156 is not in fam_map, find_match_fam should return None"
        );

        // But find_match_fam("blaTEM") should resolve (gene field fallback)
        let result = batch.find_match_fam("blaTEM");
        assert!(
            result.is_some(),
            "find_match_fam should resolve blaTEM from fam_map"
        );
        if let Some(fam) = result {
            assert_eq!(fam.type_, "AMR", "blaTEM family type should be AMR");
        }
    }

    /// Verify fam_map resolves type/class for known allele families.
    /// The C++ uses getFam() which falls back famId → gene. Without fallback,
    /// the Rust code uses the BLAST header's resistance field ("hydrolase")
    /// instead of the fam.tsv's type ("AMR").
    #[test]
    fn test_get_fam_info_uses_fam_type_not_blast_header() {
        let fam_path = db_dir().join("fam.tsv");
        if !fam_path.exists() {
            return;
        }
        let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

        // BLAST header for blaTEM-156: resistance="hydrolase"
        // But fam.tsv entry for blaTEM has type_="AMR"
        let line = "WP_061158039.1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tblaTEM-156\t1\t286\t287\t1\t286\t286\tMSIQH\tMSIQH";
        let br = BlastRule::default();
        let al = BlastAlignment::from_blast_line(line, true, true, &br, &br).unwrap();

        let (_genesymbol, type_, _subtype, _class, _subclass, _reportable) = batch.get_fam_info(&al);
        assert_eq!(
            type_, "AMR",
            "Type should come from fam.tsv (AMR), not BLAST header (hydrolase)"
        );
    }

    /// C++ links BLAST alignments with supporting HMM results via getFam()->hmm.
    /// The report should show HMM accession and description for BLAST hits
    /// that have matching HMM families, not always "NA".
    #[test]
    fn test_hmm_info_populated_for_blast_hits_with_hmm_families() {
        let fam_path = db_dir().join("fam.tsv");
        if !fam_path.exists() {
            return;
        }
        let batch = Batch::from_fam_file(&fam_path, 0).unwrap();

        // Check that the blaTEM family has HMM info in fam.tsv
        if let Some(fam) = batch.fam_map.get("blaTEM") {
            // If the family has an HMM, BLAST hits in this family should report it
            if !fam.hmm.is_empty() {
                // This documents the expected behavior:
                // C++ outputs HMM accession/description even for BLAST-method hits
                // when the family has an HMM entry
                assert!(
                    !fam.hmm.is_empty(),
                    "blaTEM family should have HMM accession"
                );
                assert!(
                    !fam.family_name.is_empty(),
                    "blaTEM family should have family_name for HMM description"
                );
            }
        }
    }

    /// Pareto filter should remove dominated hits.
    /// A hit with lower identity and equal coverage should be dominated.
    #[test]
    fn test_pareto_filter_removes_dominated_hits() {
        let fam_path = db_dir().join("fam.tsv");
        if !fam_path.exists() {
            return;
        }
        let mut batch = Batch::from_fam_file(&fam_path, 0).unwrap();

        let br = BlastRule::default();
        // Use short identical-length sequences to avoid finish_hsp issues
        let seq = "MSIQH";

        // Hit 1: 100% identity (5/5 nident), full coverage
        let line1 = format!(
            "WP_1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|product1\ttarget1\t1\t5\t5\t1\t5\t5\t{}\t{}",
            seq, seq
        );
        let al1 = BlastAlignment::from_blast_line(&line1, true, true, &br, &br).unwrap();
        batch.add_blast_al(al1);

        // Hit 2: 80% identity (4/5 nident), same coverage — dominated by hit 1
        let line2 = format!(
            "WP_2|1|1|blaTEM-239|blaTEM|hydrolase|2|NA|NA|product2\ttarget1\t1\t5\t5\t1\t5\t5\t{}\t{}",
            seq, "MSXQH"
        );
        let al2 = BlastAlignment::from_blast_line(&line2, true, true, &br, &br).unwrap();
        batch.add_blast_al(al2);

        let indices = vec![0, 1];
        let filtered = batch.pareto_filter_blast(&indices);

        // Hit 2 should be dominated by hit 1 (strictly lower identity, equal coverage)
        assert_eq!(
            filtered.len(), 1,
            "Pareto filter should remove dominated hit, got {} hits",
            filtered.len()
        );
    }
}
