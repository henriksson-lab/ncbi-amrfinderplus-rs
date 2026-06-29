// amr_report — core reporting engine
// Port of amr_report.cpp

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::alignment::{AmrMutation, SeqChange};
use crate::gff::{Annot, GffType, Locus};
use crate::seq::{Hsp, HspMerge};
use crate::tsv::TsvOut;

const FRAC_DELTA: f64 = 1e-5;
const IDENT_MIN_DEF: f64 = 0.9;
const COMPLETE_COVERAGE_MIN_DEF: f64 = 0.9;
const PARTIAL_COVERAGE_MIN_DEF: f64 = 0.5;
const MISMATCH_TAIL_AA: usize = 10;

// --- BlastRule ---

/// BLAST identity/coverage thresholds for a family
#[derive(Debug, Clone, Default)]
pub struct BlastRule {
    pub ident: f64,
    pub target_coverage: f64,
    pub ref_coverage: f64,
}

impl BlastRule {
    pub fn new(ident: f64, ref_coverage: f64) -> Self {
        let normalize = |value: f64| {
            if value > 1.0 {
                value / 100.0
            } else {
                value
            }
        };
        BlastRule {
            ident: normalize(ident),
            target_coverage: 0.0,
            ref_coverage: normalize(ref_coverage),
        }
    }

    pub fn empty(&self) -> bool {
        self.ident == 0.0
    }

    pub fn save_text<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
        write!(out, "{} {}", self.ident, self.ref_coverage)
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

    pub fn save_text<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
        write!(
            out,
            "{} {} {} {} {}",
            self.hmm, self.tc1, self.tc2, self.family_name, self.reportable
        )
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

impl Susceptible {
    pub fn save_text<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            self.genesymbol, self.cutoff, self.class, self.subclass, self.name
        )
    }
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
    pub fn save_text<W: Write>(
        &self,
        out: &mut W,
        fam_map: &HashMap<String, Fam>,
    ) -> std::io::Result<()> {
        let hmm = fam_map
            .get(&self.fam_id)
            .map(|fam| fam.hmm.as_str())
            .unwrap_or("");
        write!(
            out,
            "{} {} {} {}",
            self.sseqid, self.score1, self.score2, hmm
        )
    }

    pub fn good(&self, fam_map: &HashMap<String, Fam>) -> bool {
        if let Some(fam) = fam_map.get(&self.fam_id) {
            self.score1 >= fam.tc1 && self.score2 >= fam.tc2
        } else {
            false
        }
    }

    fn better_eq(
        &self,
        other: &HmmAlignment,
        criterion: u8,
        fam_map: &HashMap<String, Fam>,
        equidistant: bool,
    ) -> bool {
        if self.sseqid != other.sseqid {
            return false;
        }
        let Some(fam) = fam_map.get(&self.fam_id) else {
            return false;
        };
        let Some(other_fam) = fam_map.get(&other.fam_id) else {
            return false;
        };
        match criterion {
            0 => fam.descendant_of(&other_fam.id, fam_map),
            1 => {
                if other.score1 < self.score1 {
                    return true;
                }
                if self.score1 < other.score1 {
                    return false;
                }
                if other_fam.tc1 < fam.tc1 {
                    return true;
                }
                if fam.tc1 < other_fam.tc1 {
                    return false;
                }
                if !equidistant {
                    if fam.id < other_fam.id {
                        return true;
                    }
                    if other_fam.id < fam.id {
                        return false;
                    }
                }
                true
            }
            _ => false,
        }
    }

    fn better(
        &self,
        other: &HmmAlignment,
        criterion: u8,
        fam_map: &HashMap<String, Fam>,
        equidistant: bool,
    ) -> bool {
        self.better_eq(other, criterion, fam_map, equidistant)
            && !other.better_eq(self, criterion, fam_map, equidistant)
    }

    pub fn better_blast(&self, other: &BlastAlignment, fam_map: &HashMap<String, Fam>) -> bool {
        if !other.hsp.s_prot || self.sseqid != other.hsp.sseqid {
            return false;
        }
        let Some(fam) = fam_map.get(&self.fam_id) else {
            return false;
        };
        let Some(other_fam) = other.get_match_fam(fam_map) else {
            return false;
        };
        fam.id != other_fam.id && fam.descendant_of(&other_fam.id, fam_map)
    }
}

// --- BlastAlignment ---

/// BLAST alignment result with full AMR annotation
#[derive(Debug, Clone)]
pub struct BlastAlignment {
    pub hsp: Hsp,
    pub partial_dna: bool,
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
    pub br_fam_id: Option<String>,
    pub br_fam_checked: bool,
    pub complete_br: BlastRule,
    pub partial_br: BlastRule,
    pub cdss: Vec<Locus>,
    pub hmm_al_idx: Option<usize>,
    pub susceptible_idx: Option<usize>,
    pub seq_changes: Vec<SeqChange>,
    pub ref_mutations: Vec<AmrMutation>,
    pub ref_mutation: AmrMutation,
    // Fusion protein info
    pub fusion_ids: Vec<usize>,
    pub fusion_redundant: bool,
}

impl BlastAlignment {
    fn from_hmm_alignment(
        hmm_al: &HmmAlignment,
        fam: Option<&Fam>,
        best: Option<&BlastAlignment>,
        domain: &HmmDomain,
    ) -> Self {
        let mut hsp = if let Some(best) = best {
            let mut hsp = best.hsp.clone();
            hsp.sseqid = hmm_al.sseqid.clone();
            hsp
        } else {
            let align_len = domain.hmm_stop - domain.hmm_start;
            let mut hsp = Hsp::default();
            hsp.sseqid = hmm_al.sseqid.clone();
            hsp.qseqid = hmm_al.fam_id.clone();
            hsp.slen = domain.seq_len;
            hsp.s_int = crate::seq::Interval::new(domain.seq_start, domain.seq_stop, 0);
            hsp.qlen = domain.hmm_len;
            hsp.q_int = crate::seq::Interval::new(domain.hmm_start, domain.hmm_stop, 0);
            hsp.q_prot = true;
            hsp.s_prot = true;
            hsp.a_prot = true;
            hsp.length = align_len;
            hsp.nident = (align_len as f64 * 0.7) as usize;
            hsp
        };
        hsp.s_prot = true;

        BlastAlignment {
            hsp,
            partial_dna: false,
            from_hmm: true,
            ref_accession: best
                .map(|best| best.ref_accession.clone())
                .unwrap_or_default(),
            part: 1,
            parts: 1,
            fam_id: hmm_al.fam_id.clone(),
            gene: hmm_al.fam_id.clone(),
            resistance: String::new(),
            class: fam.map(|f| f.class.clone()).unwrap_or_default(),
            subclass: fam.map(|f| f.subclass.clone()).unwrap_or_default(),
            product: fam.map(|f| f.family_name.clone()).unwrap_or_default(),
            reportable: fam.map(|f| f.reportable).unwrap_or(0),
            genesymbol: fam.map(|f| f.genesymbol.clone()).unwrap_or_default(),
            method: "HMM".to_string(),
            complete_br: BlastRule::default(),
            partial_br: BlastRule::default(),
            br_fam_id: None,
            br_fam_checked: false,
            cdss: best.map(|best| best.cdss.clone()).unwrap_or_default(),
            hmm_al_idx: None,
            susceptible_idx: None,
            seq_changes: Vec::new(),
            ref_mutations: Vec::new(),
            ref_mutation: crate::alignment::AmrMutation::default(),
            fusion_ids: Vec::new(),
            fusion_redundant: false,
        }
    }

    /// Parse from BLAST qseqid format:
    /// accession|part|parts|famId|gene|resistance|reportable|subclass|classS|product
    pub fn from_blast_line(
        line: &str,
        q_prot: bool,
        s_prot: bool,
        default_complete_br: &BlastRule,
        default_partial_br: &BlastRule,
    ) -> Result<Self> {
        let a_prot = q_prot || s_prot;
        let hsp = Hsp::from_blast_line(line, q_prot, s_prot, a_prot, q_prot, true)?;

        let (
            mut ref_accession,
            part,
            parts_count,
            fam_id,
            gene,
            resistance,
            reportable,
            subclass,
            class,
            genesymbol,
            product,
        ) = (|| -> Result<_> {
            // C++ parses qseqid by repeated rfindSplit('|') from the right.
            let mut rest = hsp.qseqid.as_str();
            let mut fields = Vec::with_capacity(9);
            for _ in 0..9 {
                let Some(pos) = rest.rfind('|') else {
                    bail!("missing qseqid field")
                };
                fields.push(&rest[pos + 1..]);
                rest = &rest[..pos];
            }
            let product = fields[0].to_string().replace('_', " ");
            let class = fields[1].to_string().replace('_', " ");
            let subclass = fields[2].to_string().replace('_', " ");
            let reportable = fields[3].parse::<u8>()?;
            let resistance = fields[4].to_string();
            let gene = fields[5].to_string();
            let fam_id = fields[6].to_string();
            let parts_count = fields[7].parse::<usize>()?;
            let part = fields[8].parse::<usize>()?;
            let ref_accession = rest.to_string();
            Ok((
                ref_accession,
                part,
                parts_count,
                fam_id,
                gene,
                resistance,
                reportable,
                subclass,
                class,
                String::new(),
                product,
            ))
        })()
        .with_context(|| format!("Bad AMRFinder database\n{line}"))?;

        let ref_mutation =
            parse_declarative_ref_mutation(&mut ref_accession, resistance == "mutation")
                .with_context(|| format!("Bad AMRFinder database\n{line}"))?;
        let ref_mutations = if ref_mutation.empty() {
            Vec::new()
        } else {
            vec![ref_mutation.clone()]
        };

        let mut al = BlastAlignment {
            hsp,
            partial_dna: false,
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
            br_fam_id: None,
            br_fam_checked: false,
            complete_br: default_complete_br.clone(),
            partial_br: default_partial_br.clone(),
            cdss: Vec::new(),
            hmm_al_idx: None,
            susceptible_idx: None,
            seq_changes: Vec::new(),
            ref_mutations,
            ref_mutation,
            fusion_ids: Vec::new(),
            fusion_redundant: false,
        };
        if s_prot {
            al.finish()?;
        }
        Ok(al)
    }

    /// Check if this alignment meets BLAST thresholds.
    pub fn good(&self, susceptible: Option<&Susceptible>) -> bool {
        if self.ref_accession.is_empty() {
            return true;
        }
        if self.is_mutation_prot() && self.seq_changes.is_empty() {
            return false;
        }
        if self.is_susceptible_prot() && susceptible.is_none() {
            return false;
        }
        if let Some(susceptible) = susceptible {
            if susceptible.cutoff != 0.0 {
                if self.hsp.rel_identity() > susceptible.cutoff + FRAC_DELTA {
                    return false;
                }
            } else if self.hsp.disrs.is_empty() {
                return false;
            }
        }
        if self.partial()
            && (self.parts > 1 || (!self.truncated_cds() && self.hsp.q_int.len() <= 35))
        {
            return false;
        }

        if self.br_fam_id.is_some() {
            return true;
        }
        if self.br_fam_checked && self.in_fam() {
            return false;
        }
        let br = if self.partial() {
            &self.partial_br
        } else {
            &self.complete_br
        };
        if !self.br_fam_checked && br.empty() {
            return true;
        }
        self.pass_blast_rule(br)
    }

    pub fn is_mutation_prot(&self) -> bool {
        self.resistance == "mutation"
    }

    pub fn is_susceptible_prot(&self) -> bool {
        self.resistance == "susceptible"
    }

    fn is_strong_susceptible_prot(&self, batch: &Batch) -> bool {
        self.is_susceptible_prot()
            && batch
                .accession2susceptible
                .get(&self.ref_accession)
                .is_some_and(|susceptible| susceptible.cutoff == 0.0)
    }

    pub fn in_fam(&self) -> bool {
        !self.is_mutation_prot() && !self.is_susceptible_prot()
    }

    pub fn allele(&self) -> bool {
        self.fam_id != self.gene && self.parts == 1
    }

    pub fn ref_effective_len(&self) -> usize {
        if self.partial_dna {
            self.hsp.q_int.len()
        } else {
            self.hsp.qlen
        }
    }

    pub fn partial(&self) -> bool {
        self.hsp.q_rel_coverage() < COMPLETE_COVERAGE_MIN_DEF - FRAC_DELTA
    }

    pub fn ref_prot_exactly_matched(&self, target_complete: bool) -> bool {
        self.hsp.q_prot
            && self.hsp.qlen != 0
            && self.hsp.nident == self.hsp.qlen
            && self.hsp.nident == self.hsp.sseq.len()
            && self.hsp.c_complete != Some(false)
            && (!self.hsp.s_prot
                || !target_complete
                || self.hsp.qlen + usize::from(self.hsp.c_complete == Some(true)) == self.hsp.slen)
            && self.hsp.disrs.is_empty()
    }

    pub fn allele_match(&self) -> bool {
        self.ref_prot_exactly_matched(true) && self.allele()
    }

    pub fn get_cds_start(&self) -> usize {
        self.cdss.iter().map(|cds| cds.start).min().unwrap_or(0)
    }

    pub fn get_cds_stop(&self) -> usize {
        self.cdss.iter().map(|cds| cds.stop).max().unwrap_or(0)
    }

    pub fn same_match(&self, other: &BlastAlignment) -> bool {
        self.hsp.sseqid == other.hsp.sseqid
            && self.hsp.s_int == other.hsp.s_int
            && self.get_cds_start() == other.get_cds_start()
            && self.get_cds_stop() == other.get_cds_stop()
            && self.ref_accession == other.ref_accession
    }

    fn less(&self, other: &BlastAlignment) -> bool {
        self.hsp
            .sseqid
            .cmp(&other.hsp.sseqid)
            .then(self.hsp.s_int.start.cmp(&other.hsp.s_int.start))
            .then(self.hsp.s_int.stop.cmp(&other.hsp.s_int.stop))
            .then(self.get_cds_start().cmp(&other.get_cds_start()))
            .then(self.get_cds_stop().cmp(&other.get_cds_stop()))
            .then(self.ref_accession.cmp(&other.ref_accession))
            .then(self.part.cmp(&other.part))
            .is_lt()
    }

    fn matches_cds(&self, other: &BlastAlignment) -> bool {
        assert!(self.hsp.s_prot);
        assert!(!other.hsp.s_prot);
        assert!(!self.cdss.is_empty());

        for cds in &self.cdss {
            if cds.contig == other.hsp.sseqid
                && cds.strand == (other.hsp.s_int.strand == 1)
                && !cds.cross_origin
            {
                let mut prot_start = cds.start;
                let mut prot_stop = cds.stop;
                if cds.strand {
                    assert!(prot_stop > 3);
                    prot_stop -= 3;
                } else {
                    prot_start += 3;
                }
                assert!(prot_start < prot_stop);

                let dna_start = other.hsp.s_int.start;
                let dna_stop = other.hsp.s_int.stop;
                assert!(dna_start < dna_stop);

                let intersection_start = prot_start.max(dna_start);
                let intersection_stop = prot_stop.min(dna_stop);
                let union_start = prot_start.min(dna_start);
                let union_stop = prot_stop.max(dna_stop);

                if (intersection_start < intersection_stop
                    && (intersection_stop - intersection_start) as f64
                        / (prot_stop - prot_start) as f64
                        > 0.75)
                    || (intersection_start - union_start) + (union_stop - intersection_stop)
                        <= 60 * 3
                    || (prot_start <= dna_start + MISMATCH_TAIL_AA * 3
                        && dna_stop <= prot_stop + MISMATCH_TAIL_AA * 3
                        && other.partial())
                    || (dna_start <= prot_start + MISMATCH_TAIL_AA * 3
                        && prot_stop <= dna_stop + MISMATCH_TAIL_AA * 3
                        && self.partial())
                {
                    return true;
                }
            }
        }

        false
    }

    fn inside_eq(&self, other: &BlastAlignment, batch: &Batch) -> bool {
        if self.hsp.s_prot != other.hsp.s_prot || self.hsp.sseqid != other.hsp.sseqid {
            return false;
        }
        if self.hsp.s_int.strand != other.hsp.s_int.strand {
            return false;
        }
        let mismatch_tail = MISMATCH_TAIL_AA * self.hsp.a2s;
        if self.hsp.s_int.start + mismatch_tail >= other.hsp.s_int.start
            && self.hsp.s_int.stop <= other.hsp.s_int.stop + mismatch_tail
        {
            return true;
        }
        if !self.partial_pseudo()
            || self.is_mutation_prot()
            || other.is_mutation_prot()
            || batch.fusion_2_gene_symbols(self) != batch.fusion_2_gene_symbols(other)
        {
            return false;
        }
        let pseudo_overlap = MISMATCH_TAIL_AA * self.hsp.a2s * 2;
        (self.hsp.s_int.start < other.hsp.s_int.start
            && self.hsp.s_int.stop >= other.hsp.s_int.start + pseudo_overlap)
            || (self.hsp.s_int.stop > other.hsp.s_int.stop
                && self.hsp.s_int.start + pseudo_overlap <= other.hsp.s_int.stop)
    }

    fn fusion_overrides(&self, other: &BlastAlignment, batch: &Batch) -> bool {
        self.hsp.s_prot == other.hsp.s_prot
            && self.parts > 1
            && other.parts == 1
            && other.inside_eq(self, batch)
            && (self.hsp.s_int.start + crate::seq::CDS_PEPTIDE_SIZE_MIN * self.hsp.a2s
                <= other.hsp.s_int.start
                || self.hsp.s_int.stop
                    >= other.hsp.s_int.stop + crate::seq::CDS_PEPTIDE_SIZE_MIN * self.hsp.a2s)
    }

    fn get_fam<'a>(&self, fam_map: &'a HashMap<String, Fam>) -> Option<&'a Fam> {
        fam_map.get(&self.fam_id).or_else(|| {
            if self.gene.is_empty() {
                None
            } else {
                fam_map.get(&self.gene)
            }
        })
    }

    fn get_match_fam<'a>(&self, fam_map: &'a HashMap<String, Fam>) -> Option<&'a Fam> {
        if self.from_hmm {
            self.get_fam(fam_map)
        } else {
            self.br_fam_id
                .as_ref()
                .and_then(|fam_id| fam_map.get(fam_id))
        }
    }

    fn better_hmm(
        &self,
        other: &HmmAlignment,
        fam_map: &HashMap<String, Fam>,
        blast_als: &[BlastAlignment],
    ) -> bool {
        if self.hsp.s_prot {
            if self.hsp.sseqid != other.sseqid {
                return false;
            }
        } else {
            let Some(hmm_blast_idx) = other.blast_al_idx else {
                return false;
            };
            let Some(hmm_blast_al) = blast_als.get(hmm_blast_idx) else {
                return false;
            };
            if !hmm_blast_al.matches_cds(self) {
                return false;
            }
        }

        self.ref_prot_exactly_matched(true)
            || self
                .get_match_fam(fam_map)
                .zip(fam_map.get(&other.fam_id))
                .is_some_and(|(fam, other_fam)| fam.descendant_of(&other_fam.id, fam_map))
    }

    fn has_mutation(&self) -> bool {
        self.seq_changes.iter().any(SeqChange::has_mutation)
    }

    fn get_mutation_symbols(&self) -> HashSet<String> {
        let mut mutation_symbols = HashSet::new();
        for seq_change in &self.seq_changes {
            if !seq_change.empty() {
                for &mutation_idx in &seq_change.mutations {
                    if let Some(mutation) = self.ref_mutations.get(mutation_idx) {
                        mutation_symbols.insert(mutation.gene_mutation.clone());
                    } else if mutation_idx == 0 && !self.ref_mutation.gene_mutation.is_empty() {
                        mutation_symbols.insert(self.ref_mutation.gene_mutation.clone());
                    }
                }
            }
        }
        mutation_symbols
    }

    fn get_target_strand(&self, cds: &Locus) -> bool {
        if self.hsp.s_prot {
            if cds.empty() {
                true
            } else {
                cds.strand
            }
        } else {
            self.hsp.s_int.strand == 1
        }
    }

    fn missed_dna_start(&self, cds: &Locus) -> usize {
        3 * if self.get_target_strand(cds) {
            self.hsp.q_int.start
        } else {
            self.hsp.qlen - self.hsp.q_int.stop
        }
    }

    fn missed_dna_stop(&self, cds: &Locus) -> usize {
        3 * if self.get_target_strand(cds) {
            self.hsp.qlen - self.hsp.q_int.stop
        } else {
            self.hsp.q_int.start
        }
    }

    fn truncated(&self, cds: &Locus) -> bool {
        (self.missed_dna_start(cds) > 0
            && if self.hsp.s_prot {
                !cds.empty() && cds.at_contig_start()
            } else {
                self.hsp.s_int.start <= 3
            })
            || (self.missed_dna_stop(cds) > 0
                && if self.hsp.s_prot {
                    !cds.empty() && cds.at_contig_stop()
                } else {
                    self.hsp.slen - self.hsp.s_int.stop <= 3
                })
    }

    fn truncated_cds(&self) -> bool {
        self.cdss.iter().any(|cds| self.truncated(cds))
    }

    fn partial_pseudo(&self) -> bool {
        self.partial() && !self.cdss.is_empty() && !self.truncated_cds()
    }

    fn get_contigs(&self) -> Vec<String> {
        let mut contigs = self
            .cdss
            .iter()
            .map(|locus| locus.contig.clone())
            .collect::<Vec<_>>();
        contigs.sort();
        contigs.dedup();
        contigs
    }

    fn fusion_2_fam_ids(&self, blast_als: &[BlastAlignment]) -> String {
        if self.is_mutation_prot() || self.fusion_ids.is_empty() {
            return self.fam_id.clone();
        }
        self.fusion_ids
            .iter()
            .map(|&idx| blast_als[idx].fam_id.as_str())
            .collect::<Vec<_>>()
            .join("/")
    }

    fn fusion_2_type(&self, batch: &Batch) -> String {
        if self.fusion_ids.is_empty() {
            return batch.get_type(self);
        }
        let mut values = self
            .fusion_ids
            .iter()
            .map(|&idx| batch.get_type(&batch.blast_als[idx]))
            .collect::<Vec<_>>();
        values.sort();
        values.dedup();
        values.join("/")
    }

    fn get_subtype(&self, batch: &Batch) -> String {
        if self.is_susceptible_prot()
            && (self.susceptible_idx.is_some()
                || batch
                    .accession2susceptible
                    .contains_key(&self.ref_accession))
        {
            "AMR".to_string()
        } else {
            let (_, _, subtype, _, _, _) = batch.get_fam_info(self);
            subtype
        }
    }

    fn fusion_2_subtype(&self, batch: &Batch) -> String {
        if self.fusion_ids.is_empty() {
            return self.get_subtype(batch);
        }
        let mut values = self
            .fusion_ids
            .iter()
            .map(|&idx| batch.blast_als[idx].get_subtype(batch))
            .collect::<Vec<_>>();
        values.sort();
        values.dedup();
        values.join("/")
    }

    fn get_class(&self, batch: &Batch) -> String {
        if self.allele_match() && !self.subclass.is_empty() {
            return self.class.clone();
        }
        if self.is_susceptible_prot() {
            if let Some(susceptible) = batch.accession2susceptible.get(&self.ref_accession) {
                return susceptible.class.clone();
            }
        }
        let (_, _, _, class, _, _) = batch.get_fam_info(self);
        class
    }

    fn fusion_2_class(&self, batch: &Batch) -> String {
        if self.fusion_ids.is_empty() {
            return self.get_class(batch);
        }
        let mut values = Vec::new();
        for &idx in &self.fusion_ids {
            for value in batch.blast_als[idx].get_class(batch).split('/') {
                if !value.is_empty() {
                    values.push(value.to_string());
                }
            }
        }
        values.sort();
        values.dedup();
        values.join("/")
    }

    fn get_subclass(&self, batch: &Batch) -> String {
        if self.allele_match() && !self.subclass.is_empty() {
            return self.subclass.clone();
        }
        if self.is_susceptible_prot() {
            if let Some(susceptible) = batch.accession2susceptible.get(&self.ref_accession) {
                return susceptible.subclass.clone();
            }
        }
        let (_, _, _, _, subclass, _) = batch.get_fam_info(self);
        subclass
    }

    fn fusion_2_subclass(&self, batch: &Batch) -> String {
        if self.fusion_ids.is_empty() {
            return self.get_subclass(batch);
        }
        let mut values = Vec::new();
        for &idx in &self.fusion_ids {
            for value in batch.blast_als[idx].get_subclass(batch).split('/') {
                if !value.is_empty() {
                    values.push(value.to_string());
                }
            }
        }
        values.sort();
        values.dedup();
        values.join("/")
    }

    fn fusion_2_hmm_al<'a>(&self, batch: &'a Batch) -> Option<&'a HmmAlignment> {
        if self.fusion_ids.is_empty() {
            return self.hmm_al_idx.map(|idx| &batch.hmm_als[idx]);
        }
        self.fusion_ids.iter().find_map(|&idx| {
            batch.blast_als[idx]
                .hmm_al_idx
                .map(|hmm_idx| &batch.hmm_als[hmm_idx])
        })
    }

    fn get_method(&self, batch: &Batch, cds: Option<&Locus>) -> String {
        let mut method = if self.from_hmm {
            "HMM".to_string()
        } else if self.is_mutation_prot() || self.is_strong_susceptible_prot(batch) {
            "POINT".to_string()
        } else if self.ref_prot_exactly_matched(true) {
            if self.allele() {
                "ALLELE".to_string()
            } else {
                "EXACT".to_string()
            }
        } else if self.partial() {
            if cds.is_some_and(|cds| self.truncated(cds)) {
                "PARTIAL_CONTIG_END".to_string()
            } else {
                "PARTIAL".to_string()
            }
        } else {
            "BLAST".to_string()
        };

        let mut suffix = true;
        if matches!(method.as_str(), "BLAST" | "PARTIAL" | "PARTIAL_CONTIG_END")
            && self.hsp.s_internal_stop
        {
            method = "INTERNAL_STOP".to_string();
            suffix = false;
        }
        if suffix && method != "HMM" {
            method.push(if self.hsp.s_prot { 'P' } else { 'X' });
        }
        method
    }

    fn pass_blast_rule(&self, br: &BlastRule) -> bool {
        !br.empty()
            && self.hsp.rel_identity() + FRAC_DELTA >= br.ident
            && self.hsp.q_rel_coverage() + FRAC_DELTA >= br.ref_coverage
    }

    fn finish(&mut self) -> Result<()> {
        self.partial_dna = false;
        if !self.hsp.s_prot && self.hsp.s_abs_coverage() >= 30 {
            const MISMATCH_TAIL_DNA: usize = 10;
            if self.hsp.q_int.start > 0 && self.hsp.s_tail(true) <= MISMATCH_TAIL_DNA {
                self.partial_dna = true;
            } else if self.hsp.q_int.stop < self.hsp.qlen
                && self.hsp.s_tail(false) <= MISMATCH_TAIL_DNA
            {
                self.partial_dna = true;
            }
        }

        if !self.hsp.s_prot {
            self.cdss.clear();
            self.cdss.push(Locus::new(
                0,
                &self.hsp.sseqid,
                self.hsp.s_int.start,
                self.hsp.s_int.stop,
                self.hsp.s_int.strand == 1,
                self.partial_dna,
                0,
                String::new(),
                String::new(),
            )?);
        }

        Ok(())
    }
}

fn parse_declarative_ref_mutation(
    ref_accession: &mut String,
    is_mutation_prot: bool,
) -> Result<AmrMutation> {
    let Some(gene_mutation_pos) = ref_accession.rfind(':') else {
        return Ok(AmrMutation::default());
    };
    if !is_mutation_prot {
        bail!("declarative mutation accession on non-mutation protein");
    }
    let gene_mutation = ref_accession[gene_mutation_pos + 1..].to_string();
    let accession_and_pos = &ref_accession[..gene_mutation_pos];
    let Some(pos_pos) = accession_and_pos.rfind(':') else {
        bail!("missing mutation position")
    };
    let pos = accession_and_pos[pos_pos + 1..].parse::<usize>()?;
    let accession = accession_and_pos[..pos_pos].to_string();
    *ref_accession = accession;
    if pos == 0 {
        bail!("mutation position must be positive");
    }
    Ok(AmrMutation::new(
        pos,
        &gene_mutation,
        &gene_mutation,
        "X",
        "X",
        "X",
    ))
}

// --- Batch ---

/// Core processing container for AMR report generation
pub struct Batch {
    pub fam_map: HashMap<String, Fam>,
    pub hmm2fam: HashMap<String, String>, // hmm_id -> fam_id
    pub reportable_min: u8,
    pub ident_min: f64,
    pub coverage_min: f64,
    pub report_all_equal: bool,
    pub print_node_raw: bool,
    pub name: String,
    pub cds_exist: bool,
    pub suppress_prots: Vec<String>,
    pub alien_prots: Vec<String>,
    pub dna_lens: HashMap<String, usize>,
    pub blast_als: Vec<BlastAlignment>,
    pub hmm_als: Vec<HmmAlignment>,
    pub hmm_exist: bool,
    pub domains: HashMap<(String, String), HmmDomain>, // (sseqid, fam_id) -> domain
    pub accession2mutations: HashMap<String, Vec<AmrMutation>>,
    pub accession2susceptible: HashMap<String, Susceptible>,
    // Computed by process()
    pub target2blast_als: BTreeMap<String, Vec<usize>>, // sseqid -> indices into blast_als
    pub target2good_blast_als: BTreeMap<String, Vec<usize>>,
    pub target2hmm_als: BTreeMap<String, Vec<usize>>,
    pub target2good_hmm_als: BTreeMap<String, Vec<usize>>,
}

impl Batch {
    /// Create a new Batch by loading the FAM file
    pub fn from_fam_file(fam_path: &Path, reportable_min: u8) -> Result<Self> {
        let mut fam_map: HashMap<String, Fam> = HashMap::new();
        let mut hmm2fam: HashMap<String, String> = HashMap::new();
        let default_complete_br = BlastRule::new(IDENT_MIN_DEF, COMPLETE_COVERAGE_MIN_DEF);
        let default_partial_br = BlastRule::new(IDENT_MIN_DEF, PARTIAL_COVERAGE_MIN_DEF);

        let file = File::open(fam_path)?;
        let reader = BufReader::new(file);

        // Pass 1: Create Fam objects
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 13 {
                bail!(
                    "Cannot read {}, line {}: expected at least 13 tab-delimited fields",
                    fam_path.display(),
                    line_no + 1
                );
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
            let tc1: f64 = fields[4].parse().with_context(|| {
                format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
            })?;
            let tc2: f64 = fields[5].parse().with_context(|| {
                format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
            })?;

            let mut complete_br = BlastRule {
                ident: fields[6].parse().with_context(|| {
                    format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
                })?,
                target_coverage: fields[7].parse().with_context(|| {
                    format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
                })?,
                ref_coverage: fields[8].parse().with_context(|| {
                    format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
                })?,
            };
            let mut partial_br = BlastRule {
                ident: fields[9].parse().with_context(|| {
                    format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
                })?,
                target_coverage: fields[10].parse().with_context(|| {
                    format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
                })?,
                ref_coverage: fields[11].parse().with_context(|| {
                    format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
                })?,
            };
            for value in [
                &mut complete_br.ident,
                &mut complete_br.target_coverage,
                &mut complete_br.ref_coverage,
                &mut partial_br.ident,
                &mut partial_br.target_coverage,
                &mut partial_br.ref_coverage,
            ] {
                if !(0.0..=100.0).contains(value) {
                    bail!(
                        "Cannot read {}, line {}: BLAST rule percentage must be between 0 and 100",
                        fam_path.display(),
                        line_no + 1
                    );
                }
                *value /= 100.0;
            }
            if complete_br.empty() != partial_br.empty() {
                bail!(
                    "Cannot read {}, line {}: complete and partial BLAST rules must both be empty or both be defined",
                    fam_path.display(),
                    line_no + 1
                );
            }
            if complete_br.empty() {
                complete_br = default_complete_br.clone();
                partial_br = default_partial_br.clone();
            } else {
                complete_br.ref_coverage = default_complete_br.ref_coverage;
                partial_br.ref_coverage = default_partial_br.ref_coverage;
            }

            let reportable: u8 = fields[12].parse().with_context(|| {
                format!("Cannot read {}, line {}", fam_path.display(), line_no + 1)
            })?;
            let type_ = fields.get(13).copied().unwrap_or("").to_string();
            let subtype = fields.get(14).copied().unwrap_or("").to_string();
            let class = fields.get(15).copied().unwrap_or("").to_string();
            let subclass = fields.get(16).copied().unwrap_or("").to_string();
            let family_name = if fields.len() > 17 {
                fields[17..].join("\t")
            } else {
                String::new()
            };
            let family_name = if family_name == "NULL" {
                String::new()
            } else {
                family_name
            };

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

            if fam_map.insert(fam_id.clone(), fam).is_some() {
                bail!("Family {fam_id} is duplicated");
            }
        }

        // Pass 2: validate the family tree parent links after all families exist.
        let parent_links: Vec<(String, String)> = fam_map
            .values()
            .filter(|fam| !fam.parent_id.is_empty())
            .map(|fam| (fam.id.clone(), fam.parent_id.clone()))
            .collect();
        for (child_id, parent_id) in parent_links {
            if !fam_map.contains_key(&parent_id) {
                bail!(
                    "parentFamId \"{}\" is not found in famId2fam for child \"{}\"",
                    parent_id,
                    child_id
                );
            }
        }

        Ok(Batch {
            fam_map,
            hmm2fam,
            reportable_min,
            ident_min: -1.0,
            coverage_min: 0.5,
            report_all_equal: false,
            print_node_raw: false,
            name: String::new(),
            cds_exist: false,
            suppress_prots: Vec::new(),
            alien_prots: Vec::new(),
            dna_lens: HashMap::new(),
            blast_als: Vec::new(),
            hmm_als: Vec::new(),
            hmm_exist: false,
            domains: HashMap::new(),
            accession2mutations: HashMap::new(),
            accession2susceptible: HashMap::new(),
            target2blast_als: BTreeMap::new(),
            target2good_blast_als: BTreeMap::new(),
            target2hmm_als: BTreeMap::new(),
            target2good_hmm_als: BTreeMap::new(),
        })
    }

    fn load_dna_lens(&mut self, path: &Path) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            let mut fields = line.split_whitespace();
            let Some(contig) = fields.next() else {
                continue;
            };
            let Some(len) = fields.next() else {
                continue;
            };
            let len = len.parse::<usize>()?;
            self.dna_lens.insert(contig.to_string(), len);
        }

        Ok(())
    }

    fn target_len(&self, al: &BlastAlignment) -> usize {
        if al.hsp.s_prot {
            al.hsp.slen
        } else {
            al.hsp.s_abs_coverage() / 3
        }
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

            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() < 8 {
                bail!("Reading {} line:\n{}", path.display(), line);
            }

            let organism_field = fields[0].replace('_', " ");
            let accession = fields[1].to_string();
            if organism_field != organism {
                self.alien_prots.push(accession);
                continue;
            }

            let pos: usize = fields[2].parse().unwrap_or(0);
            if pos == 0 {
                bail!("Mutation position must be positive: {line}");
            }
            let gene_mutation_std = fields[3];
            let gene_mutation = fields[4];
            let class = fields[5];
            let subclass = fields[6];
            let name = fields[7];

            let mutation =
                AmrMutation::new(pos, gene_mutation_std, gene_mutation, class, subclass, name);

            self.accession2mutations
                .entry(accession)
                .or_default()
                .push(mutation);
        }
        for (accession, mutations) in &mut self.accession2mutations {
            mutations.sort();
            if mutations.windows(2).any(|pair| pair[0] == pair[1]) {
                bail!("Duplicate mutations for {accession}");
            }
        }
        self.alien_prots.sort();
        self.alien_prots.dedup();

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

            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() < 7 {
                bail!("Reading {} line:\n{}", path.display(), line);
            }

            let organism_field = fields[0].replace('_', " ");
            let accession = fields[2].to_string();
            if organism_field != organism {
                self.alien_prots.push(accession);
                continue;
            }

            let cutoff = fields[3]
                .parse::<f64>()
                .with_context(|| format!("Reading {} line:\n{}", path.display(), line))?
                / 100.0;
            if !(0.0..1.0).contains(&cutoff) {
                bail!("Reading {} line:\n{}", path.display(), line);
            }

            let susceptible = Susceptible {
                genesymbol: fields[1].to_string(),
                cutoff,
                class: fields[4].to_string(),
                subclass: fields[5].to_string(),
                name: fields[6].replace('_', " "),
            };

            if self
                .accession2susceptible
                .insert(accession.clone(), susceptible)
                .is_some()
            {
                bail!(
                    "Duplicate protein accession {} in {}",
                    accession,
                    path.display()
                );
            }
        }
        self.alien_prots.sort();
        self.alien_prots.dedup();

        Ok(())
    }

    pub fn load_suppress(&mut self, path: &Path, _organism: &str) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            for accver in trimmed.split_whitespace() {
                self.suppress_prots.push(accver.to_string());
            }
        }
        self.suppress_prots.sort();
        self.suppress_prots.dedup();
        Ok(())
    }

    /// Add a BLAST alignment and index it
    pub fn add_blast_al(&mut self, al: BlastAlignment) {
        let idx = self.blast_als.len();
        let sseqid = al.hsp.sseqid.clone();
        self.blast_als.push(al);
        self.target2blast_als.entry(sseqid).or_default().push(idx);
    }

    /// Add an HMM alignment and index it
    pub fn add_hmm_al(&mut self, al: HmmAlignment) {
        let idx = self.hmm_als.len();
        let sseqid = al.sseqid.clone();
        self.hmm_als.push(al);
        self.target2hmm_als.entry(sseqid).or_default().push(idx);
    }

    fn annotate_protein_mutations(&mut self) {
        for al in &mut self.blast_als {
            if !al.is_mutation_prot() {
                continue;
            }
            let Some(mutations) = self.accession2mutations.get(&al.ref_accession) else {
                continue;
            };
            al.ref_mutations = mutations.clone();

            if !al.ref_mutation.empty() {
                if al.ref_mutation.frameshift != crate::seq::NO_INDEX && !al.hsp.blastx() {
                    continue;
                }
                let Some(mutation_idx) = mutations
                    .iter()
                    .position(|mutation| mutation == &al.ref_mutation)
                else {
                    continue;
                };
                al.ref_mutation = mutations[mutation_idx].clone();
                let start = if al.ref_mutation.pos_real < al.hsp.q_int.start
                    || al.ref_mutation.pos_real + al.ref_mutation.allele.len() > al.hsp.q_int.stop
                {
                    crate::seq::NO_INDEX
                } else {
                    let qseq = al.hsp.qseq.as_bytes();
                    let sseq = al.hsp.sseq.as_bytes();
                    let allele = al.ref_mutation.allele.as_bytes();
                    let mut pos = al.hsp.q_int.start;
                    let mut frameshift_i = crate::seq::NO_INDEX;
                    let mut i = 0usize;

                    let mut found = crate::seq::NO_INDEX;
                    while i < qseq.len() {
                        debug_assert!(pos <= al.ref_mutation.pos_real);
                        if qseq[i] != b'-' {
                            if pos == al.ref_mutation.pos_real {
                                let align_stop = i + allele.len();
                                if align_stop <= qseq.len()
                                    && &qseq[i..align_stop] == allele
                                    && &sseq[i..align_stop] == allele
                                    && (pos == 0 || (i > 0 && qseq[i - 1] == sseq[i - 1]))
                                {
                                    if al.ref_mutation.frameshift != crate::seq::NO_INDEX {
                                        frameshift_i = i;
                                        i += allele.len();
                                        pos = al.ref_mutation.get_stop();
                                        break;
                                    }
                                    let ref_stop = pos + allele.len();
                                    debug_assert!(align_stop <= qseq.len());
                                    if ref_stop == al.hsp.qlen
                                        || (align_stop < qseq.len()
                                            && qseq[align_stop] == sseq[align_stop])
                                    {
                                        found = i;
                                    }
                                }
                                break;
                            }
                            pos += 1;
                        }
                        i += 1;
                    }

                    while found == crate::seq::NO_INDEX
                        && frameshift_i != crate::seq::NO_INDEX
                        && i < qseq.len()
                    {
                        debug_assert_ne!(frameshift_i, crate::seq::NO_INDEX);
                        if qseq[i] != b'-' {
                            debug_assert!(al.ref_mutation.get_stop() > 0);
                            if pos == al.ref_mutation.get_stop() - 1 + al.ref_mutation.frameshift {
                                if qseq[i] == b'*' && sseq[i] == b'*' {
                                    found = frameshift_i;
                                }
                                break;
                            }
                            pos += 1;
                        }
                        i += 1;
                    }
                    found
                };
                if start != crate::seq::NO_INDEX {
                    let mut seq_change = SeqChange {
                        start,
                        len: al.ref_mutation.allele.len(),
                        reference: al.ref_mutation.reference.clone(),
                        allele: al.ref_mutation.allele.clone(),
                        ..SeqChange::default()
                    };
                    if seq_change.finish_pos(&al.hsp, 0) {
                        seq_change.mutations.push(mutation_idx);
                        al.seq_changes.push(seq_change);
                    }
                }
                al.seq_changes.sort();
                continue;
            }

            let qseq = al.hsp.qseq.as_bytes();
            let sseq = al.hsp.sseq.as_bytes();
            debug_assert_eq!(qseq.len(), sseq.len());

            let mut seq_change = SeqChange::default();
            let mut in_seq_change = false;
            for i in 0..=qseq.len() {
                if in_seq_change {
                    if i == qseq.len() || sseq[i] == qseq[i] {
                        if seq_change.finish(&al.hsp, &al.hsp.qseq, 0) {
                            al.seq_changes.push(seq_change);
                        }
                        seq_change = SeqChange::default();
                        in_seq_change = false;
                    } else {
                        seq_change.len += 1;
                    }
                } else if i < qseq.len() && sseq[i] != qseq[i] {
                    seq_change.start = i;
                    seq_change.len = 1;
                    in_seq_change = true;
                }
            }

            let mut mutation_idx = 0usize;
            while mutation_idx < al.ref_mutations.len()
                && al.ref_mutations[mutation_idx].pos_real < al.hsp.q_int.start
            {
                mutation_idx += 1;
            }

            let mut start_ref_prev = crate::seq::NO_INDEX;
            for seq_change in &mut al.seq_changes {
                seq_change.qc(&al.hsp, &al.ref_mutations);
                if start_ref_prev != crate::seq::NO_INDEX {
                    debug_assert!(start_ref_prev <= seq_change.start_ref);
                }
                while mutation_idx < al.ref_mutations.len() {
                    let mutation = &al.ref_mutations[mutation_idx];
                    if mutation.pos_real >= seq_change.stop_ref {
                        break;
                    }
                    if seq_change
                        .matches_mutation(mutation)
                        .expect("mutation reference sequence matches alignment")
                    {
                        seq_change.mutations.push(mutation_idx);
                    }
                    mutation_idx += 1;
                }
                start_ref_prev = seq_change.start_ref;
            }

            let mut ref_pos = al.hsp.q_int.start;
            let mut i = 0usize;
            for (mutation_idx, mutation) in al.ref_mutations.iter().enumerate() {
                while ref_pos < mutation.pos_real {
                    if i >= qseq.len() {
                        break;
                    }
                    debug_assert_ne!(qseq[i], b'-');
                    i += 1;
                    while i < qseq.len() && qseq[i] == b'-' {
                        i += 1;
                    }
                    if i >= qseq.len() {
                        break;
                    }
                    ref_pos += 1;
                }
                if i >= qseq.len() {
                    break;
                }
                if ref_pos == mutation.pos_real
                    && al.hsp.qseq[i..].starts_with(&mutation.reference)
                    && al.hsp.sseq[i..].starts_with(&mutation.reference)
                {
                    al.seq_changes.push(SeqChange {
                        mutations: vec![mutation_idx],
                        ..SeqChange::default()
                    });
                }
            }
            al.seq_changes.sort();
        }
    }

    fn process_disruptions(&mut self, orig_indices: &[usize]) -> Vec<usize> {
        let mut indices = Vec::new();
        let mut orig_hsps = Vec::new();
        let mut orig_al_indices = Vec::new();

        for &idx in orig_indices {
            let al = &self.blast_als[idx];
            if al.hsp.blastx() && !al.is_mutation_prot() && al.is_strong_susceptible_prot(self) {
                orig_hsps.push(al.hsp.clone());
                orig_al_indices.push(idx);
            } else {
                indices.push(idx);
            }
        }

        let mut merge = HspMerge::new(orig_hsps, None, 20.0, true);
        loop {
            let (hsp, orig_hsp, _score) = merge.get();
            if hsp.empty() {
                break;
            }
            let orig_hsp = orig_hsp.expect("merged BLASTX HSP has an original alignment");
            let mut al = self.blast_als[orig_al_indices[orig_hsp]].clone();
            al.hsp = hsp;
            al.ref_mutation = AmrMutation::default();
            al.seq_changes.clear();
            let idx = self.blast_als.len();
            self.blast_als.push(al);
            indices.push(idx);
        }

        indices
    }

    /// Process all alignments: filter, merge, resolve
    pub fn process(&mut self, retain_blasts: bool, skip_hmm_check: bool) {
        for target in self.target2blast_als.keys().cloned().collect::<Vec<_>>() {
            let Some(mut indices) = self.target2blast_als.remove(&target) else {
                continue;
            };
            indices.sort_by(|&a, &b| Hsp::less(&self.blast_als[a].hsp, &self.blast_als[b].hsp));

            let mut merged = Vec::new();
            let mut group = Vec::new();
            let mut prev: Option<usize> = None;
            for idx in indices {
                if let Some(prev_idx) = prev {
                    let al = &self.blast_als[idx];
                    let prev_al = &self.blast_als[prev_idx];
                    if al.hsp.sseqid != prev_al.hsp.sseqid
                        || al.hsp.s_int.strand != prev_al.hsp.s_int.strand
                        || al.hsp.qseqid != prev_al.hsp.qseqid
                    {
                        merged.extend(self.process_disruptions(&group));
                        group.clear();
                    }
                }
                group.push(idx);
                prev = Some(idx);
            }
            if !group.is_empty() {
                merged.extend(self.process_disruptions(&group));
            }

            for &idx in &merged {
                if self.blast_als[idx].is_strong_susceptible_prot(self) {
                    if let Some(disr) = self.blast_als[idx].hsp.disrs.first().cloned() {
                        self.blast_als[idx].seq_changes.push(SeqChange {
                            disr: Some(disr),
                            ..SeqChange::default()
                        });
                    }
                }
            }

            if !merged.is_empty() {
                self.target2blast_als.insert(target, merged);
            }
        }

        for indices in self.target2blast_als.values() {
            for &idx in indices {
                if !self.blast_als[idx].from_hmm && !self.blast_als[idx].hsp.s_prot {
                    self.blast_als[idx]
                        .finish()
                        .expect("BLASTX alignment has a valid target locus");
                }
            }
        }

        self.annotate_protein_mutations();

        // Set BlastRules from FAM hierarchy for each alignment
        for al in &mut self.blast_als {
            if !al.from_hmm {
                al.br_fam_id = None;
                al.br_fam_checked = false;
            }
            if al.is_susceptible_prot() {
                if let Some(susceptible) = self.accession2susceptible.get(&al.ref_accession) {
                    al.genesymbol = susceptible.genesymbol.clone();
                }
            }
            if !al.from_hmm && al.in_fam() && !al.fam_id.is_empty() {
                let mut fam = self.fam_map.get(&al.fam_id).or_else(|| {
                    if al.gene.is_empty() {
                        None
                    } else {
                        self.fam_map.get(&al.gene)
                    }
                });
                let ident = al.hsp.rel_identity();
                let ref_cov = al.hsp.q_rel_coverage();
                let partial = al.partial();
                while let Some(candidate) = fam {
                    let br = if partial {
                        &candidate.partial_br
                    } else {
                        &candidate.complete_br
                    };
                    al.br_fam_checked = true;
                    if ident + FRAC_DELTA >= br.ident && ref_cov + FRAC_DELTA >= br.ref_coverage {
                        al.br_fam_id = Some(candidate.id.clone());
                        if !candidate.genesymbol.is_empty() {
                            al.genesymbol = candidate.genesymbol.clone();
                        }
                        al.complete_br = candidate.complete_br.clone();
                        al.partial_br = candidate.partial_br.clone();
                        break;
                    }
                    fam = if candidate.parent_id.is_empty() {
                        None
                    } else {
                        self.fam_map.get(&candidate.parent_id)
                    };
                }
            }
        }

        // Step 1: Good BLAST filtering — keep only hits above thresholds
        for target in self.target2blast_als.keys().cloned().collect::<Vec<_>>() {
            let Some(indices) = self.target2blast_als.get_mut(&target) else {
                continue;
            };
            indices.retain(|&idx| {
                self.alien_prots
                    .binary_search(&self.blast_als[idx].ref_accession)
                    .is_err()
            });
            if indices.is_empty() {
                self.target2blast_als.remove(&target);
            }
        }

        let mut filtered_target2blast_als = BTreeMap::new();
        for (target, indices) in &self.target2blast_als {
            let good: Vec<usize> = indices
                .iter()
                .filter(|&&idx| {
                    let al = &self.blast_als[idx];
                    let susceptible = if al.is_susceptible_prot() {
                        self.accession2susceptible.get(&al.ref_accession)
                    } else {
                        None
                    };
                    al.from_hmm || al.good(susceptible)
                })
                .copied()
                .collect();
            if !good.is_empty() {
                filtered_target2blast_als.insert(target.clone(), good.clone());
            }
        }
        self.target2blast_als = filtered_target2blast_als;

        self.mark_protein_internal_stops_from_blastx();

        // Step 2: Good HMM filtering
        for (target, indices) in &self.target2hmm_als {
            let good: Vec<usize> = indices
                .iter()
                .filter(|&&idx| self.hmm_als[idx].good(&self.fam_map))
                .copied()
                .collect();
            if !good.is_empty() {
                let filtered = self.pareto_filter_hmm(&good);
                self.target2good_hmm_als.insert(target.clone(), filtered);
            }
        }

        if retain_blasts {
            self.target2good_blast_als = self.target2blast_als.clone();
        } else {
            // Step 3: Pareto-better BLAST filtering
            // For each target, keep only alignments that are not dominated by another
            let targets: Vec<String> = self.target2blast_als.keys().cloned().collect();
            for target in &targets {
                if let Some(indices) = self.target2blast_als.get(target) {
                    let filtered = self.pareto_filter_blast(indices);
                    self.target2good_blast_als.insert(target.clone(), filtered);
                }
            }
        }

        // Cf. dna_mutation.cpp: mark weaker colocated mutation calls as replaced.
        for indices in self.target2good_blast_als.values() {
            let mut replacements = Vec::new();
            for &idx1 in indices {
                let al1 = &self.blast_als[idx1];
                for (seq_change1_idx, seq_change1) in al1.seq_changes.iter().enumerate() {
                    for &idx2 in indices {
                        if idx1 == idx2 {
                            continue;
                        }
                        let al2 = &self.blast_als[idx2];
                        if al2.hsp.s_prot != al1.hsp.s_prot
                            || al2.hsp.sseqid != al1.hsp.sseqid
                            || al2.hsp.s_int.strand != al1.hsp.s_int.strand
                        {
                            continue;
                        }
                        for (seq_change2_idx, seq_change2) in al2.seq_changes.iter().enumerate() {
                            if seq_change1.start_target == seq_change2.start_target
                                && seq_change1.has_frameshift(&al1.ref_mutations)
                                    == seq_change2.has_frameshift(&al2.ref_mutations)
                                && seq_change1.better(
                                    seq_change2,
                                    al1.hsp.rel_identity(),
                                    al2.hsp.rel_identity(),
                                )
                            {
                                replacements.push((idx2, seq_change2_idx, seq_change1_idx));
                            }
                        }
                    }
                }
            }
            for (idx, seq_change_idx, replacement_idx) in replacements {
                self.blast_als[idx].seq_changes[seq_change_idx].replacement = Some(replacement_idx);
            }
        }

        if !skip_hmm_check {
            self.remove_blast_hits_without_required_hmm();
        }

        // C++ PD-2783: first let stronger BLAST rows consume HMM rows.
        let hmm_targets: Vec<String> = self.target2good_hmm_als.keys().cloned().collect();
        for target in &hmm_targets {
            let Some(hmm_indices) = self.target2good_hmm_als.get(target).cloned() else {
                continue;
            };
            let blast_indices = self
                .target2good_blast_als
                .get(target)
                .cloned()
                .unwrap_or_default();
            let mut remaining = Vec::new();
            for hmm_idx in hmm_indices {
                let suppressing_blast = blast_indices.iter().copied().find(|&blast_idx| {
                    let blast_al = &self.blast_als[blast_idx];
                    blast_al.in_fam()
                        && blast_al.better_hmm(
                            &self.hmm_als[hmm_idx],
                            &self.fam_map,
                            &self.blast_als,
                        )
                });
                if let Some(blast_idx) = suppressing_blast {
                    self.blast_als[blast_idx].hmm_al_idx = Some(hmm_idx);
                } else {
                    remaining.push(hmm_idx);
                }
            }
            if remaining.is_empty() {
                self.target2good_hmm_als.remove(target);
            } else {
                self.target2good_hmm_als.insert(target.clone(), remaining);
            }
        }

        // Then let more specific HMM rows remove weaker BLAST rows.
        let blast_targets: Vec<String> = self.target2good_blast_als.keys().cloned().collect();
        for target in blast_targets {
            let Some(indices) = self.target2good_blast_als.get(&target).cloned() else {
                continue;
            };
            let hmm_indices = self
                .target2good_hmm_als
                .get(&target)
                .cloned()
                .unwrap_or_default();
            let filtered: Vec<usize> = indices
                .into_iter()
                .filter(|&idx| {
                    let al = &self.blast_als[idx];
                    !al.in_fam()
                        || !hmm_indices
                            .iter()
                            .any(|&hmm_idx| self.hmm_als[hmm_idx].better_blast(al, &self.fam_map))
                })
                .collect();
            if filtered.is_empty() {
                self.target2good_blast_als.remove(&target);
            } else {
                self.target2good_blast_als.insert(target, filtered);
            }
        }

        // Surviving HMM rows are reported through their synthetic BlastAlignment.
        let hmm_targets: Vec<String> = self.target2good_hmm_als.keys().cloned().collect();
        for target in hmm_targets {
            let Some(hmm_indices) = self.target2good_hmm_als.get(&target).cloned() else {
                continue;
            };
            for hmm_idx in hmm_indices {
                let blast_idx = match self.hmm_als[hmm_idx].blast_al_idx {
                    Some(idx) if self.blast_als[idx].from_hmm => idx,
                    Some(best_idx) => {
                        let best = self.blast_als[best_idx].clone();
                        let fam = self.fam_map.get(&self.hmm_als[hmm_idx].fam_id);
                        let Some(domain) = self.hmm_als[hmm_idx].domain.as_ref() else {
                            continue;
                        };
                        let mut al = BlastAlignment::from_hmm_alignment(
                            &self.hmm_als[hmm_idx],
                            fam,
                            Some(&best),
                            domain,
                        );
                        al.hmm_al_idx = Some(hmm_idx);
                        let new_idx = self.blast_als.len();
                        self.blast_als.push(al);
                        self.hmm_als[hmm_idx].blast_al_idx = Some(new_idx);
                        new_idx
                    }
                    None => {
                        let fam = self.fam_map.get(&self.hmm_als[hmm_idx].fam_id);
                        let Some(domain) = self.hmm_als[hmm_idx].domain.as_ref() else {
                            continue;
                        };
                        let mut al = BlastAlignment::from_hmm_alignment(
                            &self.hmm_als[hmm_idx],
                            fam,
                            None,
                            domain,
                        );
                        al.hmm_al_idx = Some(hmm_idx);
                        let new_idx = self.blast_als.len();
                        self.blast_als.push(al);
                        self.hmm_als[hmm_idx].blast_al_idx = Some(new_idx);
                        new_idx
                    }
                };
                let entry = self
                    .target2good_blast_als
                    .entry(target.clone())
                    .or_default();
                if !entry.contains(&blast_idx) {
                    entry.push(blast_idx);
                }
            }
        }

        // C++ Batch::process final pass: BlastAlignment::{fusions,fusionRedundant}.
        for indices in self.target2good_blast_als.values_mut() {
            indices.sort_by(|&a_idx, &b_idx| {
                let a = &self.blast_als[a_idx];
                let b = &self.blast_als[b_idx];
                if a.less(b) {
                    std::cmp::Ordering::Less
                } else if b.less(a) {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Equal
                }
            });

            let mut prev_idx: Option<usize> = None;
            let mut fusion_main_idx: Option<usize> = None;
            for &idx in indices.iter() {
                let al = &self.blast_als[idx];
                if !al.in_fam() || al.from_hmm {
                    continue;
                }

                if let Some(prev) = prev_idx {
                    if self.blast_als[prev].same_match(al) {
                        debug_assert_eq!(al.parts, self.blast_als[prev].parts);
                        debug_assert!(al.part >= self.blast_als[prev].part);
                        if al.part == self.blast_als[prev].part {
                            fusion_main_idx = None;
                        } else {
                            let main = if let Some(main) = fusion_main_idx {
                                main
                            } else {
                                fusion_main_idx = Some(prev);
                                if self.blast_als[prev].fusion_ids.is_empty() {
                                    self.blast_als[prev].fusion_ids.push(prev);
                                }
                                prev
                            };
                            debug_assert!(self.blast_als[idx].parts >= 2);
                            self.blast_als[main].fusion_ids.push(idx);
                            self.blast_als[idx].fusion_redundant = true;
                        }
                    } else {
                        fusion_main_idx = None;
                    }
                }
                prev_idx = Some(idx);
            }
        }
    }

    fn remove_blast_hits_without_required_hmm(&mut self) {
        if !self.hmm_exist {
            return;
        }

        let targets: Vec<String> = self.target2good_blast_als.keys().cloned().collect();
        for target in targets {
            let Some(indices) = self.target2good_blast_als.get(&target) else {
                continue;
            };
            let filtered: Vec<usize> = indices
                .iter()
                .copied()
                .filter(|&idx| {
                    let al = &self.blast_als[idx];
                    if !al.in_fam()
                        || !al.hsp.s_prot
                        || al.partial()
                        || al.hsp.rel_identity() >= 0.98 - FRAC_DELTA
                    {
                        return true;
                    }

                    let Some(hmm_fam) = al
                        .get_match_fam(&self.fam_map)
                        .and_then(|fam| fam.get_hmm_fam(&self.fam_map))
                    else {
                        return true;
                    };

                    self.target2good_hmm_als
                        .get(&target)
                        .is_some_and(|hmm_indices| {
                            hmm_indices.iter().any(|&hmm_idx| {
                                let hmm_al = &self.hmm_als[hmm_idx];
                                hmm_al.sseqid == al.hsp.sseqid && hmm_al.fam_id == hmm_fam.id
                            })
                        })
                })
                .collect();
            if filtered.is_empty() {
                self.target2good_blast_als.remove(&target);
            } else {
                self.target2good_blast_als.insert(target, filtered);
            }
        }
    }

    fn mark_protein_internal_stops_from_blastx(&mut self) {
        let protein_indices: Vec<usize> = self
            .target2blast_als
            .values()
            .flatten()
            .copied()
            .chain(self.hmm_als.iter().filter_map(|hmm_al| hmm_al.blast_al_idx))
            .filter(|&idx| {
                let al = &self.blast_als[idx];
                al.hsp.s_prot && !al.cdss.is_empty()
            })
            .collect();

        for idx in protein_indices {
            for cds in &self.blast_als[idx].cdss {
                let Some(blastx_indices) = self.target2blast_als.get(&cds.contig) else {
                    continue;
                };
                let has_internal_stop = blastx_indices.iter().any(|&blastx_idx| {
                    let blastx_al = &self.blast_als[blastx_idx];
                    !blastx_al.hsp.s_prot
                        && blastx_al.hsp.s_internal_stop
                        && self.blast_better(&self.blast_als[idx], blastx_al)
                });
                if has_internal_stop {
                    self.blast_als[idx].hsp.s_internal_stop = true;
                    break;
                }
            }
        }
    }

    /// Get the reportability level for an alignment
    fn get_reportable(&self, al: &BlastAlignment) -> u8 {
        let fam = if al.from_hmm {
            self.fam_map.get(&al.fam_id).or_else(|| {
                if !al.gene.is_empty() {
                    self.fam_map.get(&al.gene)
                } else {
                    None
                }
            })
        } else {
            al.br_fam_id
                .as_ref()
                .and_then(|fam_id| self.fam_map.get(fam_id))
        };
        if let Some(fam) = fam {
            return fam.reportable;
        }
        0
    }

    fn fusion_2_reportable(&self, al: &BlastAlignment) -> u8 {
        if al.is_susceptible_prot()
            && (al.susceptible_idx.is_some()
                || self.accession2susceptible.contains_key(&al.ref_accession))
        {
            return 2;
        }
        if al.fusion_ids.is_empty() {
            return self.get_reportable(al);
        }
        al.fusion_ids
            .iter()
            .map(|&idx| self.get_reportable(&self.blast_als[idx]))
            .max()
            .unwrap_or(0)
    }

    fn get_gene_symbol(&self, al: &BlastAlignment) -> String {
        if al.is_susceptible_prot() {
            if let Some(susceptible) = self.accession2susceptible.get(&al.ref_accession) {
                return susceptible.genesymbol.clone();
            }
        }
        if al.allele_match() {
            return al.fam_id.clone();
        }
        al.get_match_fam(&self.fam_map)
            .map(|fam| fam.genesymbol.as_str())
            .filter(|genesymbol| !genesymbol.is_empty())
            .unwrap_or(crate::columns::NA)
            .to_string()
    }

    fn fusion_2_gene_symbols(&self, al: &BlastAlignment) -> String {
        if al.fusion_ids.is_empty() {
            return self.get_gene_symbol(al);
        }
        al.fusion_ids
            .iter()
            .map(|&idx| self.get_gene_symbol(&self.blast_als[idx]))
            .collect::<Vec<_>>()
            .join("/")
    }

    fn is_core(&self, al: &BlastAlignment) -> bool {
        self.fusion_2_reportable(al) >= 2 || (al.allele_match() && al.reportable >= 2)
    }

    fn fusion_2_core(&self, al: &BlastAlignment) -> bool {
        if al.fusion_ids.is_empty() {
            return self.is_core(al);
        }
        al.fusion_ids
            .iter()
            .any(|&idx| self.is_core(&self.blast_als[idx]))
    }

    fn get_type(&self, al: &BlastAlignment) -> String {
        if al.is_susceptible_prot()
            && (al.susceptible_idx.is_some()
                || self.accession2susceptible.contains_key(&al.ref_accession))
        {
            "AMR".to_string()
        } else {
            let (_, type_, _, _, _, _) = self.get_fam_info(al);
            type_
        }
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

    fn pareto_filter_hmm(&self, indices: &[usize]) -> Vec<usize> {
        let mut current = indices.to_vec();
        let mut good: Vec<usize> = Vec::new();

        for criterion in 0..2 {
            good.clear();
            for &idx in &current {
                let hmm_al = &self.hmm_als[idx];
                if good.iter().any(|&good_idx| {
                    self.hmm_als[good_idx].better(
                        hmm_al,
                        criterion,
                        &self.fam_map,
                        self.report_all_equal,
                    )
                }) {
                    continue;
                }
                good.retain(|&good_idx| {
                    !hmm_al.better(
                        &self.hmm_als[good_idx],
                        criterion,
                        &self.fam_map,
                        self.report_all_equal,
                    )
                });
                good.push(idx);
            }
            if criterion < 1 {
                current = good.clone();
            }
        }

        good
    }

    /// C++ BlastAlignment::better(): is `a` strictly better than `b`?
    fn blast_better(&self, a: &BlastAlignment, b: &BlastAlignment) -> bool {
        self.blast_better_eq(a, b)
            && (!self.blast_better_eq(b, a)
                || (!a.is_mutation_prot()
                    && !self.report_all_equal
                    && a.ref_accession < b.ref_accession))
    }

    /// C++ BlastAlignment::betterEq(), bounded to currently translated fields.
    fn blast_better_eq(&self, a: &BlastAlignment, b: &BlastAlignment) -> bool {
        if std::ptr::eq(a, b) {
            return true;
        }
        if a.is_mutation_prot() != b.is_mutation_prot() {
            return false;
        }
        if a.ref_mutation.empty() != b.ref_mutation.empty() {
            return false;
        }
        if a.is_susceptible_prot() != b.is_susceptible_prot() {
            return false;
        }

        if a.hsp.s_prot == b.hsp.s_prot {
            if a.hsp.sseqid != b.hsp.sseqid {
                return false;
            }
            if !b.inside_eq(a, self) && !a.inside_eq(b, self) {
                if !a.hsp.s_prot
                    || (a.in_fam()
                        && b.in_fam()
                        && self.fusion_2_gene_symbols(a) != self.fusion_2_gene_symbols(b))
                {
                    return false;
                }
            }

            if a.in_fam() && !a.ref_accession.is_empty() && a.ref_accession == b.ref_accession {
                if b.hsp.nident < a.hsp.nident {
                    return true;
                }
                if a.hsp.nident < b.hsp.nident {
                    return false;
                }
                if a.hsp.s_int.start < b.hsp.s_int.start {
                    return true;
                }
                if b.hsp.s_int.start < a.hsp.s_int.start {
                    return false;
                }
                if b.hsp.s_int.stop < a.hsp.s_int.stop {
                    return true;
                }
                if a.hsp.s_int.stop < b.hsp.s_int.stop {
                    return false;
                }
                if a.hsp.q_int.start < b.hsp.q_int.start {
                    return true;
                }
                if b.hsp.q_int.start < a.hsp.q_int.start {
                    return false;
                }
                if b.hsp.q_int.stop < a.hsp.q_int.stop {
                    return true;
                }
                if a.hsp.q_int.stop < b.hsp.q_int.stop {
                    return false;
                }
            } else {
                if b.has_mutation() != a.has_mutation() {
                    return !b.has_mutation() && a.has_mutation();
                }
                if a.fusion_overrides(b, self) {
                    return true;
                }
                if b.fusion_overrides(a, self) {
                    return false;
                }
                let a_exact = a.ref_prot_exactly_matched(false);
                let b_exact = b.ref_prot_exactly_matched(false);
                if b_exact != a_exact {
                    return !b_exact && a_exact;
                }
                if b.hsp.nident < a.hsp.nident {
                    return true;
                }
                if a.hsp.nident < b.hsp.nident {
                    return false;
                }
                if a.ref_effective_len() < b.ref_effective_len() {
                    return true;
                }
                if b.ref_effective_len() < a.ref_effective_len() {
                    return false;
                }
            }
        } else {
            if a.hsp.s_prot && !a.matches_cds(b) {
                return false;
            }
            if !a.hsp.s_prot && !b.matches_cds(a) {
                return false;
            }
            let a_exact = a.ref_prot_exactly_matched(false);
            let b_exact = b.ref_prot_exactly_matched(false);
            if b_exact != a_exact {
                return !b_exact && a_exact;
            }
            if b.allele_match() != a.allele_match() {
                return !b.allele_match() && a.allele_match();
            }
        }

        if a.is_mutation_prot() && b.is_mutation_prot() {
            let a_mutation_symbols = a.get_mutation_symbols();
            let b_mutation_symbols = b.get_mutation_symbols();
            if a_mutation_symbols == b_mutation_symbols && a.hsp.s_prot != b.hsp.s_prot {
                return a.hsp.s_prot;
            }
            if a_mutation_symbols.is_superset(&b_mutation_symbols)
                && !b_mutation_symbols.is_superset(&a_mutation_symbols)
            {
                return true;
            }
            if b_mutation_symbols.is_superset(&a_mutation_symbols)
                && !a_mutation_symbols.is_superset(&b_mutation_symbols)
            {
                return false;
            }
        }

        if b.hsp.s_prot != a.hsp.s_prot {
            return a.hsp.s_prot;
        }
        true
    }

    fn alignment_exact(&self, al: &BlastAlignment) -> bool {
        al.ref_prot_exactly_matched(true)
    }

    fn allele_match(&self, al: &BlastAlignment) -> bool {
        al.allele_match()
    }

    /// Generate TSV report to output
    pub fn report(&self, out: &mut dyn Write, print_node: bool) -> Result<()> {
        let mut tsv = TsvOut::new(Some(out));
        tsv.use_pound = false;

        self.write_header(&mut tsv, print_node)?;

        // Output rows for each target
        for (target_name, indices) in &self.target2good_blast_als {
            for &idx in indices {
                let al = &self.blast_als[idx];
                // Skip mutation proteins without detected mutations
                if al.is_mutation_prot() && al.seq_changes.is_empty() {
                    continue;
                }
                let ordinary_reportable = self.fusion_2_reportable(al) >= self.reportable_min
                    || (al.allele_match() && al.reportable >= 2);
                if al.seq_changes.is_empty()
                    && (!ordinary_reportable
                        || self.suppress_prots.contains(&al.ref_accession)
                        || al.fusion_redundant)
                {
                    continue;
                }
                self.report_alignment(al, target_name, &mut tsv, print_node, false)?;
            }
        }

        Ok(())
    }

    pub fn report_mutation_all(&self, out: &mut dyn Write, print_node: bool) -> Result<()> {
        let mut tsv = TsvOut::new(Some(out));
        tsv.use_pound = false;
        self.write_header(&mut tsv, print_node)?;

        for (target_name, indices) in &self.target2good_blast_als {
            for &idx in indices {
                let al = &self.blast_als[idx];
                if !al.is_mutation_prot() {
                    continue;
                }
                let ordinary_reportable = self.fusion_2_reportable(al) >= self.reportable_min
                    || (al.allele_match() && al.reportable >= 2);
                if al.seq_changes.is_empty()
                    && (!ordinary_reportable
                        || self.suppress_prots.contains(&al.ref_accession)
                        || al.fusion_redundant)
                {
                    continue;
                }
                self.report_alignment(al, target_name, &mut tsv, print_node, true)?;
            }
        }

        Ok(())
    }

    fn write_header(&self, tsv: &mut TsvOut, print_node: bool) -> Result<()> {
        if !self.name.is_empty() {
            tsv.write_field(&"Name")?;
        }
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
        tsv.new_ln()?;
        Ok(())
    }

    /// Get fam-derived info for an alignment
    fn get_fam_info(&self, al: &BlastAlignment) -> (String, String, String, String, String, u8) {
        // Returns: (genesymbol, type, subtype, class, subclass, reportable)
        // C++ getFam(): try famId first, then gene field as fallback
        let fam = self.fam_map.get(&al.fam_id).or_else(|| {
            if !al.gene.is_empty() {
                self.fam_map.get(&al.gene)
            } else {
                None
            }
        });
        let match_fam = if al.from_hmm {
            fam
        } else {
            al.br_fam_id
                .as_ref()
                .and_then(|fam_id| self.fam_map.get(fam_id))
        };

        let is_exact = self.alignment_exact(al);
        let allele_match = self.allele_match(al);

        let genesymbol = if al.is_mutation_prot() {
            if !al.ref_mutation.gene_mutation.is_empty() {
                al.ref_mutation.gene_mutation.clone()
            } else {
                al.genesymbol.clone()
            }
        } else if allele_match {
            al.fam_id.clone()
        } else if is_exact || al.parts >= 2 {
            // For exact/allele matches: use the fam's own genesymbol
            if !self.fam_map.contains_key(&al.fam_id) {
                al.fam_id.clone()
            } else if let Some(f) = fam {
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
        let (type_, subtype, mut class, mut subclass) =
            if al.is_mutation_prot() && !al.ref_mutation.empty() {
                (
                    "AMR".to_string(),
                    "POINT".to_string(),
                    al.ref_mutation.class.clone(),
                    al.ref_mutation.subclass.clone(),
                )
            } else if let Some(mf) = match_fam {
                (
                    mf.type_.clone(),
                    mf.subtype.clone(),
                    mf.class.clone(),
                    mf.subclass.clone(),
                )
            } else if let Some(f) = fam {
                (
                    f.type_.clone(),
                    f.subtype.clone(),
                    f.class.clone(),
                    f.subclass.clone(),
                )
            } else {
                (
                    al.resistance.clone(),
                    String::new(),
                    al.class.clone(),
                    al.subclass.clone(),
                )
            };
        if allele_match && !al.subclass.is_empty() {
            class = al.class.clone();
            subclass = al.subclass.clone();
        }

        let reportable = self.get_reportable(al);

        (genesymbol, type_, subtype, class, subclass, reportable)
    }

    fn report_alignment(
        &self,
        al: &BlastAlignment,
        target_name: &str,
        tsv: &mut TsvOut,
        print_node: bool,
        mutation_all: bool,
    ) -> Result<()> {
        let na = crate::columns::NA;
        let (mut genesymbol, mut type_, mut subtype, mut class, mut subclass, mut reportable) =
            self.get_fam_info(al);

        let is_mutation = al.is_mutation_prot();
        if !is_mutation {
            genesymbol = self.fusion_2_gene_symbols(al);
            reportable = self.fusion_2_reportable(al);
            type_ = al.fusion_2_type(self);
            subtype = al.fusion_2_subtype(self);
            class = al.fusion_2_class(self);
            subclass = al.fusion_2_subclass(self);
        }

        let is_exact = self.alignment_exact(al);

        let empty_seq_change = SeqChange::default();
        let seq_changes: Vec<&SeqChange> = if al.seq_changes.is_empty() {
            vec![&empty_seq_change]
        } else {
            al.seq_changes.iter().collect()
        };
        let cds_count = al.cdss.len().max(1);
        for cds_idx in 0..cds_count {
            let cds = al.cdss.get(cds_idx);
            if al.hsp.s_prot && target_name != al.hsp.sseqid {
                if cds.is_none() {
                    continue;
                } else if cds.is_some_and(|cds| cds.contig != target_name) {
                    continue;
                }
            } else if al.hsp.sseqid != target_name {
                continue;
            }
            let method = al.get_method(self, cds);
            for seq_change in &seq_changes {
                let mutation_indices: Vec<Option<usize>> = if seq_change.mutations.is_empty() {
                    vec![None]
                } else {
                    seq_change.mutations.iter().copied().map(Some).collect()
                };
                for mutation_idx in mutation_indices {
                    let mutation = mutation_idx.and_then(|idx| al.ref_mutations.get(idx));
                    let is_mutation_row = !seq_change.empty()
                        && mutation.is_some()
                        && seq_change.replacement.is_none();
                    if mutation_all {
                        if !is_mutation {
                            continue;
                        }
                    } else if is_mutation && !is_mutation_row {
                        continue;
                    }

                    if !self.name.is_empty() {
                        tsv.write_field(&self.name)?;
                    }
                    // Product name: susceptible rows use organism-specific metadata; exact/allele uses product.
                    let protein_name = if is_mutation {
                        al.product.clone()
                    } else if al.is_susceptible_prot() {
                        if let Some(susceptible) = self.accession2susceptible.get(&al.ref_accession)
                        {
                            susceptible.name.clone()
                        } else {
                            al.product.clone()
                        }
                    } else if al.ref_accession.is_empty() {
                        al.product.clone()
                    } else if is_exact || al.parts >= 2 {
                        al.product.clone()
                    } else if let Some(mf) = if al.from_hmm {
                        self.fam_map.get(&al.fam_id).or_else(|| {
                            if !al.gene.is_empty() {
                                self.fam_map.get(&al.gene)
                            } else {
                                None
                            }
                        })
                    } else {
                        al.br_fam_id
                            .as_ref()
                            .and_then(|fam_id| self.fam_map.get(fam_id))
                    } {
                        if !mf.family_name.is_empty() {
                            mf.family_name.clone()
                        } else {
                            na.to_string()
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
                        if let Some(cds) = cds {
                            tsv.write_field(&cds.contig)?;
                            tsv.write_field(&(cds.start + 1))?;
                            tsv.write_field(&cds.stop)?;
                            tsv.write_field(&if cds.strand { "+" } else { "-" })?;
                        } else {
                            tsv.write_field(&al.hsp.sseqid)?;
                            tsv.write_field(&(al.hsp.s_int.start + 1))?;
                            tsv.write_field(&al.hsp.s_int.stop)?;
                            tsv.write_field(
                                &crate::seq::strand2char(al.hsp.s_int.strand).to_string(),
                            )?;
                        }
                    }

                    // Element symbol
                    let element_symbol = if is_mutation {
                        if let Some(mutation) = mutation {
                            if seq_change.empty() {
                                mutation.wildtype()
                            } else {
                                mutation.gene_mutation.clone()
                            }
                        } else {
                            format!("{}_{}", al.gene, seq_change.get_mutation_str())
                        }
                    } else if seq_change.empty() {
                        genesymbol.clone()
                    } else {
                        format!(
                            "{}{}{}",
                            al.gene,
                            crate::columns::DISRUPTION_DELIM,
                            seq_change.get_mutation_str()
                        )
                    };
                    tsv.write_field(&if element_symbol.is_empty() {
                        na.to_string()
                    } else {
                        element_symbol
                    })?;
                    // Element name
                    let element_name = if is_mutation {
                        if let Some(mutation) = mutation {
                            if seq_change.empty() {
                                format!("{} [WILDTYPE]", protein_name)
                            } else {
                                mutation.name.clone()
                            }
                        } else {
                            format!("{} [UNKNOWN]", protein_name)
                        }
                    } else {
                        protein_name
                    };
                    tsv.write_field(&if element_name.is_empty() {
                        na.to_string()
                    } else {
                        element_name
                    })?;

                    // Scope
                    let scope = if is_mutation
                        || reportable >= 2
                        || (!is_mutation && self.fusion_2_core(al))
                    {
                        "core"
                    } else {
                        "plus"
                    };
                    tsv.write_field(&scope)?;

                    // Type, Subtype, Class, Subclass
                    if is_mutation {
                        tsv.write_field(&"AMR")?;
                        tsv.write_field(&"POINT")?;
                        if let Some(mutation) = mutation {
                            tsv.write_field(&if mutation.class.is_empty() {
                                na
                            } else {
                                &mutation.class
                            })?;
                            tsv.write_field(&if mutation.subclass.is_empty() {
                                na
                            } else {
                                &mutation.subclass
                            })?;
                        } else {
                            tsv.write_field(&na)?;
                            tsv.write_field(&na)?;
                        }
                    } else {
                        tsv.write_field(&if type_.is_empty() {
                            na.to_string()
                        } else {
                            type_.clone()
                        })?;
                        tsv.write_field(&if al.is_strong_susceptible_prot(self) {
                            "POINT_DISRUPT".to_string()
                        } else if subtype.is_empty() {
                            na.to_string()
                        } else {
                            subtype.clone()
                        })?;
                        tsv.write_field(&if class.is_empty() {
                            na.to_string()
                        } else {
                            class.clone()
                        })?;
                        tsv.write_field(&if subclass.is_empty() {
                            na.to_string()
                        } else {
                            subclass.clone()
                        })?;
                    }

                    // Method
                    tsv.write_field(&method)?;

                    // Target length
                    tsv.write_field(&self.target_len(al))?;

                    if al.ref_accession.is_empty() {
                        tsv.write_field(&na)?;
                        tsv.write_field(&na)?;
                        tsv.write_field(&na)?;
                        tsv.write_field(&na)?;
                        tsv.write_field(&na)?;
                        tsv.write_field(&na)?;
                    } else {
                        tsv.write_field(&al.hsp.q_effective_len())?;
                        tsv.write_field(&format!("{:.2}", al.hsp.q_rel_coverage() * 100.0))?;
                        tsv.write_field(&format!("{:.2}", al.hsp.rel_identity() * 100.0))?;
                        tsv.write_field(&al.hsp.sseq.len())?;
                        tsv.write_field(&al.ref_accession)?;
                        tsv.write_field(&al.product)?;
                    }

                    // HMM info — C++ uses getFam()->hmm which walks famId then gene fallback
                    // Try: 1) from_hmm fam, 2) hmm_al_idx, 3) fam hierarchy (famId→gene), 4) NA
                    if is_mutation {
                        tsv.write_field(&na)?;
                        tsv.write_field(&na)?;
                    } else {
                        let hmm_fam = if let Some(hmm_al) = al.fusion_2_hmm_al(self) {
                            self.fam_map.get(&hmm_al.fam_id)
                        } else if al.from_hmm {
                            self.fam_map.get(&al.fam_id)
                        } else if let Some(hmm_idx) = al.hmm_al_idx {
                            self.fam_map.get(&self.hmm_als[hmm_idx].fam_id)
                        } else {
                            None
                        };
                        if let Some(fam) = hmm_fam.filter(|f| !f.hmm.is_empty()) {
                            tsv.write_field(&fam.hmm)?;
                            tsv.write_field(&fam.family_name)?;
                        } else {
                            tsv.write_field(&na)?;
                            tsv.write_field(&na)?;
                        }
                    }

                    if print_node {
                        if !self.print_node_raw && al.allele() && !al.ref_prot_exactly_matched(true)
                        {
                            tsv.write_field(&al.gene)?;
                        } else {
                            tsv.write_field(&al.fusion_2_fam_ids(&self.blast_als))?;
                        }
                    }

                    tsv.new_ln()?;
                }
            }
        }

        Ok(())
    }
}

/// Parse BLASTP/BLASTX tabular output and add to batch
pub fn load_blast_results(
    blast_file: &Path,
    batch: &mut Batch,
    is_protein: bool, // true for blastp, false for blastx
    nosame: bool,
) -> Result<()> {
    let file = File::open(blast_file)?;
    let reader = BufReader::new(file);

    let ident_min = if batch.ident_min >= 0.0 {
        batch.ident_min
    } else {
        IDENT_MIN_DEF
    };
    let partial_coverage_min = if batch.coverage_min >= 0.0 {
        batch.coverage_min
    } else {
        PARTIAL_COVERAGE_MIN_DEF
    };
    let default_complete_br = BlastRule::new(ident_min, COMPLETE_COVERAGE_MIN_DEF);
    let default_partial_br = BlastRule::new(ident_min, partial_coverage_min);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let mut al = BlastAlignment::from_blast_line(
            trimmed,
            true,       // q_prot (reference is always protein)
            is_protein, // s_prot (target is protein for blastp, DNA for blastx)
            &default_complete_br,
            &default_partial_br,
        )?;

        al.br_fam_checked = false;
        if al.is_susceptible_prot() {
            if let Some(susceptible) = batch.accession2susceptible.get(&al.ref_accession) {
                al.genesymbol = susceptible.genesymbol.clone();
            }
        }
        if al.in_fam() {
            let mut fam = batch.fam_map.get(&al.fam_id).or_else(|| {
                if al.gene.is_empty() {
                    None
                } else {
                    batch.fam_map.get(&al.gene)
                }
            });
            if fam.is_none() {
                bail!(
                    "Cannot find hierarchy for: {} (genesymbol: {})",
                    al.fam_id,
                    al.gene
                );
            }
            let ident = al.hsp.rel_identity();
            let ref_cov = al.hsp.q_rel_coverage();
            let partial = al.partial();
            while let Some(candidate) = fam {
                if !candidate.complete_br.empty() {
                    al.br_fam_checked = true;
                    let br = if partial {
                        &candidate.partial_br
                    } else {
                        &candidate.complete_br
                    };
                    if ident + FRAC_DELTA >= br.ident && ref_cov + FRAC_DELTA >= br.ref_coverage {
                        al.br_fam_id = Some(candidate.id.clone());
                        al.complete_br = candidate.complete_br.clone();
                        al.partial_br = candidate.partial_br.clone();
                        break;
                    }
                }
                fam = if candidate.parent_id.is_empty() {
                    None
                } else {
                    batch.fam_map.get(&candidate.parent_id)
                };
            }
        }

        if nosame && al.ref_accession == al.hsp.sseqid {
            continue;
        }

        batch.add_blast_al(al);
    }

    Ok(())
}

/// Parse HMM domain table output
pub fn load_hmm_domains(dom_file: &Path, batch: &mut Batch) -> Result<()> {
    batch.hmm_exist = true;

    let file = File::open(dom_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 23 {
            continue;
        }

        // --domtblout format:
        // 0:target_name 1:accession 2:tlen 3:query_name 4:query_accession 5:qlen
        // 6:E-value 7:score 8:bias 9:# 10:of 11:c-evalue 12:i-evalue 13:score 14:bias
        // 15:hmm_from 16:hmm_to 17:ali_from 18:ali_to 19:env_from 20:env_to 21:acc
        let target_name = fields[0].to_string();
        let query_accession = fields[4].to_string();
        let seq_len: usize = fields[2].parse().unwrap_or(0);
        let hmm_len: usize = fields[5].parse().unwrap_or(0);
        let score: f64 = fields[13].parse().unwrap_or(0.0);
        let hmm_start: usize = fields[15].parse::<usize>().unwrap_or(1).saturating_sub(1);
        let hmm_stop: usize = fields[16].parse().unwrap_or(0);
        let seq_start: usize = fields[17].parse::<usize>().unwrap_or(1).saturating_sub(1);
        let seq_stop: usize = fields[18].parse().unwrap_or(0);

        if score <= 0.0 || hmm_start >= hmm_stop || seq_start >= seq_stop {
            continue;
        }

        // Look up fam by HMM accession
        let fam_id = match batch.hmm2fam.get(&query_accession) {
            Some(id) => id.clone(),
            None => continue,
        };

        let domain = HmmDomain {
            score,
            hmm_len,
            hmm_start,
            hmm_stop,
            seq_len,
            seq_start,
            seq_stop,
        };

        let key = (target_name, fam_id);
        let existing = batch.domains.get(&key);
        if existing.is_none() || existing.unwrap().score <= score {
            batch.domains.insert(key, domain);
        }
    }

    Ok(())
}

/// Parse HMM search results
pub fn load_hmm_results(hmmsearch_file: &Path, batch: &mut Batch) -> Result<()> {
    let file = File::open(hmmsearch_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 9 {
            continue;
        }

        // hmmsearch --tblout format:
        // target_name  accession  query_name  query_accession  E-value  score  bias  E-value  score  ...
        let sseqid = fields[0].to_string();
        let hmm_accession = fields[3].to_string();
        let score1: f64 = fields[5].parse().unwrap_or(0.0);
        let score2: f64 = fields[8].parse().unwrap_or(0.0);

        // Look up fam
        let fam_id = match batch.hmm2fam.get(&hmm_accession) {
            Some(id) => id.clone(),
            None => bail!("No family for HMM {hmm_accession}"),
        };

        let hmm_al = HmmAlignment {
            sseqid: sseqid.clone(),
            score1,
            score2,
            fam_id: fam_id.clone(),
            domain: None,
            blast_al_idx: None,
        };

        // Check if good
        if !hmm_al.good(&batch.fam_map) {
            continue;
        }

        // Get domain info
        let domain = batch
            .domains
            .get(&(sseqid.clone(), fam_id.clone()))
            .cloned();

        if domain.is_none() || domain.as_ref().unwrap().hmm_len == 0 {
            continue;
        }

        // Find best BLAST alignment for this target
        let mut best_blast_idx = None;
        let mut nident_max = 0;
        if let Some(indices) = batch.target2blast_als.get(&sseqid) {
            for &idx in indices {
                let nident = batch.blast_als[idx].hsp.nident;
                if nident > nident_max {
                    nident_max = nident;
                    best_blast_idx = Some(idx);
                }
            }
        }

        // Create synthetic BlastAlignment from HMM
        let domain = domain.unwrap();
        // Update domain info on HMM alignment
        let mut hmm_al_with_domain = hmm_al;
        hmm_al_with_domain.domain = Some(domain.clone());
        hmm_al_with_domain.blast_al_idx = best_blast_idx;
        let hmm_idx = batch.hmm_als.len();
        batch.add_hmm_al(hmm_al_with_domain);

        if let Some(blast_idx) = best_blast_idx {
            let susceptible = if batch.blast_als[blast_idx].is_susceptible_prot() {
                batch
                    .accession2susceptible
                    .get(&batch.blast_als[blast_idx].ref_accession)
            } else {
                None
            };
            if batch.blast_als[blast_idx].good(susceptible) {
                let replace = batch.blast_als[blast_idx]
                    .hmm_al_idx
                    .map(|old_idx| batch.hmm_als[hmm_idx].score1 > batch.hmm_als[old_idx].score1)
                    .unwrap_or(true);
                if replace {
                    batch.blast_als[blast_idx].hmm_al_idx = Some(hmm_idx);
                }
                continue;
            }

            let fam = batch.fam_map.get(&fam_id);
            let mut al = batch.blast_als[blast_idx].clone();
            al.from_hmm = true;
            al.hmm_al_idx = Some(hmm_idx);
            al.fam_id = fam_id.clone();
            al.gene = fam_id.clone();
            al.resistance.clear();
            al.class = fam.map(|f| f.class.clone()).unwrap_or_default();
            al.subclass = fam.map(|f| f.subclass.clone()).unwrap_or_default();
            al.product = fam.map(|f| f.family_name.clone()).unwrap_or_default();
            al.reportable = fam.map(|f| f.reportable).unwrap_or(0);
            al.genesymbol = fam.map(|f| f.genesymbol.clone()).unwrap_or_default();
            batch.add_blast_al(al);
            let synthetic_idx = batch.blast_als.len() - 1;
            batch.hmm_als[hmm_idx].blast_al_idx = Some(synthetic_idx);
            continue;
        }

        let hmm_blast_al = {
            // Stand-alone HMM hit — create minimal BlastAlignment
            let fam = batch.fam_map.get(&fam_id);
            let genesymbol = fam.map(|f| f.genesymbol.clone()).unwrap_or_default();
            let family_name = fam.map(|f| f.family_name.clone()).unwrap_or_default();

            let align_len = domain.hmm_stop - domain.hmm_start;
            let mut hsp = Hsp::default();
            hsp.sseqid = sseqid.clone();
            hsp.qseqid = fam_id.clone();
            hsp.slen = domain.seq_len;
            hsp.s_int = crate::seq::Interval::new(domain.seq_start, domain.seq_stop, 0);
            hsp.qlen = domain.hmm_len;
            hsp.q_int = crate::seq::Interval::new(domain.hmm_start, domain.hmm_stop, 0);
            hsp.q_prot = true;
            hsp.s_prot = true;
            hsp.a_prot = true;
            hsp.length = align_len;
            hsp.nident = (align_len as f64 * 0.7) as usize;

            BlastAlignment {
                hsp,
                partial_dna: false,
                from_hmm: true,
                ref_accession: String::new(),
                part: 1,
                parts: 1,
                fam_id: fam_id.clone(),
                gene: fam_id.clone(),
                resistance: String::new(),
                class: fam.map(|f| f.class.clone()).unwrap_or_default(),
                subclass: fam.map(|f| f.subclass.clone()).unwrap_or_default(),
                product: family_name,
                reportable: fam.map(|f| f.reportable).unwrap_or(0),
                genesymbol,
                method: "HMM".to_string(),
                complete_br: BlastRule::default(),
                partial_br: BlastRule::default(),
                br_fam_id: None,
                br_fam_checked: false,
                cdss: Vec::new(),
                hmm_al_idx: Some(hmm_idx),
                susceptible_idx: None,
                seq_changes: Vec::new(),
                ref_mutations: Vec::new(),
                ref_mutation: crate::alignment::AmrMutation::default(),
                fusion_ids: Vec::new(),
                fusion_redundant: false,
            }
        };
        batch.add_blast_al(hmm_blast_al);
        let synthetic_idx = batch.blast_als.len() - 1;
        batch.hmm_als[hmm_idx].blast_al_idx = Some(synthetic_idx);
    }

    Ok(())
}

/// Configuration for amr_report
pub struct AmrReportConfig<'a> {
    pub fam_file: &'a Path,
    pub blastp_file: Option<&'a Path>,
    pub blastx_file: Option<&'a Path>,
    pub dna_len_file: Option<&'a Path>,
    pub hmmsearch_file: Option<&'a Path>,
    pub hmmdom_file: Option<&'a Path>,
    pub gff_file: Option<&'a Path>,
    pub gff_type: &'a str,
    pub organism: &'a str,
    pub mutation_file: Option<&'a Path>,
    pub susceptible_file: Option<&'a Path>,
    pub suppress_file: Option<&'a Path>,
    pub coverage_min: f64,
    pub ident_min: f64,
    pub print_node: bool,
    pub print_node_raw: bool,
    pub mutation_all: Option<&'a Path>,
    pub name: &'a str,
    pub non_reportable: bool,
    pub report_core_only: bool,
    pub report_all_equal: bool,
    pub cds_exist: bool,
    pub nosame: bool,
    pub noblast: bool,
    pub nohmm: bool,
    pub retain_blasts: bool,
    pub skip_hmm_check: bool,
}

/// Run the amr_report processing pipeline
pub fn run_amr_report(config: &AmrReportConfig, out: &mut dyn Write) -> Result<()> {
    if config.hmmsearch_file.is_some() != config.hmmdom_file.is_some() {
        bail!("hmmsearch and hmmdom must be both present or both absent");
    }
    if config.print_node_raw && !config.print_node {
        bail!("print_node_raw requires print_node");
    }
    if !config.organism.is_empty() && config.mutation_file.is_none() {
        bail!("mutation_tab is empty");
    }

    let ident_min_user = config.ident_min != -1.0;
    let ident_min = if ident_min_user {
        config.ident_min
    } else {
        IDENT_MIN_DEF
    };
    if !(0.0..=1.0).contains(&ident_min) {
        bail!("-ident_min must be -1 or between 0 and 1");
    }
    if !(0.0..=1.0).contains(&config.coverage_min) {
        bail!("-coverage_min must be -1 or between 0 and 1");
    }
    if config.coverage_min > COMPLETE_COVERAGE_MIN_DEF {
        bail!(
            "-coverage_min must be less than {} - threshod for complete matches",
            COMPLETE_COVERAGE_MIN_DEF
        );
    }
    if config.blastp_file.is_some() && config.blastx_file.is_some() && config.gff_file.is_none() {
        bail!("If BLASTP and BLASTX files are present then a GFF file must be present");
    }

    let reportable_min = if config.non_reportable {
        0
    } else if config.report_core_only {
        2
    } else {
        1
    };
    let mut batch = Batch::from_fam_file(config.fam_file, reportable_min)?;
    batch.cds_exist = config.cds_exist || config.blastx_file.is_some() || config.gff_file.is_some();
    batch.name = config.name.to_string();
    batch.ident_min = config.ident_min;
    batch.coverage_min = config.coverage_min;
    batch.report_all_equal = config.report_all_equal;
    batch.print_node_raw = config.print_node_raw;

    let partial_coverage_min = config.coverage_min;
    for fam in batch.fam_map.values_mut() {
        if fam.complete_br.empty() {
            continue;
        }
        fam.complete_br.ref_coverage = COMPLETE_COVERAGE_MIN_DEF;
        fam.partial_br.ref_coverage = partial_coverage_min;
        if ident_min_user {
            fam.complete_br.ident = ident_min;
            fam.partial_br.ident = ident_min;
        }
    }

    if let Some(dna_len_file) = config.dna_len_file {
        if dna_len_file.exists() {
            batch.load_dna_lens(dna_len_file)?;
        }
    }

    // Load mutations and susceptible data
    if let Some(mut_file) = config.mutation_file {
        let organism_display = config.organism.replace('_', " ");
        batch.load_mutations(mut_file, &organism_display)?;
    }
    if let Some(sus_file) = config.susceptible_file {
        let organism_display = config.organism.replace('_', " ");
        batch.load_susceptible(sus_file, &organism_display)?;
    }
    if let Some(suppress_file) = config.suppress_file {
        batch.load_suppress(suppress_file, config.organism)?;
    }

    // Load BLAST results
    if !config.noblast {
        if let Some(bp_file) = config.blastp_file {
            if bp_file.exists() {
                load_blast_results(bp_file, &mut batch, true, config.nosame)?;
            }
        }
    }

    // Load HMM results
    if !config.nohmm {
        if let (Some(dom_file), Some(search_file)) = (config.hmmdom_file, config.hmmsearch_file) {
            if dom_file.exists() && search_file.exists() {
                load_hmm_domains(dom_file, &mut batch)?;
                load_hmm_results(search_file, &mut batch)?;
            }
        }
    }

    // Load BLASTX results
    if !config.noblast {
        if let Some(bx_file) = config.blastx_file {
            if bx_file.exists() {
                load_blast_results(bx_file, &mut batch, false, config.nosame)?;
            }
        }
    }

    // Load GFF annotations and assign CDSs to BLAST alignments
    if let Some(gff_file) = config.gff_file {
        let gff_type = GffType::name2type(config.gff_type)?;
        let annot = Annot::from_gff(gff_file.to_str().unwrap_or(""), gff_type, false, false)?;
        // Assign CDS loci to protein BLAST alignments
        for al in &mut batch.blast_als {
            if !al.hsp.s_prot {
                continue;
            }
            if let Ok(loci) = annot.find_loci(&al.hsp.sseqid) {
                al.cdss = loci
                    .iter()
                    .map(|l| crate::gff::Locus {
                        line_num: l.line_num,
                        contig: l.contig.clone(),
                        start: l.start,
                        stop: l.stop,
                        strand: l.strand,
                        partial: l.partial,
                        contig_len: if l.contig_len == 0 {
                            batch.dna_lens.get(&l.contig).copied().unwrap_or(0)
                        } else {
                            l.contig_len
                        },
                        cross_origin: l.cross_origin,
                        gene: l.gene.clone(),
                        product: l.product.clone(),
                    })
                    .collect();
            }
        }
    }
    if config.blastx_file.is_some() {
        batch.target2blast_als.clear();
        for (idx, al) in batch.blast_als.iter().enumerate() {
            if al.hsp.s_prot {
                for contig in al.get_contigs() {
                    batch.target2blast_als.entry(contig).or_default().push(idx);
                }
            } else {
                batch
                    .target2blast_als
                    .entry(al.hsp.sseqid.clone())
                    .or_default()
                    .push(idx);
            }
        }

        batch.target2hmm_als.clear();
        for (idx, hmm_al) in batch.hmm_als.iter().enumerate() {
            let Some(blast_idx) = hmm_al.blast_al_idx else {
                continue;
            };
            for contig in batch.blast_als[blast_idx].get_contigs() {
                batch.target2hmm_als.entry(contig).or_default().push(idx);
            }
        }
    }

    // Process
    batch.process(config.retain_blasts, config.skip_hmm_check);

    // Report
    batch.report(out, config.print_node)?;

    if let Some(path) = config.mutation_all {
        let mut mutation_all = File::create(path)?;
        batch.report_mutation_all(&mut mutation_all, config.print_node)?;
    }

    Ok(())
}
