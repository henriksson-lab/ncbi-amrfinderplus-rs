// Protein or DNA mutations library — port of alignment.hpp/cpp

use crate::seq::{Disruption, Hsp, NO_INDEX};

/// Delimiter for point mutations
pub const PM_DELIMITER: char = '_';

/// Database point mutation
#[derive(Debug, Clone)]
pub struct AmrMutation {
    pub pos_real: usize,
    pub gene_mutation_std: String,
    // Parsed from gene_mutation_std
    pub reference: String,
    pub allele: String,
    pub gene: String,
    pub pos_std: i32,
    pub frameshift: usize,      // NO_INDEX if not a frameshift
    pub frameshift_insertion: i32,
    // To be reported
    pub gene_mutation: String,
    pub class: String,
    pub subclass: String,
    pub name: String,
}

impl AmrMutation {
    pub fn new(
        pos_real: usize,
        gene_mutation_std: &str,
        gene_mutation: &str,
        class: &str,
        subclass: &str,
        name: &str,
    ) -> Self {
        let (reference, allele, gene, pos_std, frameshift, frameshift_insertion) =
            Self::parse(gene_mutation_std);

        AmrMutation {
            pos_real,
            gene_mutation_std: gene_mutation_std.to_string(),
            reference,
            allele,
            gene,
            pos_std,
            frameshift,
            frameshift_insertion,
            gene_mutation: gene_mutation.to_string(),
            class: class.to_string(),
            subclass: subclass.to_string(),
            name: name.to_string(),
        }
    }

    pub fn empty(&self) -> bool {
        self.gene_mutation_std.is_empty()
    }

    pub fn get_stop(&self) -> usize {
        self.pos_real + self.reference.len()
    }

    pub fn wildtype(&self) -> String {
        format!(
            "{}_{}{}{}",
            self.gene, self.reference, self.pos_std + 1, self.reference
        )
    }

    /// Parse gene_mutation_std into components.
    ///
    /// Format: `GENE_REF<pos>ALLELE[*frameshift_pos]`
    /// Examples:
    /// - `gyrA_S83L` → gene=gyrA, ref=S, pos=83, allele=L
    /// - `blaTEMp_G162T` → gene=blaTEMp, ref=G, pos=162, allele=T
    /// - `ampC_T-14TGT` → gene=ampC, ref=T, pos=-14, allele=TGT
    /// - `nfsA_K141Ter` → gene=nfsA, ref=K, pos=141, allele=Ter (stop codon)
    /// - `nfsA_R15C` → gene=nfsA, ref=R, pos=15, allele=C
    fn parse(gene_mutation_std: &str) -> (String, String, String, i32, usize, i32) {
        let mut reference = String::new();
        let mut allele = String::new();
        let mut gene = String::new();
        let mut pos_std: i32 = 0;
        let mut frameshift: usize = NO_INDEX;
        let mut frameshift_insertion: i32 = 0;

        if let Some(underscore_pos) = gene_mutation_std.find(PM_DELIMITER) {
            gene = gene_mutation_std[..underscore_pos].to_string();
            let rest = &gene_mutation_std[underscore_pos + 1..];
            let chars: Vec<char> = rest.chars().collect();

            // Phase 1: Collect reference residues (uppercase letters at start)
            let mut i = 0;
            while i < chars.len() && chars[i].is_ascii_uppercase() {
                reference.push(chars[i]);
                i += 1;
            }

            // Phase 2: Collect position (digits and optional leading '-')
            let pos_start = i;
            if i < chars.len() && chars[i] == '-' {
                i += 1; // negative position
            }
            while i < chars.len() && chars[i].is_ascii_digit() {
                i += 1;
            }
            if i > pos_start {
                if let Ok(p) = rest[pos_start..i].parse::<i32>() {
                    pos_std = p - 1; // Convert from 1-based to 0-based (matching C++)
                }
            }

            // Phase 3: Collect allele (remaining, possibly with * frameshift)
            let allele_part = &rest[i..];
            if let Some(star_pos) = allele_part.find('*') {
                allele = allele_part[..star_pos].to_string();
                // After * is frameshift position, possibly with +/- insertion count
                let fs_part = &allele_part[star_pos + 1..];
                if let Some(sign_pos) = fs_part.find(['+', '-']) {
                    if let Ok(fs) = fs_part[..sign_pos].parse::<usize>() {
                        frameshift = fs;
                    }
                    if let Ok(ins) = fs_part[sign_pos..].parse::<i32>() {
                        frameshift_insertion = ins;
                    }
                } else if let Ok(fs) = fs_part.parse::<usize>() {
                    frameshift = fs;
                }
            } else {
                allele = allele_part.to_string();
            }
        }

        (reference, allele, gene, pos_std, frameshift, frameshift_insertion)
    }

    pub fn apply(&self, seq: &mut String) -> anyhow::Result<()> {
        if self.pos_real >= seq.len() {
            anyhow::bail!(
                "AmrMutation position {} is outside the sequence",
                self.pos_real
            );
        }
        if self.frameshift != NO_INDEX {
            anyhow::bail!("AmrMutation is a frameshift");
        }
        let prefix = &seq[..self.pos_real];
        let suffix = &seq[self.pos_real + self.reference.len()..];
        *seq = format!("{}{}{}", prefix, self.allele, suffix);
        Ok(())
    }
}

impl Default for AmrMutation {
    fn default() -> Self {
        AmrMutation {
            pos_real: 0,
            gene_mutation_std: String::new(),
            reference: String::new(),
            allele: String::new(),
            gene: String::new(),
            pos_std: 0,
            frameshift: NO_INDEX,
            frameshift_insertion: 0,
            gene_mutation: String::new(),
            class: String::new(),
            subclass: String::new(),
            name: String::new(),
        }
    }
}

/// Observed sequence change relative to reference
#[derive(Debug, Clone)]
pub struct SeqChange {
    pub start: usize,
    pub len: usize,
    pub reference: String,
    pub allele: String,
    pub start_ref: usize,
    pub stop_ref: usize,
    pub start_target: usize,
    pub neighborhood_mismatch: f64,
    pub mutations: Vec<usize>, // indices into AmrMutation array
    pub disr: Option<Disruption>,
    pub replacement: Option<usize>, // index to replacing SeqChange
}

impl SeqChange {
    pub fn empty(&self) -> bool {
        self.len == 0 && self.disr.is_none()
    }

    pub fn has_mutation(&self) -> bool {
        !self.empty() && !self.mutations.is_empty() && self.replacement.is_none()
    }

    pub fn has_frameshift(&self, mutations: &[AmrMutation]) -> bool {
        self.has_mutation() && mutations[self.mutations[0]].frameshift != NO_INDEX
    }

    pub fn is_frameshift(&self) -> bool {
        self.reference.is_empty()
    }

    pub fn get_stop(&self) -> usize {
        self.start + self.len
    }
}

impl Default for SeqChange {
    fn default() -> Self {
        SeqChange {
            start: 0,
            len: 0,
            reference: String::new(),
            allele: String::new(),
            start_ref: 0,
            stop_ref: 0,
            start_target: 0,
            neighborhood_mismatch: 0.0,
            mutations: Vec::new(),
            disr: None,
            replacement: None,
        }
    }
}

/// Alignment extends HSP with mutation analysis
#[derive(Debug, Clone)]
pub struct Alignment {
    pub hsp: Hsp,
    pub ref_mutation: AmrMutation,
    pub seq_changes: Vec<SeqChange>,
}

impl Alignment {
    pub fn from_blast_line(
        line: &str,
        q_prot: bool,
        s_prot: bool,
    ) -> anyhow::Result<Self> {
        let a_prot = q_prot || s_prot;
        let hsp = Hsp::from_blast_line(line, q_prot, s_prot, a_prot, q_prot, true)?;

        Ok(Alignment {
            hsp,
            ref_mutation: AmrMutation::default(),
            seq_changes: Vec::new(),
        })
    }

    pub fn has_mutation(&self) -> bool {
        self.seq_changes.iter().any(|sc| sc.has_mutation())
    }

    pub fn has_declarative_frameshift(&self, mutations: &[AmrMutation]) -> bool {
        self.seq_changes.len() == 1 && self.seq_changes[0].has_frameshift(mutations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_amr_mutation_parse_simple() {
        let m = AmrMutation::new(100, "gyrA_S83L", "gyrA_S83L", "QUINOLONE", "FLUOROQUINOLONE", "name");
        assert_eq!(m.gene, "gyrA");
        assert_eq!(m.reference, "S");
        assert_eq!(m.allele, "L");
        assert_eq!(m.pos_std, 82); // 0-based: 83 - 1
        assert_eq!(m.wildtype(), "gyrA_S83S"); // wildtype uses pos_std + 1
    }

    #[test]
    fn test_amr_mutation_parse_negative_pos() {
        let m = AmrMutation::new(40, "ampC_T-14TGT", "ampC_T-14TGT", "BETA-LACTAM", "CEPH", "name");
        assert_eq!(m.gene, "ampC");
        assert_eq!(m.reference, "T");
        assert_eq!(m.allele, "TGT");
        assert_eq!(m.pos_std, -15); // 0-based: -14 - 1
    }

    #[test]
    fn test_amr_mutation_parse_stop_codon() {
        let m = AmrMutation::new(141, "nfsA_K141Ter", "nfsA_K141Ter", "NITRO", "NITRO", "name");
        assert_eq!(m.gene, "nfsA");
        assert_eq!(m.reference, "K");
        assert_eq!(m.allele, "Ter");
        assert_eq!(m.pos_std, 140); // 0-based: 141 - 1
    }

    #[test]
    fn test_amr_mutation_parse_promoter() {
        let m = AmrMutation::new(162, "blaTEMp_G162T", "blaTEMp_G162T", "BL", "BL", "name");
        assert_eq!(m.gene, "blaTEMp");
        assert_eq!(m.reference, "G");
        assert_eq!(m.allele, "T");
        assert_eq!(m.pos_std, 161); // 0-based: 162 - 1
    }

    #[test]
    fn test_amr_mutation_empty() {
        let m = AmrMutation::default();
        assert!(m.empty());
    }
}
