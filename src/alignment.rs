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

    /// Parse gene_mutation_std into components
    fn parse(gene_mutation_std: &str) -> (String, String, String, i32, usize, i32) {
        // Format: GENE_REFposALLELE or GENE_REFposALLELE*frameshift_pos
        // Simple parsing — the C++ version is more complex
        let mut reference = String::new();
        let mut allele = String::new();
        let mut gene = String::new();
        let mut pos_std: i32 = 0;
        let mut frameshift: usize = NO_INDEX;
        let frameshift_insertion: i32 = 0;

        if let Some(underscore_pos) = gene_mutation_std.find(PM_DELIMITER) {
            gene = gene_mutation_std[..underscore_pos].to_string();
            let rest = &gene_mutation_std[underscore_pos + 1..];

            // Find where the reference sequence ends and position begins
            let mut ref_end = 0;
            for (i, c) in rest.chars().enumerate() {
                if c.is_ascii_digit() || c == '-' {
                    ref_end = i;
                    break;
                }
            }

            reference = rest[..ref_end].to_string();

            // Find position
            let mut pos_end = ref_end;
            for (i, c) in rest[ref_end..].chars().enumerate() {
                if !c.is_ascii_digit() && c != '-' {
                    pos_end = ref_end + i;
                    break;
                }
                pos_end = ref_end + i + 1;
            }

            if let Ok(p) = rest[ref_end..pos_end].parse::<i32>() {
                pos_std = p - 1; // Convert from 1-based
            }

            // Remaining is the allele (possibly with frameshift info)
            let allele_str = &rest[pos_end..];
            if let Some(star_pos) = allele_str.find('*') {
                allele = allele_str[..star_pos].to_string();
                if let Ok(fs) = allele_str[star_pos + 1..].parse::<usize>() {
                    frameshift = fs;
                }
            } else {
                allele = allele_str.to_string();
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
    fn test_amr_mutation_parse() {
        let m = AmrMutation::new(100, "gyrA_S83L", "gyrA_S83L", "QUINOLONE", "FLUOROQUINOLONE", "Escherichia coli gyrA S83L");
        assert_eq!(m.gene, "gyrA");
        assert_eq!(m.reference, "S");
        assert_eq!(m.allele, "L");
        assert_eq!(m.pos_std, 82); // 0-based
        assert!(!m.empty());
    }

    #[test]
    fn test_amr_mutation_empty() {
        let m = AmrMutation::default();
        assert!(m.empty());
    }
}
