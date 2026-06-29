// Protein or DNA mutations library — port of alignment.hpp/cpp

use crate::seq::{Disruption, Hsp, NO_INDEX, TERMINATOR_WORD};
use std::cmp::Ordering;
use std::io::Write;

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
    pub frameshift: usize, // NO_INDEX if not a frameshift
    pub frameshift_insertion: i32,
    // To be reported
    pub gene_mutation: String,
    pub class: String,
    pub subclass: String,
    pub name: String,
}

impl AmrMutation {
    pub fn new(
        pos_real_arg: usize,
        gene_mutation_std: &str,
        gene_mutation: &str,
        class: &str,
        subclass: &str,
        name: &str,
    ) -> Self {
        assert!(pos_real_arg > 0);
        let pos_real = pos_real_arg - 1;

        let (reference, mut allele, gene, pos_std, frameshift, frameshift_insertion) =
            Self::parse(gene_mutation_std);

        assert!(!name.is_empty());
        assert!(!name.contains('\t'));
        let name = name.replace('_', " ");
        assert!(!name.contains("  "));

        if allele == TERMINATOR_WORD {
            allele = "*".to_string();
        } else if allele == "del" {
            allele.clear();
        }

        let mutation = AmrMutation {
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
            name,
        };
        mutation.qc();
        mutation
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
            self.gene,
            self.reference,
            self.pos_std + 1,
            self.reference
        )
    }

    pub fn qc(&self) {
        if self.empty() {
            return;
        }
        assert!(!self.gene_mutation.is_empty());
        assert!(self
            .reference
            .chars()
            .all(|c| !c.is_ascii_alphabetic() || c.is_ascii_uppercase()));
        assert!(self
            .allele
            .chars()
            .all(|c| !c.is_ascii_alphabetic() || c.is_ascii_uppercase()));
        assert_ne!(self.reference, self.allele);
        assert!(!self.reference.contains('-'));
        assert!(!self.allele.contains('-'));
        assert!(!self.gene.is_empty());
        assert!(!self.reference.is_empty());
        if self.frameshift != NO_INDEX {
            assert!(!self.allele.is_empty());
            assert!(self.pos_std >= 0);
        }
        assert_eq!(self.frameshift == NO_INDEX, self.frameshift_insertion == 0);
        if self.frameshift_insertion != 0 {
            assert_ne!(self.frameshift_insertion % 3, 0);
        }
    }

    pub fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        if self.empty() {
            write!(os, "empty")?;
        } else {
            write!(
                os,
                "{} {} {} {}",
                self.pos_real + 1,
                self.gene_mutation,
                self.frameshift_insertion,
                self.name
            )?;
        }
        Ok(())
    }

    fn parse(gene_mutation_std: &str) -> (String, String, String, i32, usize, i32) {
        assert!(!gene_mutation_std.is_empty());

        let mut reference = String::new();
        let mut allele = String::new();
        let mut frameshift: usize = NO_INDEX;
        let mut frameshift_insertion: i32 = 0;
        let mut pure_mutation = gene_mutation_std.to_string();

        let fs_infix = "fsTer";
        if let Some(fs_infix_pos) = gene_mutation_std.find(fs_infix) {
            pure_mutation.truncate(fs_infix_pos);
            let mut suffix = gene_mutation_std[fs_infix_pos + fs_infix.len()..].to_string();
            let (indel_pos, ins) = if let Some(pos) = suffix.find("ins") {
                (pos, true)
            } else {
                (
                    suffix.find("del").expect("frameshift without ins/del"),
                    false,
                )
            };
            frameshift_insertion = suffix[indel_pos + 3..]
                .parse::<i32>()
                .expect("bad frameshift insertion length");
            if !ins {
                frameshift_insertion *= -1;
            }
            assert_ne!(frameshift_insertion % 3, 0);
            suffix.truncate(indel_pos);
            frameshift = suffix.parse::<usize>().expect("bad frameshift stop");
            assert_ne!(frameshift, NO_INDEX);
        }

        enum Type {
            InAllele,
            InPos,
            InRef,
        }

        let mut type_ = Type::InAllele;
        let mut pos_std_s = String::new();
        for (i, c) in pure_mutation.char_indices().rev() {
            match type_ {
                Type::InAllele => {
                    if c.is_ascii_alphabetic() {
                        allele.push(c);
                    } else {
                        pos_std_s.push(c);
                        type_ = Type::InPos;
                    }
                }
                Type::InPos => {
                    if c.is_ascii_alphabetic() {
                        reference.push(c);
                        type_ = Type::InRef;
                    } else {
                        pos_std_s.push(c);
                    }
                }
                Type::InRef => {
                    if c.is_ascii_alphabetic() {
                        reference.push(c);
                    } else {
                        assert_eq!(c, PM_DELIMITER);
                        assert!(i > 0);
                        allele = allele.chars().rev().collect();
                        reference = reference.chars().rev().collect();
                        pos_std_s = pos_std_s.chars().rev().collect();
                        let gene = pure_mutation[..i].to_string();
                        let pos_std = pos_std_s.parse::<i32>().expect("bad mutation position") - 1;
                        return (
                            reference,
                            allele,
                            gene,
                            pos_std,
                            frameshift,
                            frameshift_insertion,
                        );
                    }
                }
            }
        }
        panic!("bad mutation: {gene_mutation_std}");
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

impl PartialEq for AmrMutation {
    fn eq(&self, other: &Self) -> bool {
        self.gene_mutation_std == other.gene_mutation_std
    }
}

impl Eq for AmrMutation {}

impl PartialOrd for AmrMutation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AmrMutation {
    fn cmp(&self, other: &Self) -> Ordering {
        self.gene
            .cmp(&other.gene)
            .then(self.pos_real.cmp(&other.pos_real))
            .then(self.gene_mutation_std.cmp(&other.gene_mutation_std))
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
    pub fn qc(&self, hsp: &Hsp, mutations: &[AmrMutation]) {
        if self.empty() {
            assert!(!self.mutations.is_empty());
            assert!(self.disr.is_none());
            return;
        }

        if self.disr.is_none() {
            assert!(self.len > 0);
            assert!(self.start_ref < self.stop_ref);
            assert!(self.start_target <= hsp.s_int.stop);
            assert!(!self.reference.contains('-'));
            assert!(!self.allele.contains('-'));
            assert!(self
                .reference
                .chars()
                .all(|c| !c.is_ascii_alphabetic() || c.is_ascii_uppercase()));
            assert!(self
                .allele
                .chars()
                .all(|c| !c.is_ascii_alphabetic() || c.is_ascii_uppercase()));
            if self.reference.is_empty() {
                assert!(self.is_frameshift());
                assert!(self.allele.is_empty());
                assert_eq!(self.len, 1);
                assert_eq!(self.stop_ref, self.start_ref + 1);
                assert_eq!(self.neighborhood_mismatch, 0.0);
                assert!(self.mutations.is_empty());
            } else {
                assert!(self.stop_ref <= hsp.q_int.stop);
                assert!(self.reference.len() * hsp.a2q <= self.stop_ref - self.start_ref);
                assert_ne!(self.reference, self.allele);
            }
            for &mutation_idx in &self.mutations {
                let mutation = &mutations[mutation_idx];
                assert!(
                    mutation.pos_real >= hsp.q_int.start && mutation.pos_real <= hsp.q_int.stop
                );
                if mutation.frameshift != NO_INDEX {
                    assert!(hsp.q_prot);
                    assert!(!hsp.s_prot);
                }
            }
        } else {
            assert!(self.mutations.is_empty());
        }
    }

    pub fn save_text(&self, os: &mut dyn Write, mutations: &[AmrMutation]) -> anyhow::Result<()> {
        write!(
            os,
            "{} {} {:?} -> {:?} {}..{} {} {}",
            self.start + 1,
            self.len,
            self.reference,
            self.allele,
            self.start_ref + 1,
            self.stop_ref,
            self.start_target + 1,
            self.neighborhood_mismatch
        )?;
        if let Some(disr) = &self.disr {
            if !disr.empty() {
                write!(os, " {}", disr.genesymbol_raw())?;
            }
        }
        for &mutation_idx in &self.mutations {
            write!(os, " ")?;
            mutations[mutation_idx].save_text(os)?;
        }
        writeln!(os)?;
        Ok(())
    }

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

    pub fn get_mutation_str(&self) -> String {
        if let Some(disr) = &self.disr {
            return disr.genesymbol_raw();
        }
        let allele = if self.allele.is_empty() {
            "DEL".to_string()
        } else if self.allele == "*" {
            TERMINATOR_WORD.to_string()
        } else {
            self.allele.clone()
        };
        format!("{}{}{}", self.reference, self.start_ref + 1, allele)
    }

    pub fn get_stop(&self) -> usize {
        self.start + self.len
    }

    pub fn finish(&mut self, hsp: &Hsp, ref_seq: &str, flanking_len: usize) -> bool {
        self.set_seq(hsp);
        if ref_seq.as_bytes().get(self.start) == Some(&b'-') {
            self.start -= 1;
            self.len += 1;
            self.set_seq(hsp);
        }
        assert!(!self.reference.is_empty());
        self.finish_pos(hsp, flanking_len)
    }

    pub fn finish_pos(&mut self, hsp: &Hsp, flanking_len: usize) -> bool {
        self.set_start_stop_ref(hsp);
        self.set_start_target(hsp);
        self.set_neighborhood_mismatch(hsp, flanking_len);
        self.neighborhood_mismatch <= 0.04
    }

    pub fn set_seq(&mut self, hsp: &Hsp) {
        let stop = self.get_stop();
        assert!(stop <= hsp.qseq.len());
        assert!(stop <= hsp.sseq.len());
        self.reference = hsp.qseq[self.start..stop].replace('-', "");
        self.allele = hsp.sseq[self.start..stop].replace('-', "");
        self.reference.make_ascii_uppercase();
        self.allele.make_ascii_uppercase();
        assert_ne!(self.reference, self.allele);
    }

    pub fn set_start_stop_ref(&mut self, hsp: &Hsp) {
        self.start_ref = hsp.q_int.start;
        for c in hsp.qseq[..self.start].bytes() {
            if c != b'-' {
                self.start_ref += hsp.a2q;
            }
        }
        self.stop_ref = self.start_ref;
        for c in hsp.qseq[self.start..self.get_stop()].bytes() {
            if c != b'-' {
                self.stop_ref += hsp.a2q;
            }
        }
    }

    pub fn set_start_target(&mut self, hsp: &Hsp) {
        self.start_target = hsp.s_int.start;
        if hsp.s_int.strand == 1 {
            for c in hsp.sseq[..self.start].bytes() {
                if c != b'-' {
                    self.start_target += hsp.a2s;
                }
            }
        } else {
            for c in hsp.sseq[self.get_stop()..].bytes().rev() {
                if c != b'-' {
                    self.start_target += hsp.a2s;
                }
            }
        }
    }

    pub fn set_neighborhood_mismatch(&mut self, hsp: &Hsp, flanking_len: usize) {
        self.neighborhood_mismatch = 0.0;
        if flanking_len == 0 {
            return;
        }

        let qseq = hsp.qseq.as_bytes();
        let sseq = hsp.sseq.as_bytes();
        assert_eq!(qseq.len(), sseq.len());

        let mut span = 0usize;
        let mut mismatches = 0usize;

        let mut j = 0usize;
        if self.start != 0 {
            j = self.start - 1;
            while self.start - j <= flanking_len {
                span += 1;
                if sseq[j] != qseq[j] {
                    mismatches += 1;
                }
                if j == 0 {
                    break;
                }
                j -= 1;
            }
        }
        assert!(self.start >= j);
        let missing_left = hsp
            .s_tail(true)
            .min(hsp.q_int.start)
            .min(flanking_len + 1 - (self.start - j));
        span += missing_left;
        mismatches += missing_left;

        let stop = self.get_stop();
        let mut j = stop + 1;
        while j < sseq.len() && j - stop <= flanking_len {
            span += 1;
            if sseq[j] != qseq[j] {
                mismatches += 1;
            }
            j += 1;
        }
        let missing_right = hsp
            .s_tail(false)
            .min(hsp.qlen.saturating_sub(hsp.q_int.stop))
            .min(flanking_len + 1 - (j - stop));
        span += missing_right;
        mismatches += missing_right;

        if span != 0 {
            self.neighborhood_mismatch = mismatches as f64 / span as f64;
        }
    }

    pub fn matches_mutation(&self, mutation: &AmrMutation) -> anyhow::Result<bool> {
        if self.empty() {
            return Ok(false);
        }
        if mutation.pos_real < self.start_ref || mutation.get_stop() > self.stop_ref {
            return Ok(false);
        }

        let head = mutation.pos_real - self.start_ref;
        let tail = self.stop_ref - mutation.get_stop();
        let reference = &self.reference[head..self.reference.len() - tail];
        if reference != mutation.reference {
            anyhow::bail!(
                "Reference sequence has '{}', but mutation is: {}",
                reference,
                mutation.gene_mutation_std
            );
        }
        if head + mutation.allele.len() + tail != self.allele.len() {
            return Ok(false);
        }
        let allele = &self.allele[head..self.allele.len() - tail];
        Ok(allele == mutation.allele)
    }

    pub fn better(&self, other: &SeqChange, rel_identity: f64, other_rel_identity: f64) -> bool {
        if other.has_mutation() != self.has_mutation() {
            return self.has_mutation();
        }
        if self.neighborhood_mismatch != other.neighborhood_mismatch {
            return self.neighborhood_mismatch < other.neighborhood_mismatch;
        }
        rel_identity > other_rel_identity
    }
}

impl PartialEq for SeqChange {
    fn eq(&self, other: &Self) -> bool {
        self.start_target == other.start_target
    }
}

impl Eq for SeqChange {}

impl PartialOrd for SeqChange {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SeqChange {
    fn cmp(&self, other: &Self) -> Ordering {
        self.start_target.cmp(&other.start_target)
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
    pub fn from_blast_line(line: &str, q_prot: bool, s_prot: bool) -> anyhow::Result<Self> {
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

    pub fn qc(&self, mutations: &[AmrMutation]) {
        self.hsp.qc();
        self.ref_mutation.qc();
        for seq_change in &self.seq_changes {
            seq_change.qc(&self.hsp, mutations);
        }
    }

    pub fn save_text(&self, os: &mut dyn Write, mutations: &[AmrMutation]) -> anyhow::Result<()> {
        self.hsp.save_text(os)?;
        if !self.ref_mutation.empty() {
            write!(os, " ")?;
            self.ref_mutation.save_text(os)?;
        }
        write!(os, " #seqChanges:{}", self.seq_changes.len())?;
        if !self.seq_changes.is_empty() {
            writeln!(os)?;
            for seq_change in &self.seq_changes {
                seq_change.save_text(os, mutations)?;
            }
        }
        Ok(())
    }

    pub fn set_seq_changes(
        &mut self,
        ref_mutations: &[AmrMutation],
        flanking_len: usize,
    ) -> anyhow::Result<()> {
        assert!(self.seq_changes.is_empty());
        assert!(!ref_mutations.is_empty());

        if !self.ref_mutation.empty() {
            if self.ref_mutation.frameshift != NO_INDEX && !self.hsp.blastx() {
                return Ok(());
            }
            let index = ref_mutations
                .iter()
                .position(|mutation| mutation == &self.ref_mutation)
                .expect("declarative mutation must exist in reference mutations");
            self.ref_mutation = ref_mutations[index].clone();

            let start = self.ref_mutation2ref_seq_pos();
            if start != NO_INDEX {
                let mut seq_change = SeqChange {
                    start,
                    len: self.ref_mutation.allele.len(),
                    reference: self.ref_mutation.reference.clone(),
                    allele: self.ref_mutation.allele.clone(),
                    ..Default::default()
                };
                if seq_change.finish_pos(&self.hsp, flanking_len) {
                    seq_change.mutations.push(index);
                    self.seq_changes.push(seq_change);
                }
            }
            return Ok(());
        }

        let qseq = self.hsp.qseq.as_bytes();
        let sseq = self.hsp.sseq.as_bytes();
        assert_eq!(qseq.len(), sseq.len());

        let mut seq_change = SeqChange::default();
        let mut in_seq_change = false;
        for i in 0..=qseq.len() {
            if in_seq_change {
                if i == qseq.len() || sseq[i] == qseq[i] {
                    if seq_change.finish(&self.hsp, &self.hsp.qseq, flanking_len) {
                        self.seq_changes.push(seq_change);
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

        let mut j = 0usize;
        while j < ref_mutations.len() && ref_mutations[j].pos_real < self.hsp.q_int.start {
            j += 1;
        }

        let mut start_ref_prev = NO_INDEX;
        for seq_change in &mut self.seq_changes {
            seq_change.qc(&self.hsp, ref_mutations);
            if start_ref_prev != NO_INDEX {
                assert!(start_ref_prev <= seq_change.start_ref);
            }
            while j < ref_mutations.len() {
                let mutation = &ref_mutations[j];
                if mutation.pos_real >= seq_change.stop_ref {
                    break;
                }
                if seq_change.matches_mutation(mutation)? {
                    seq_change.mutations.push(j);
                }
                j += 1;
            }
            start_ref_prev = seq_change.start_ref;
        }

        let mut ref_pos = self.hsp.q_int.start;
        let mut i = 0usize;
        for (mutation_idx, mutation) in ref_mutations.iter().enumerate() {
            while ref_pos < mutation.pos_real {
                assert!(i < qseq.len());
                assert_ne!(qseq[i], b'-');
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
                && self
                    .hsp
                    .qseq
                    .get(i..)
                    .is_some_and(|seq| seq.starts_with(&mutation.reference))
                && self
                    .hsp
                    .sseq
                    .get(i..)
                    .is_some_and(|seq| seq.starts_with(&mutation.reference))
            {
                self.seq_changes.push(SeqChange {
                    mutations: vec![mutation_idx],
                    ..Default::default()
                });
            }
        }

        self.seq_changes.sort();
        Ok(())
    }

    pub fn ref_mutation2ref_seq_pos(&self) -> usize {
        assert!(!self.ref_mutation.empty());

        if self.ref_mutation.pos_real < self.hsp.q_int.start
            || self.ref_mutation.pos_real + self.ref_mutation.allele.len() > self.hsp.q_int.stop
        {
            return NO_INDEX;
        }

        let qseq = self.hsp.qseq.as_bytes();
        let sseq = self.hsp.sseq.as_bytes();
        let allele = self.ref_mutation.allele.as_bytes();
        let mut pos = self.hsp.q_int.start;
        let mut frameshift_i = NO_INDEX;
        let mut i = 0usize;

        while i < qseq.len() {
            assert!(pos <= self.ref_mutation.pos_real);
            if qseq[i] != b'-' {
                if pos == self.ref_mutation.pos_real {
                    let align_stop = i + allele.len();
                    if align_stop <= qseq.len()
                        && &qseq[i..align_stop] == allele
                        && &sseq[i..align_stop] == allele
                        && (pos == 0 || (i > 0 && qseq[i - 1] == sseq[i - 1]))
                    {
                        if self.ref_mutation.frameshift != NO_INDEX {
                            frameshift_i = i;
                            i += allele.len();
                            pos = self.ref_mutation.get_stop();
                            break;
                        }
                        let ref_stop = pos + allele.len();
                        assert!(align_stop <= qseq.len());
                        if ref_stop == self.hsp.qlen
                            || (align_stop < qseq.len() && qseq[align_stop] == sseq[align_stop])
                        {
                            return i;
                        }
                    }
                    return NO_INDEX;
                }
                pos += 1;
            }
            i += 1;
        }

        while i < qseq.len() {
            assert_ne!(frameshift_i, NO_INDEX);
            if qseq[i] != b'-' {
                assert!(self.ref_mutation.get_stop() > 0);
                if pos == self.ref_mutation.get_stop() - 1 + self.ref_mutation.frameshift {
                    if qseq[i] == b'*' && sseq[i] == b'*' {
                        return frameshift_i;
                    }
                    return NO_INDEX;
                }
                pos += 1;
            }
            i += 1;
        }

        NO_INDEX
    }
}
