// Sequence analysis types — port of seq.hpp
// Focus on types used by amr_report, amrfinder, dna_mutation

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Frame: -3, -2, -1, 1, 2, 3 (0 = unknown)
pub type Frame = i8;

/// Strand: 1 (top/+), -1 (bottom/-), 0 (unknown)
pub type Strand = i8;

/// Returns true if `frame` is a valid reading frame (-3..=-1 or 1..=3).
pub fn is_frame(frame: Frame) -> bool {
    (-3..=-1).contains(&frame) || (1..=3).contains(&frame)
}

/// Returns true if `strand` is a valid strand value (-1 or 1).
pub fn is_strand(strand: Strand) -> bool {
    strand == -1 || strand == 1
}

/// Converts a strand value to '+', '-', or '?' for display.
pub fn strand2char(strand: Strand) -> char {
    match strand {
        -1 => '-',
        1 => '+',
        _ => '?',
    }
}

/// No-index sentinel (matches C++ no_index)
pub const NO_INDEX: usize = usize::MAX;

// --- Nucleotide/amino acid matching ---

pub const DNA_ALPHABET: &str = "acgt";
pub const EXT_SPARSE_DNA_ALPHABET: &str = "-acgtbdhkmnrsvwyACGTBDHKMNRSVWY";
pub const EXT_DNA_ALPHABET: &str = "acgtbdhkmnrsvwyACGTBDHKMNRSVWY";
pub const DNA_WILDCARDS: &str = "bdhkmnrsvwyACGTBDHKMNRSVWY";

pub const PEPTIDE_ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWY";
pub const EXT_PEPTIDE_ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWYXBZJUO";
pub const EXT_SPARSE_TERM_PEPTIDE_ALPHABET: &str = "-ACDEFGHIKLMNPQRSTVWYXBZJUO*";
pub const TERMINATOR: &str = "*";
pub const EXT_TERM_PEPTIDE_ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWYXBZJUO*";
pub const PEPTIDE_WILDCARDS: &str = "XBZJUO";

pub const TERMINATOR_WORD: &str = "Ter";

pub type Gencode = u8;

pub fn count_insertions(seq: &str) -> usize {
    seq.chars().filter(|&c| c == '-').count()
}

pub fn sparse_seq_len(seq: &str) -> usize {
    seq.len() - count_insertions(seq)
}

pub fn codon2aa(
    codon: &str,
    gencode: Gencode,
    lowercase_possible_start_codon: bool,
) -> anyhow::Result<char> {
    assert!(codon.len() >= 3);
    let c = codon.to_ascii_lowercase();
    let b = c.as_bytes();
    let c0 = b[0] as char;
    let c1 = b[1] as char;
    let c2 = b[2] as char;

    let mut aa = match c0 {
        't' => match c1 {
            't' => match c2 {
                't' | 'c' | 'y' => 'F',
                'a' | 'g' | 'r' => 'L',
                _ => 'X',
            },
            'c' => 'S',
            'a' => match c2 {
                't' | 'c' | 'y' => 'Y',
                'a' | 'g' | 'r' => '*',
                _ => 'X',
            },
            'g' => match c2 {
                't' | 'c' | 'y' => 'C',
                'a' if gencode == 4 || gencode == 25 => 'W',
                'a' => '*',
                'g' => 'W',
                _ => 'X',
            },
            'r' if c2 == 'a' => '*',
            _ => 'X',
        },
        'c' => match c1 {
            't' => 'L',
            'c' => 'P',
            'a' => match c2 {
                't' | 'c' | 'y' => 'H',
                'a' | 'g' | 'r' => 'Q',
                _ => 'X',
            },
            'g' => 'R',
            _ => 'X',
        },
        'a' => match c1 {
            't' => match c2 {
                't' | 'c' | 'a' | 'y' | 'w' | 'm' | 'h' => 'I',
                'g' => 'M',
                _ => 'X',
            },
            'c' => 'T',
            'a' => match c2 {
                't' | 'c' | 'y' => 'N',
                'a' | 'g' | 'r' => 'K',
                _ => 'X',
            },
            'g' => match c2 {
                't' | 'c' | 'y' => 'S',
                'a' | 'g' | 'r' => 'R',
                _ => 'X',
            },
            _ => 'X',
        },
        'g' => match c1 {
            't' => 'V',
            'c' => 'A',
            'a' => match c2 {
                't' | 'c' | 'y' => 'D',
                'a' | 'g' | 'r' => 'E',
                _ => 'X',
            },
            'g' => 'G',
            _ => 'X',
        },
        'm' if c1 == 'g' && matches!(c2, 'a' | 'g' | 'r') => 'R',
        'r' if c1 == 'a' && matches!(c2, 't' | 'c' | 'y') => 'B',
        's' if c1 == 'a' && matches!(c2, 'a' | 'g' | 'r') => 'Z',
        'y' if c1 == 't' && matches!(c2, 'a' | 'g' | 'r') => 'L',
        _ => 'X',
    };
    assert!(EXT_TERM_PEPTIDE_ALPHABET.contains(aa));

    if lowercase_possible_start_codon {
        let startcodon = match gencode {
            1 => c0 == 'a' && c1 == 't' && c2 == 'g',
            4 => {
                (c1 == 't' && c2 == 'g')
                    || (c0 == 'a' && c1 == 't')
                    || (c0 == 't' && c1 == 't' && c2 == 'a')
            }
            11 => (c1 == 't' && c2 == 'g') || (c0 == 'a' && c1 == 't'),
            25 => c0 != 'c' && c1 == 't' && c2 == 'g',
            _ => anyhow::bail!("codon2aa: Genetic code {} is not implemented", gencode),
        };
        if startcodon {
            aa = aa.to_ascii_lowercase();
        }
    }

    Ok(aa)
}

pub fn alphabet2_pos(alphabet: &str, c: char) -> usize {
    alphabet
        .find(c)
        .unwrap_or_else(|| panic!("Bad character {c} for alphabet: {alphabet}"))
}

pub fn is_ambig_nucl(c: char) -> bool {
    DNA_WILDCARDS.contains(c)
}

pub fn is_ambig_aa(c: char) -> bool {
    PEPTIDE_WILDCARDS.contains(c)
}

pub fn is_ambig(c: char, aa: bool) -> bool {
    if aa {
        is_ambig_aa(c)
    } else {
        is_ambig_nucl(c)
    }
}

pub fn nuc2num(wild_nucleotide: char) -> usize {
    match wild_nucleotide.to_ascii_lowercase() {
        'a' => 0,
        'c' => 1,
        'g' => 2,
        't' => 3,
        c if is_ambig_nucl(c) => 4,
        _ => panic!("{wild_nucleotide:?} is not a nucleotide"),
    }
}

pub fn wild2nucleotides(wild_nucleotide: char, acgtb: &mut [bool; 5]) -> u8 {
    *acgtb = [false; 5];

    if wild_nucleotide == ' ' {
        return 0;
    }
    if wild_nucleotide == '-' {
        acgtb[4] = true;
        return 1;
    }

    acgtb[4] = wild_nucleotide.is_ascii_uppercase();
    match wild_nucleotide.to_ascii_lowercase() {
        'a' => acgtb[0] = true,
        'c' => acgtb[1] = true,
        'g' => acgtb[2] = true,
        't' => acgtb[3] = true,
        'm' => acgtb[0..2].fill(true),
        'r' => {
            acgtb[0] = true;
            acgtb[2] = true;
        }
        'w' => {
            acgtb[0] = true;
            acgtb[3] = true;
        }
        's' => acgtb[1..3].fill(true),
        'y' => {
            acgtb[1] = true;
            acgtb[3] = true;
        }
        'k' => acgtb[2..4].fill(true),
        'v' => acgtb[0..3].fill(true),
        'h' => {
            acgtb[0] = true;
            acgtb[1] = true;
            acgtb[3] = true;
        }
        'd' => {
            acgtb[0] = true;
            acgtb[2] = true;
            acgtb[3] = true;
        }
        'b' => acgtb[1..4].fill(true),
        'n' => acgtb[0..4].fill(true),
        _ => panic!("Bad wildNucleotide: {}", wild_nucleotide as u32),
    }

    acgtb.iter().filter(|present| **present).count() as u8
}

pub fn nucleotides2wild(acgtb: &[bool; 5]) -> char {
    let [a, c, g, t, blank] = *acgtb;
    let res = match (a, c, g, t) {
        (true, true, true, true) => 'n',
        (true, true, true, false) => 'v',
        (true, true, false, true) => 'h',
        (true, true, false, false) => 'm',
        (true, false, true, true) => 'd',
        (true, false, true, false) => 'r',
        (true, false, false, true) => 'w',
        (true, false, false, false) => 'a',
        (false, true, true, true) => 'b',
        (false, true, true, false) => 's',
        (false, true, false, true) => 'y',
        (false, true, false, false) => 'c',
        (false, false, true, true) => 'k',
        (false, false, true, false) => 'g',
        (false, false, false, true) => 't',
        (false, false, false, false) => '-',
    };

    if blank {
        res.to_ascii_uppercase()
    } else if res == '-' {
        ' '
    } else {
        res
    }
}

pub fn complementary_nucleotide(wild_nucleotide: char) -> char {
    let mut r = match wild_nucleotide.to_ascii_lowercase() {
        'a' => 't',
        'c' => 'g',
        'g' => 'c',
        't' => 'a',
        'm' => 'k',
        'r' => 'y',
        'w' => 'w',
        's' => 's',
        'y' => 'r',
        'k' => 'm',
        'v' => 'b',
        'h' => 'd',
        'd' => 'h',
        'b' => 'v',
        'n' => 'n',
        _ => panic!("Bad wild nucleotide {}", wild_nucleotide as u32),
    };
    if wild_nucleotide.is_ascii_uppercase() {
        r = r.to_ascii_uppercase();
    }
    r
}

pub fn get_union_nucleotide(char_set: &str) -> char {
    if char_set.is_empty() {
        return ' ';
    }

    let mut acgtb = [false; 5];
    for c in char_set.chars() {
        let mut c_acgtb = [false; 5];
        wild2nucleotides(c, &mut c_acgtb);
        for i in 0..5 {
            acgtb[i] |= c_acgtb[i];
        }
    }
    nucleotides2wild(&acgtb)
}

pub fn get_intersect_nucleotide(char_set: &str) -> char {
    if char_set.is_empty() {
        return ' ';
    }

    let mut acgtb = [true; 5];
    for c in char_set.chars() {
        let mut c_acgtb = [false; 5];
        wild2nucleotides(c, &mut c_acgtb);
        for i in 0..5 {
            acgtb[i] &= c_acgtb[i];
        }
    }
    nucleotides2wild(&acgtb)
}

/// Check if two nucleotides match (considering IUPAC ambiguity codes)
pub fn nucleotide_match(wild_nucleotide1: char, wild_nucleotide2: char) -> bool {
    assert_ne!(wild_nucleotide1, '\0');
    assert_ne!(wild_nucleotide2, '\0');

    wild_nucleotide1 == wild_nucleotide2
        || get_intersect_nucleotide(&format!("{wild_nucleotide1}{wild_nucleotide2}")) != ' '
}

pub fn dna2codons_len(dna_len: usize) -> usize {
    if dna_len == 0 {
        0
    } else {
        (dna_len - 1) / 3 + 1
    }
}

pub fn aa2num(wild_aminoacid: char) -> usize {
    let c = wild_aminoacid.to_ascii_uppercase();
    if let Some(i) = PEPTIDE_ALPHABET.find(c) {
        return i;
    }
    if wild_aminoacid == TERMINATOR.chars().next().unwrap() {
        return 20;
    }
    if is_ambig_aa(c) {
        return 21;
    }
    panic!("Bad wildAminoacid: {wild_aminoacid}");
}

pub fn more_general_aminoacid(wild_aminoacid1: char, wild_aminoacid2: char) -> bool {
    if wild_aminoacid1 == wild_aminoacid2 {
        return true;
    }

    match wild_aminoacid1 {
        'x' => true,
        'b' => matches!(wild_aminoacid2, 'd' | 'n'),
        'z' => matches!(wild_aminoacid2, 'e' | 'q'),
        'j' => matches!(wild_aminoacid2, 'i' | 'l'),
        'u' | 'o' => wild_aminoacid2 == '*',
        '*' => matches!(wild_aminoacid2, 'u' | 'o'),
        _ => false,
    }
}

/// Check if two amino acids match
pub fn aa_match(a: char, b: char) -> bool {
    more_general_aminoacid(a, b) || more_general_aminoacid(b, a)
}

pub const FASTA_LINE_LEN: usize = 80;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Seq {
    pub name: String,
    pub seq: String,
    pub sparse: bool,
}

impl Seq {
    pub fn new(name: &str, seq: &str, sparse: bool) -> Self {
        let mut this = Seq {
            name: name.to_string(),
            seq: seq.to_string(),
            sparse,
        };
        this.qc_name().unwrap();
        if !this.sparse {
            this.un_sparse();
        }
        this
    }

    pub fn with_len(name: &str, seq_len: usize, sparse: bool) -> Self {
        let this = Seq {
            name: name.to_string(),
            seq: "\0".repeat(seq_len),
            sparse,
        };
        this.qc_name().unwrap();
        this
    }

    pub fn from_multifasta(
        fasta: &mut Multifasta,
        reserve_len: usize,
        sparse: bool,
        make_upper: bool,
    ) -> anyhow::Result<Self> {
        assert!(!fasta.eof);
        let header = fasta.line.as_ref().expect("Multifasta current line exists");
        assert!(!header.is_empty());
        assert!(header.starts_with('>'));

        let name = header[1..].replace('\t', " ");
        let mut this = Seq {
            name,
            seq: String::with_capacity(reserve_len),
            sparse,
        };
        this.qc_name()?;

        while fasta.next_line()? {
            let line = fasta.line.as_ref().expect("Multifasta current line exists");
            if line.trim().is_empty() || line.starts_with('>') {
                break;
            }
            this.seq
                .push_str(line.trim_end_matches(|c: char| c.is_ascii_whitespace()));
        }

        while !fasta.eof
            && fasta
                .line
                .as_ref()
                .is_some_and(|line| line.trim().is_empty())
        {
            fasta.next_line()?;
        }

        if make_upper {
            this.seq = this.seq.to_ascii_uppercase();
        } else {
            this.seq = this.seq.to_ascii_lowercase();
        }
        if !this.sparse {
            this.un_sparse();
        }
        Ok(this)
    }

    pub fn qc_name(&self) -> anyhow::Result<()> {
        anyhow::ensure!(!self.name.trim().is_empty(), "Blank sequence name");
        anyhow::ensure!(!self.name.contains('\t'), "Sequence name contains tab");
        Ok(())
    }

    pub fn qc_alphabet(&self, alphabet: &str) -> anyhow::Result<()> {
        if let Some((i, c)) = self
            .seq
            .chars()
            .enumerate()
            .find(|(_, c)| !alphabet.contains(*c))
        {
            anyhow::bail!(
                "Bad sequence character in {:?}: ASCII={} pos={}",
                self.get_id(),
                c as u32,
                i + 1
            );
        }
        Ok(())
    }

    pub fn get_id_size(&self) -> usize {
        self.name.find(' ').unwrap_or(self.name.len())
    }

    pub fn get_taxon_start_in(s: &str) -> Option<usize> {
        let bytes = s.as_bytes();
        let mut brackets = 0usize;
        for i in (0..bytes.len()).rev() {
            match bytes[i] {
                b']' => brackets += 1,
                b'[' => {
                    if brackets == 1 {
                        return Some(i);
                    }
                    if brackets == 0 {
                        return None;
                    }
                    brackets -= 1;
                }
                _ => {}
            }
            if brackets == 0 {
                return None;
            }
        }
        None
    }

    pub fn get_taxon_start(&self) -> Option<usize> {
        Self::get_taxon_start_in(&self.name)
    }

    pub fn get_id(&self) -> String {
        self.name[..self.get_id_size()].to_string()
    }

    pub fn get_gi(&self) -> anyhow::Result<i64> {
        let prefix = "gi|";
        let (start, end) = if let Some(mut start) = self.name.find(prefix) {
            start += prefix.len();
            let end = self.name[start..]
                .find('|')
                .map(|pos| start + pos)
                .expect("gi| id has closing pipe");
            (start, end)
        } else {
            let end = self.name.find(' ').unwrap_or(self.name.len());
            (0, end)
        };
        Ok(self.name[start..end].parse()?)
    }

    pub fn append_id(&mut self, suffix: &str) {
        let pos = self.get_id_size();
        if !suffix.is_empty() {
            self.name.insert_str(pos, &format!("-{suffix}"));
        }
    }

    pub fn get_description(&self, trim_taxon: bool) -> String {
        let id_size = self.get_id_size();
        if id_size == self.name.len() {
            return String::new();
        }
        assert_eq!(self.name.as_bytes()[id_size], b' ');
        let mut desc = self.name[id_size + 1..].to_string();
        if trim_taxon {
            if let Some(taxon_start) = self.get_taxon_start() {
                desc.truncate(taxon_start.saturating_sub(id_size + 1));
            }
        }
        desc.trim().to_string()
    }

    pub fn get_xs<F>(&self, is_ambiguous: F) -> usize
    where
        F: Fn(char) -> bool,
    {
        self.seq.chars().filter(|&c| is_ambiguous(c)).count()
    }

    pub fn get_contiguous_xs<F>(&self, is_ambiguous: F) -> usize
    where
        F: Fn(char) -> bool,
    {
        let mut max_len = 0usize;
        let mut len = 0usize;
        for c in self.seq.chars() {
            if is_ambiguous(c) {
                len += 1;
            } else {
                max_len = max_len.max(len);
                len = 0;
            }
        }
        max_len.max(len)
    }

    pub fn un_sparse(&mut self) {
        self.seq.retain(|c| c != '-');
        self.sparse = false;
    }

    pub fn get_char_count(&self) -> BTreeMap<char, usize> {
        let mut counts = BTreeMap::new();
        for c in self.seq.chars() {
            *counts.entry(c).or_insert(0) += 1;
        }
        counts
    }
}

#[derive(Debug, Clone)]
pub struct Multifasta {
    lines: Vec<String>,
    index: usize,
    pub line: Option<String>,
    pub line_num: usize,
    pub eof: bool,
    pub aa: bool,
    pub display_period: usize,
}

impl Multifasta {
    pub fn new<P: AsRef<Path>>(fname: P, aa: bool, display_period: usize) -> anyhow::Result<Self> {
        let file = File::open(fname)?;
        let reader = BufReader::new(file);
        let mut lines = Vec::new();
        for line in reader.lines() {
            lines.push(line?);
        }
        let mut this = Multifasta {
            lines,
            index: 0,
            line: None,
            line_num: 0,
            eof: false,
            aa,
            display_period,
        };
        this.next_line()?;
        this.qc_new_seq()?;
        Ok(this)
    }

    pub fn qc_new_seq(&self) -> anyhow::Result<()> {
        if !self.eof
            && !self
                .line
                .as_ref()
                .is_some_and(|line| !line.is_empty() && line.starts_with('>'))
        {
            anyhow::bail!("Error in Multifasta, line {}", self.line_num);
        }
        Ok(())
    }

    pub fn next(&self) -> anyhow::Result<bool> {
        self.qc_new_seq()?;
        Ok(!self.eof)
    }

    fn next_line(&mut self) -> anyhow::Result<bool> {
        if self.index >= self.lines.len() {
            self.line = None;
            self.eof = true;
            return Ok(false);
        }
        self.line = Some(self.lines[self.index].clone());
        self.index += 1;
        self.line_num = self.index;
        self.eof = false;
        Ok(true)
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Dna {
    pub name: String,
    pub seq: String,
    pub sparse: bool,
}

impl Dna {
    pub const STD_MIN_COMPLEXITY: f64 = 2.0;

    pub fn new(name: &str, seq: &str, sparse: bool) -> Self {
        let base = Seq::new(name, seq, sparse);
        Dna {
            name: base.name,
            seq: base.seq,
            sparse: base.sparse,
        }
    }

    pub fn with_len(name: &str, seq_len: usize, sparse: bool) -> Self {
        let base = Seq::with_len(name, seq_len, sparse);
        Dna {
            name: base.name,
            seq: base.seq,
            sparse: base.sparse,
        }
    }

    pub fn from_multifasta(
        fasta: &mut Multifasta,
        reserve_len: usize,
        sparse: bool,
    ) -> anyhow::Result<Self> {
        assert!(!fasta.aa);
        let base = Seq::from_multifasta(fasta, reserve_len, sparse, false)?;
        let this = Dna {
            name: base.name,
            seq: base.seq,
            sparse: base.sparse,
        };
        this.qc()?;
        Ok(this)
    }

    pub fn get_seq_alphabet(&self) -> &'static str {
        if self.sparse {
            EXT_SPARSE_DNA_ALPHABET
        } else {
            EXT_DNA_ALPHABET
        }
    }

    pub fn qc(&self) -> anyhow::Result<()> {
        let base = Seq {
            name: self.name.clone(),
            seq: self.seq.clone(),
            sparse: self.sparse,
        };
        base.qc_name()?;
        base.qc_alphabet(self.get_seq_alphabet())
    }

    pub fn is_ambiguous(c: char) -> bool {
        !DNA_ALPHABET.contains(c)
    }

    pub fn get_id_size(&self) -> usize {
        self.name.find(' ').unwrap_or(self.name.len())
    }

    pub fn get_id(&self) -> String {
        self.name[..self.get_id_size()].to_string()
    }

    pub fn get_taxon_start(&self) -> Option<usize> {
        Seq::get_taxon_start_in(&self.name)
    }

    pub fn get_description(&self, trim_taxon: bool) -> String {
        Seq {
            name: self.name.clone(),
            seq: self.seq.clone(),
            sparse: self.sparse,
        }
        .get_description(trim_taxon)
    }

    pub fn get_xs(&self) -> usize {
        self.seq.chars().filter(|&c| Self::is_ambiguous(c)).count()
    }

    pub fn get_contiguous_xs(&self) -> usize {
        let mut max_len = 0usize;
        let mut len = 0usize;
        for c in self.seq.chars() {
            if Self::is_ambiguous(c) {
                len += 1;
            } else {
                max_len = max_len.max(len);
                len = 0;
            }
        }
        max_len.max(len)
    }

    pub fn get_char_count(&self) -> BTreeMap<char, usize> {
        Seq {
            name: self.name.clone(),
            seq: self.seq.clone(),
            sparse: self.sparse,
        }
        .get_char_count()
    }

    pub fn un_sparse(&mut self) {
        self.seq.retain(|c| c != '-');
        self.sparse = false;
    }

    pub fn get_complexity_int(&self, start: usize, mut end: usize) -> f64 {
        if self.seq.is_empty() {
            return 0.0;
        }

        let mut dinuc_freq = [[0usize; 4]; 4];
        let mut n = 0usize;
        end = end.min(self.seq.len());
        let bytes = self.seq.as_bytes();
        for i in start + 1..end {
            let a = nuc2num(bytes[i - 1] as char);
            let b = nuc2num(bytes[i] as char);
            if a < 4 && b < 4 {
                dinuc_freq[a][b] += 1;
                n += 1;
            }
        }

        let mut s = 0.0;
        for row in dinuc_freq {
            for count in row {
                if count != 0 {
                    let p = count as f64 / n as f64;
                    s -= p * p.ln();
                }
            }
        }
        assert!(s >= 0.0);
        s
    }

    pub fn get_complexity(&self) -> f64 {
        self.get_complexity_int(0, self.seq.len())
    }

    pub fn mono_nuc2n(&mut self, repeat_min: usize) -> usize {
        let mut n = 0usize;
        let mut chars: Vec<char> = self.seq.chars().collect();
        let mut acgt_size = [0usize; 4];

        for i in 0..chars.len() {
            let mut acgtb = [false; 5];
            wild2nucleotides(chars[i], &mut acgtb);
            for j in 0..4 {
                if acgtb[j] {
                    acgt_size[j] += 1;
                } else {
                    let len = acgt_size[j];
                    if len >= repeat_min {
                        for c in &mut chars[i - len..i] {
                            *c = 'n';
                        }
                        n += len;
                    }
                    acgt_size[j] = 0;
                }
            }
        }

        for len in acgt_size {
            if len >= repeat_min {
                let start = chars.len() - len;
                for c in &mut chars[start..] {
                    *c = 'n';
                }
                n += len;
            }
        }

        self.seq = chars.into_iter().collect();
        n
    }

    pub fn make_complementary(&self) -> Dna {
        let seq = self
            .seq
            .chars()
            .rev()
            .map(complementary_nucleotide)
            .collect();
        Dna {
            name: format!("{}.rev", self.get_id()),
            seq,
            sparse: self.sparse,
        }
    }

    pub fn reverse(&mut self) {
        let dna_rev = self.make_complementary();
        assert_eq!(self.seq.len(), dna_rev.seq.len());
        self.seq = dna_rev.seq;
    }

    pub fn make_peptide(
        &self,
        frame: Frame,
        gencode: Gencode,
        lowercase_possible_start_codon: bool,
        first_start_codon2m: bool,
    ) -> anyhow::Result<(Peptide, usize)> {
        assert!(is_frame(frame));
        assert!(!first_start_codon2m || lowercase_possible_start_codon);

        let peptide_name = format!("{}.fr{}", self.get_id(), frame);
        let frame_offset = frame.unsigned_abs() as usize - 1;
        assert!(frame_offset <= 2);
        let aa_seq_len = if self.seq.len() >= 3 {
            (self.seq.len() - frame_offset) / 3
        } else {
            0
        };

        let comp_dna;
        let dna = if frame < 0 {
            comp_dna = self.make_complementary();
            &comp_dna
        } else {
            self
        };
        assert_eq!(dna.seq.len(), self.seq.len());

        let translation_start = if frame > 0 {
            frame_offset
        } else {
            dna.seq.len() - frame_offset
        };

        let mut peptide = Peptide::with_len(&peptide_name, aa_seq_len, false);
        peptide.pseudo = true;
        for j in 0..aa_seq_len {
            let i = frame_offset + j * 3;
            peptide.seq.replace_range(
                j..j + 1,
                &codon2aa(&dna.seq[i..i + 3], gencode, lowercase_possible_start_codon)?.to_string(),
            );
        }

        if !peptide.seq.is_empty()
            && first_start_codon2m
            && peptide.seq.as_bytes()[0].is_ascii_lowercase()
        {
            peptide.seq.replace_range(0..1, "m");
        }

        Ok((peptide, translation_start))
    }

    pub fn cds2prot(
        &self,
        gencode: Gencode,
        trunc5: bool,
        trunc3: bool,
        has_stop_codon: bool,
        allow_extra_stop_codon: bool,
    ) -> anyhow::Result<Peptide> {
        let (mut pep, _) = self.make_peptide(1, gencode, true, true)?;
        if let Some(name) = pep.name.strip_suffix(".fr1") {
            pep.name = name.to_string();
        }

        if pep.seq.is_empty() {
            anyhow::bail!("cds2prot: empty sequence");
        }
        if !trunc5 && !trunc3 && pep.seq.len() * 3 != self.seq.len() {
            anyhow::bail!("cds2prot: incomplete codon");
        }
        if !trunc5 && !pep.seq.starts_with('m') {
            anyhow::bail!("cds2prot: no start codon");
        }

        let star_pos = pep.seq.find('*');
        if has_stop_codon {
            if let Some(star_pos) = star_pos {
                if star_pos == pep.seq.len() - 1 {
                    pep.trim_stop();
                } else if allow_extra_stop_codon {
                    pep.seq.remove(star_pos);
                } else {
                    anyhow::bail!(
                        "cds2prot: stop codon at peptide position {}\n{}",
                        star_pos + 1,
                        &self.seq[star_pos * 3..]
                    );
                }
            } else if !trunc3 {
                anyhow::bail!("cds2prot: no stop codon");
            }
        } else if star_pos.is_some() {
            anyhow::bail!("cds2prot: has a stop codon");
        }

        Ok(pep)
    }

    pub fn get_peptides(
        &self,
        frame: Frame,
        gencode: Gencode,
        len_min: usize,
    ) -> anyhow::Result<Vec<Peptide>> {
        assert!(is_frame(frame));
        assert!(len_min > 0);

        let mut peps = Vec::with_capacity(self.seq.len() / 300 + 1);
        let mut dna_seq = self.seq.clone();
        if frame < 0 {
            reverse_dna(&mut dna_seq);
        }

        let mut start = frame.unsigned_abs() as usize - 1;
        let mut pep_seq = String::new();

        for i in (start..dna_seq.len().saturating_sub(2)).step_by(3) {
            let aa = codon2aa(&dna_seq[i..i + 3], gencode, true)?;
            let stop = i;
            if aa == TERMINATOR.chars().next().unwrap() {
                if pep_seq.len() >= len_min {
                    let (start_, stop_) = if frame < 0 {
                        (self.seq.len() - stop, self.seq.len() - start)
                    } else {
                        (start, stop)
                    };
                    peps.push(Peptide::new(
                        &format!("{}:{}..{}", self.get_id(), start_ + 1, stop_),
                        &pep_seq,
                        false,
                    ));
                }
                pep_seq.clear();
            } else if !pep_seq.is_empty() {
                pep_seq.push(aa.to_ascii_uppercase());
            } else if Peptide::is_start_aa(aa) {
                pep_seq = "M".to_string();
                start = i;
            }
        }

        Ok(peps)
    }

    pub fn contains_ambiguity(&self) -> bool {
        self.seq.chars().any(|c| DNA_WILDCARDS.contains(c))
    }

    pub fn get_ambiguous_prefix_end(&self) -> usize {
        let mut i = 0usize;
        while i < self.seq.len() && DNA_WILDCARDS.contains(self.seq.as_bytes()[i] as char) {
            i += 1;
        }
        i
    }

    pub fn get_ambiguous_suffix_start(&self) -> usize {
        let mut i = self.seq.len();
        while i > 0 && DNA_WILDCARDS.contains(self.seq.as_bytes()[i - 1] as char) {
            i -= 1;
        }
        i
    }

    pub fn delete_ambiguous_prefix(&mut self) -> bool {
        let ambiguous_prefix_end = self.get_ambiguous_prefix_end();
        self.seq.drain(0..ambiguous_prefix_end);
        ambiguous_prefix_end > 0
    }

    pub fn delete_ambiguous_suffix(&mut self) -> bool {
        let ambiguous_suffix_start = self.get_ambiguous_suffix_start();
        let trimmed = ambiguous_suffix_start < self.seq.len();
        self.seq.truncate(ambiguous_suffix_start);
        trimmed
    }

    pub fn trim_n(&mut self, gap_coded: bool) {
        while self.seq.ends_with(|c| c == 'N' || (!gap_coded && c == 'n')) {
            self.seq.pop();
        }

        let mut k = 0usize;
        while k < self.seq.len() {
            let c = self.seq.as_bytes()[k] as char;
            if c == 'N' || (!gap_coded && c == 'n') {
                k += 1;
            } else {
                break;
            }
        }
        self.seq.drain(0..k);
    }

    pub fn poly_nuc_window(&self, start: usize, window_len: usize, nucleotide: char) -> bool {
        assert!(start <= self.seq.len());
        assert!(DNA_ALPHABET.contains(nucleotide));

        let end = (start + window_len).min(self.seq.len());
        assert!(start <= end);
        if start == end {
            return false;
        }

        let mismatch = self.seq[start..end]
            .chars()
            .filter(|&c| !nucleotide_match(c, nucleotide))
            .count();
        mismatch as f64 / (end - start) as f64 <= 0.2
    }

    pub fn find_poly_nuc_end(&self, window_len: usize, nucleotide: char) -> usize {
        assert!(window_len > 0);
        assert!(DNA_ALPHABET.contains(nucleotide));

        let mut last_init_start = self.seq.len().saturating_sub(window_len);
        let mut not_poly_nuc = 0usize;
        for init_start in (0..=last_init_start).rev() {
            if self.poly_nuc_window(init_start, window_len, nucleotide) {
                last_init_start = init_start;
                not_poly_nuc = 0;
            } else {
                not_poly_nuc += 1;
                if not_poly_nuc as f64 > window_len as f64 * 1.3 {
                    break;
                }
            }
        }

        for start in last_init_start..self.seq.len() {
            if self.poly_nuc_window(start, window_len.min(self.seq.len() - start), nucleotide) {
                let mut nucleotide_num = 0usize;
                let mut i = start;
                while i < self.seq.len() {
                    if nucleotide_match(self.seq.as_bytes()[i] as char, nucleotide) {
                        nucleotide_num += 1;
                        if nucleotide_num >= 2 {
                            return i + 1 - nucleotide_num;
                        }
                    } else {
                        nucleotide_num = 0;
                    }
                    i += 1;
                }
                return i - nucleotide_num;
            }
        }

        self.seq.len()
    }

    pub fn poly_a_window(&self, start: usize, window_len: usize) -> bool {
        self.poly_nuc_window(start, window_len, 'a')
    }

    pub fn find_poly_a(&self, window_len: usize) -> usize {
        self.find_poly_nuc_end(window_len, 'a')
    }

    pub fn remove_poly_a(&mut self, window_len: usize) -> usize {
        let poly_a_start = self.find_poly_a(window_len);
        let poly_a_len = self.seq.len() - poly_a_start;
        self.seq.truncate(poly_a_start);
        poly_a_len
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Cds {
    pub start: usize,
    pub stop: usize,
    pub ref_prot: String,
    pub positives: f64,
    pub prevs: Vec<usize>,
    pub best_prev: Option<usize>,
    pub sum_len: isize,
}

impl Cds {
    pub const PEPTIDE_SIZE_MIN: usize = 20;
    pub const PROMOTER_MIN: usize = 60;
    pub const SIZE_MIN: usize = 3 * (Self::PEPTIDE_SIZE_MIN + 1);

    pub fn new(start_arg: usize, stop_arg: usize) -> Self {
        Cds {
            start: start_arg,
            stop: stop_arg,
            positives: 1.0,
            ..Cds::default()
        }
    }

    pub fn new_fields(
        start_arg: usize,
        stop_arg: usize,
        ref_prot_arg: String,
        positives_arg: f64,
    ) -> Self {
        Cds {
            start: start_arg,
            stop: stop_arg,
            ref_prot: ref_prot_arg,
            positives: positives_arg,
            ..Cds::default()
        }
    }

    pub fn qc(&self) {
        assert!(!self.ref_prot.is_empty());
        assert_ne!(self.start, self.stop);
        assert_eq!(self.size() % 3, 0);
        assert!(self.size_effective() > 0);
        assert!(self.positives > 0.5);
        assert!(self.positives <= 1.0);
    }

    pub fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        write!(
            os,
            "{}\t{}\t{}\t{}",
            self.left(),
            self.right(),
            self.strand() as i32,
            (self.size() - 1) / 3
        )?;
        Ok(())
    }

    pub fn strand(&self) -> bool {
        self.start < self.stop
    }

    pub fn left(&self) -> usize {
        self.start.min(self.stop)
    }

    pub fn right(&self) -> usize {
        self.start.max(self.stop)
    }

    pub fn start_human(&self) -> usize {
        (if self.strand() { self.start } else { self.stop }) + 1
    }

    pub fn stop_human(&self) -> usize {
        if self.strand() {
            self.stop
        } else {
            self.start
        }
    }

    pub fn size(&self) -> usize {
        self.right() - self.left()
    }

    pub fn size_effective(&self) -> isize {
        self.size() as isize - (Self::SIZE_MIN as isize - 3)
    }

    pub fn get_overlap(&self, other: &Cds) -> usize {
        if other.left() >= self.right() {
            return 0;
        }
        if self.left() >= other.right() {
            return 0;
        }
        if self.left() <= other.left() {
            return self.get_overlap_(other);
        }
        other.get_overlap_(self)
    }

    fn get_overlap_(&self, other: &Cds) -> usize {
        assert!(self.left() <= other.left());
        assert!(other.left() <= self.right());
        if other.right() <= self.right() {
            return other.size();
        }
        self.right() - other.left()
    }

    pub fn coexists(&self, next: &Cds) -> bool {
        assert!(self.left() <= next.left());
        if !self.strand() && next.strand() {
            return self.right() + Self::PROMOTER_MIN <= next.left();
        }
        let overlap_min =
            (Self::SIZE_MIN - 3).min((self.size().min(next.size()) as f64 * 0.2) as usize);
        self.get_overlap(next) <= overlap_min
    }

    pub fn get_len_increase(&self, prev: &Cds) -> isize {
        assert!(prev.left() <= self.left());
        self.size_effective()
    }

    pub fn worse(&self, than: &Cds) -> bool {
        than.left() == self.left()
            && than.size() == self.size()
            && (than.positives > self.positives
                || (than.positives == self.positives && than.ref_prot < self.ref_prot))
    }
}

impl Eq for Cds {}

impl PartialOrd for Cds {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cds {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.left()
            .cmp(&other.left())
            .then(other.right().cmp(&self.right()))
            .then_with(|| {
                other
                    .positives
                    .partial_cmp(&self.positives)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }
}

#[derive(Debug, Clone, Default)]
pub struct DnaAnnot {
    pub cdss: Vec<Cds>,
}

impl DnaAnnot {
    pub fn new() -> Self {
        DnaAnnot {
            cdss: Vec::with_capacity(10000),
        }
    }

    pub fn run(&mut self) -> Option<usize> {
        self.cdss.sort();

        let mut nexts: Vec<usize> = Vec::new();
        for i in 0..self.cdss.len() {
            if i > 0 {
                assert!(self.cdss[i].left() >= self.cdss[i - 1].left());
            }

            assert!(self.cdss[i].best_prev.is_none());
            self.cdss[i].sum_len = self.cdss[i].size_effective();
            let prevs = self.cdss[i].prevs.clone();
            for prev in prevs {
                let sum_len_new =
                    self.cdss[prev].sum_len + self.cdss[i].get_len_increase(&self.cdss[prev]);
                if self.cdss[i].sum_len < sum_len_new
                    || (self.cdss[i].sum_len == sum_len_new
                        && self.cdss[i]
                            .best_prev
                            .is_some_and(|best_prev| prev < best_prev))
                {
                    self.cdss[i].sum_len = sum_len_new;
                    self.cdss[i].best_prev = Some(prev);
                }
            }
            if !self.cdss[i].prevs.is_empty() {
                assert!(self.cdss[i].best_prev.is_some() && self.cdss[i].sum_len > 0);
            }

            nexts.clear();
            for j in i + 1..self.cdss.len() {
                if !self.cdss[i].coexists(&self.cdss[j]) {
                    continue;
                }

                let mut too_far = false;
                let next_len_increase = self.cdss[j].get_len_increase(&self.cdss[i]);
                for &next_old in &nexts {
                    if self.cdss[next_old].coexists(&self.cdss[j])
                        && self.cdss[j].get_len_increase(&self.cdss[next_old])
                            + self.cdss[next_old].get_len_increase(&self.cdss[i])
                            > next_len_increase
                    {
                        too_far = true;
                        break;
                    }
                }
                if too_far {
                    break;
                }

                nexts.push(j);
                self.cdss[j].prevs.push(i);
            }
        }

        let mut cds_last = None;
        let mut sum_len_max = 0;
        for (i, cds) in self.cdss.iter().enumerate() {
            if sum_len_max < cds.sum_len {
                sum_len_max = cds.sum_len;
                cds_last = Some(i);
            }
        }

        cds_last
    }
}

#[derive(Debug, Clone, Default, Eq)]
pub struct Mutation {
    pub prot: bool,
    pub gene_name: String,
    pub pos: usize,
    pub ref_: String,
    pub allele: String,
    pub frameshift: bool,
    pub ambig: bool,
}

impl Mutation {
    pub fn from_line(prot_arg: bool, line: &str) -> anyhow::Result<Self> {
        let mut this = Mutation {
            prot: prot_arg,
            pos: NO_INDEX,
            ..Mutation::default()
        };
        let s = line.trim();
        let dash_pos = s.find('-');

        if let Some(dash_pos) = dash_pos {
            this.gene_name = s[..dash_pos].to_string();
            assert!(!this.gene_name.is_empty());
        }

        let ref_start = dash_pos.map_or(0, |pos| pos + 1);
        let mut pos_start = NO_INDEX;
        let mut allele_start = NO_INDEX;
        for (i, c) in s.char_indices().skip_while(|(i, _)| *i < ref_start) {
            if c.is_ascii_digit() {
                if pos_start == NO_INDEX {
                    pos_start = i;
                }
            } else if pos_start != NO_INDEX && allele_start == NO_INDEX {
                allele_start = i;
                break;
            }
        }
        assert!(ref_start < pos_start);
        assert!(pos_start < allele_start);
        assert!(allele_start < s.len());

        this.ref_ = s[ref_start..pos_start].to_string();
        this.pos = s[pos_start..allele_start].parse::<usize>()?;
        this.allele = s[allele_start..].to_string();
        assert!(this.pos > 0);
        this.pos -= 1;

        if this.prot {
            if this.ref_ == "ins" {
                this.ref_.clear();
            }
            if this.allele == "del" {
                this.allele.clear();
            }
            if this.allele == "fs" {
                this.frameshift = true;
                this.allele.clear();
            }
        } else {
            if this.ref_ == "INS" {
                this.ref_.clear();
            }
            if this.allele == "DEL" {
                this.allele.clear();
            }
        }

        this.set_ambig();
        Ok(this)
    }

    pub fn new(gene_name_arg: String, pos_arg: usize, ref_arg: String, allele_arg: String) -> Self {
        let mut this = Mutation {
            prot: !gene_name_arg.is_empty(),
            gene_name: gene_name_arg,
            pos: pos_arg,
            ref_: ref_arg,
            allele: allele_arg,
            frameshift: false,
            ambig: false,
        };
        this.set_ambig();
        this
    }

    pub fn new_frameshift(
        gene_name_arg: String,
        pos_arg: usize,
        ref_arg: String,
        allele_arg: String,
        frameshift_arg: bool,
    ) -> Self {
        assert!(!gene_name_arg.is_empty());
        let mut this = Mutation {
            prot: true,
            gene_name: gene_name_arg,
            pos: pos_arg,
            ref_: ref_arg,
            allele: allele_arg,
            frameshift: frameshift_arg,
            ambig: false,
        };
        this.set_ambig();
        this
    }

    fn set_ambig(&mut self) {
        self.ambig = false;
        for c in self.allele.chars() {
            if is_ambig(c, self.prot) {
                self.ambig = true;
                return;
            }
        }
    }

    pub fn qc(&self) -> anyhow::Result<()> {
        if !self.gene_name.is_empty() {
            assert!(self
                .gene_name
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '-'));
        }
        assert_ne!(self.pos, NO_INDEX);

        if self.prot {
            for c in self.ref_.chars() {
                if c != TERMINATOR.chars().next().unwrap() && !PEPTIDE_ALPHABET.contains(c) {
                    anyhow::bail!(
                        "Protein mutation cannot have ambiguities in the reference sequence"
                    );
                }
            }
            for c in self.allele.chars() {
                assert!(EXT_TERM_PEPTIDE_ALPHABET.contains(c));
            }
        } else {
            for c in self.ref_.chars() {
                if !DNA_ALPHABET.contains(c) {
                    anyhow::bail!("DNA mutation cannot have ambiguities in the reference sequence");
                }
            }
            for c in self.allele.chars() {
                assert!(EXT_DNA_ALPHABET.contains(c));
            }
        }

        if !self.frameshift {
            assert_ne!(self.ref_, self.allele);
        }
        if self.frameshift {
            assert!(self.prot && self.allele.is_empty());
        }
        Ok(())
    }

    pub fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        if self.prot {
            write!(
                os,
                "{}-{}{}{}",
                self.gene_name,
                if self.ref_.is_empty() {
                    "ins"
                } else {
                    &self.ref_
                },
                self.pos + 1,
                if self.frameshift {
                    "fs"
                } else if self.allele.is_empty() {
                    "del"
                } else {
                    &self.allele
                }
            )?;
        } else {
            write!(
                os,
                "{}{}{}",
                if self.ref_.is_empty() {
                    "INS"
                } else {
                    &self.ref_
                },
                self.pos + 1,
                if self.allele.is_empty() {
                    "DEL"
                } else {
                    &self.allele
                }
            )?;
        }
        Ok(())
    }

    pub fn stop(&self) -> usize {
        self.pos + self.ref_.len()
    }

    pub fn replace(&self, ref_dna: &mut Dna) {
        assert!(!self.prot);
        assert!(self.stop() <= ref_dna.seq.len());
        assert_eq!(
            &ref_dna.seq[self.pos..self.pos + self.ref_.len()],
            self.ref_
        );
        ref_dna
            .seq
            .replace_range(self.pos..self.pos + self.ref_.len(), &self.allele);
    }
}

impl PartialEq for Mutation {
    fn eq(&self, other: &Self) -> bool {
        self.prot == other.prot
            && self.gene_name == other.gene_name
            && self.pos == other.pos
            && self.ref_ == other.ref_
            && self.allele == other.allele
            && self.frameshift == other.frameshift
    }
}

impl PartialOrd for Mutation {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Mutation {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.prot
            .cmp(&other.prot)
            .then(self.gene_name.cmp(&other.gene_name))
            .then(self.pos.cmp(&other.pos))
            .then(self.ref_.cmp(&other.ref_))
            .then(self.allele.cmp(&other.allele))
            .then(self.frameshift.cmp(&other.frameshift))
    }
}

impl Hash for Mutation {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.prot.hash(state);
        self.gene_name.hash(state);
        self.pos.hash(state);
        self.ref_.hash(state);
        self.allele.hash(state);
        self.frameshift.hash(state);
    }
}

pub type AlignScore = f64;

#[derive(Debug, Clone)]
pub struct SubstMat {
    pub sim: [[AlignScore; Self::SIM_SIZE]; Self::SIM_SIZE],
}

impl Default for SubstMat {
    fn default() -> Self {
        SubstMat {
            sim: [[AlignScore::NAN; Self::SIM_SIZE]; Self::SIM_SIZE],
        }
    }
}

impl SubstMat {
    pub const CHARS: &'static str = "ARNDCQEGHILKMFPSTWYVBZX*";
    pub const SIM_SIZE: usize = 128;

    pub fn new(f_name: &str) -> anyhow::Result<Self> {
        for c in Self::CHARS.chars() {
            assert!(EXT_TERM_PEPTIDE_ALPHABET.contains(c));
        }

        let mut mat = SubstMat::default();
        let mut char_num = 0usize;
        for line in std::fs::read_to_string(f_name)?.lines() {
            let c1 = line
                .chars()
                .next()
                .ok_or_else(|| anyhow::anyhow!("Empty substitution matrix line"))?;
            if c1 == '#' {
                continue;
            }
            if c1 == ' ' {
                assert_eq!(char_num, 0);
                let mut header = line.to_string();
                header.retain(|c| c != ' ');
                assert_eq!(header, Self::CHARS);
            } else {
                assert!(char_num < Self::CHARS.len());
                assert_eq!(c1, Self::CHARS.as_bytes()[char_num] as char);

                let row = c1 as usize;
                assert!(row < Self::SIM_SIZE);
                let mut fields = line.split_whitespace();
                assert_eq!(fields.next(), Some(&line[..c1.len_utf8()]));
                for c in Self::CHARS.chars() {
                    let col = c as usize;
                    assert!(col < Self::SIM_SIZE);
                    let num_s = fields
                        .next()
                        .ok_or_else(|| anyhow::anyhow!("Missing substitution score"))?;
                    let r: AlignScore = num_s.parse()?;
                    assert!(!r.is_nan());
                    mat.sim[row][col] = r;
                }
                char_num += 1;
            }
        }
        assert_eq!(char_num, Self::CHARS.len());
        Ok(mat)
    }

    pub fn qc(&self) {
        for i in 0..Self::SIM_SIZE {
            if !self.good_index(i) {
                continue;
            }
            if !is_ambig_aa(i as u8 as char) {
                assert!(self.sim[i][i] >= 0.0);
            }
            for j in 0..Self::SIM_SIZE {
                if !self.good_index(j) {
                    continue;
                }
                assert!((self.sim[i][j] - self.sim[j][i]).abs() < 1e-6);
            }
        }
    }

    pub fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        for i in 0..Self::SIM_SIZE {
            if self.good_index(i) {
                write!(os, " {}", i as u8 as char)?;
            }
        }
        writeln!(os)?;

        for i in 0..Self::SIM_SIZE {
            if self.good_index(i) {
                write!(os, "{}", i as u8 as char)?;
                for j in 0..Self::SIM_SIZE {
                    if self.good_index(j) {
                        write!(os, " {}", self.sim[i][j])?;
                    }
                }
                writeln!(os)?;
            }
        }
        Ok(())
    }

    pub fn good_index(&self, i: usize) -> bool {
        self.sim[i][i] == self.sim[i][i]
    }

    pub fn print_anomalies(&self) {
        for i in 0..Self::SIM_SIZE {
            if !self.good_index(i) {
                continue;
            }
            for j in 0..Self::SIM_SIZE {
                if !self.good_index(j) {
                    continue;
                }
                if self.sim[i][i] < self.sim[i][j] {
                    println!(
                        "sim({},{}) < sim({},{})",
                        i as u8 as char, i as u8 as char, i as u8 as char, j as u8 as char
                    );
                }
            }
        }
    }

    pub fn get_substitution_dist(&self, row: usize, col: usize) -> AlignScore {
        self.sim[row][row] + self.sim[col][col] - 2.0 * self.sim[row][col]
    }

    pub fn get_deletion_dist(&self, row: usize, gap_sim: AlignScore) -> AlignScore {
        self.sim[row][row] - 2.0 * gap_sim
    }

    pub fn char2score(&self, mut c1: char, mut c2: char) -> anyhow::Result<AlignScore> {
        assert!((c1 as u32) < Self::SIM_SIZE as u32);
        assert!((c2 as u32) < Self::SIM_SIZE as u32);

        if PEPTIDE_WILDCARDS.contains(c1) {
            c1 = 'X';
        }
        if PEPTIDE_WILDCARDS.contains(c2) {
            c2 = 'X';
        }

        let i1 = c1 as usize;
        let i2 = c2 as usize;
        assert!(i1 < Self::SIM_SIZE);
        assert!(i2 < Self::SIM_SIZE);

        if c1 == '*' || c2 == '*' {
            return Ok(-10.0);
        }
        if c1 == '-' || c2 == '-' {
            return Ok(-1.0);
        }

        if !self.good_index(i1) {
            anyhow::bail!("Bad amino acid: {} ({})", c1, i1);
        }
        if !self.good_index(i2) {
            anyhow::bail!("Bad amino acid: {} ({})", c2, i2);
        }

        let s = self.sim[i1][i2];
        assert!(!s.is_nan());
        Ok(s)
    }

    pub fn char2score_static(
        sm: Option<&SubstMat>,
        c1: char,
        c2: char,
    ) -> anyhow::Result<AlignScore> {
        Ok(if let Some(sm) = sm {
            sm.char2score(c1, c2)?
        } else {
            aa_match(c1, c2) as u8 as AlignScore
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PeptideOrf {
    pub translation_start: usize,
    pub strand: Strand,
    pub start: usize,
    pub start_m: bool,
    pub stop: usize,
    pub stop_terminator: bool,
}

impl PeptideOrf {
    pub fn new(translation_start: usize, strand: Strand, peptide: &Peptide, start: usize) -> Self {
        assert!(is_strand(strand));
        assert!(!peptide.seq.is_empty());
        assert!(!peptide.sparse);
        assert!(start < peptide.seq.len());

        let mut orf = PeptideOrf {
            translation_start,
            strand,
            start,
            start_m: Peptide::is_start_aa(peptide.seq.as_bytes()[start] as char),
            stop: peptide.seq.len(),
            stop_terminator: false,
        };

        for i in start..peptide.seq.len() {
            if peptide.seq.as_bytes()[i] as char == TERMINATOR.chars().next().unwrap() {
                orf.stop = i;
                orf.stop_terminator = true;
                break;
            }
        }
        assert!(orf.stop <= peptide.seq.len());
        assert!(!orf.empty());
        orf
    }

    pub fn new_fields(
        translation_start: usize,
        strand: Strand,
        start: usize,
        start_m: bool,
        stop: usize,
        stop_terminator: bool,
    ) -> Self {
        PeptideOrf {
            translation_start,
            strand,
            start,
            start_m,
            stop,
            stop_terminator,
        }
    }

    pub fn from_text(s: &str) -> anyhow::Result<Self> {
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() != 6 {
            anyhow::bail!("PeptideOrf: expected 6 fields");
        }
        let translation_start = fields[0].parse()?;
        let strand = fields[1].parse::<i8>()?;
        let start = fields[2].parse()?;
        let stop = fields[3].parse()?;
        let start_m = fields[4].parse::<usize>()? != 0;
        let stop_terminator = fields[5].parse::<usize>()? != 0;
        Ok(PeptideOrf {
            translation_start,
            strand,
            start,
            start_m,
            stop,
            stop_terminator,
        })
    }

    pub fn qc(&self) {
        if self.empty() {
            return;
        }

        assert!(is_strand(self.strand));
        if self.strand == 0 {
            assert_eq!(self.translation_start, 0);
        }
        assert!(self.start <= self.stop);
        if self.start == self.stop {
            assert!(!self.start_m && self.stop_terminator);
        }
        if self.start != 0 {
            assert!(self.start_m);
        }
    }

    pub fn save_text<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.translation_start,
            self.strand,
            self.start,
            self.stop,
            self.start_m as u8,
            self.stop_terminator as u8
        )
    }

    pub fn empty(&self) -> bool {
        self.start == NO_INDEX
    }

    pub fn good(&self, size_min: usize) -> bool {
        !self.empty() && self.start_m && self.stop_terminator && self.size() >= size_min
    }

    pub fn size(&self) -> usize {
        self.stop - self.start
    }

    pub fn dna_pos(&self, pos: usize) -> usize {
        let dna_len = 3 * pos;
        if self.strand == 1 {
            self.translation_start + dna_len
        } else {
            self.translation_start - dna_len
        }
    }

    pub fn cds_start(&self) -> usize {
        self.dna_pos(self.start)
    }

    pub fn cds_stop(&self) -> usize {
        self.dna_pos(self.stop + self.stop_terminator as usize)
    }

    pub fn to_peptide(&self, peptide: &Peptide) -> Peptide {
        assert!(!peptide.sparse);
        assert!(!self.empty());
        assert!(self.stop <= peptide.seq.len());
        if self.stop < peptide.seq.len() {
            assert!(self.stop_terminator);
            assert_eq!(
                peptide.seq.as_bytes()[self.stop] as char,
                TERMINATOR.chars().next().unwrap()
            );
        } else {
            assert!(!self.stop_terminator);
        }

        let mut p = Peptide::new(&peptide.name, &peptide.seq, false);
        p.seq.truncate(self.stop + self.stop_terminator as usize);
        p.seq.drain(0..self.start);
        p
    }
}

impl Default for PeptideOrf {
    fn default() -> Self {
        PeptideOrf {
            translation_start: NO_INDEX,
            strand: 0,
            start: NO_INDEX,
            start_m: false,
            stop: NO_INDEX,
            stop_terminator: false,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Peptide {
    pub name: String,
    pub seq: String,
    pub sparse: bool,
    pub pseudo: bool,
}

impl Peptide {
    pub const STD_AVE_LEN: usize = 400;
    pub const STD_MIN_COMPLEXITY: f64 = 2.5;

    pub fn new(name: &str, seq: &str, sparse: bool) -> Self {
        let base = Seq::new(name, seq, sparse);
        Peptide {
            name: base.name,
            seq: base.seq,
            sparse: base.sparse,
            pseudo: false,
        }
    }

    pub fn with_len(name: &str, seq_len: usize, sparse: bool) -> Self {
        let base = Seq::with_len(name, seq_len, sparse);
        Peptide {
            name: base.name,
            seq: base.seq,
            sparse: base.sparse,
            pseudo: false,
        }
    }

    pub fn from_multifasta(
        fasta: &mut Multifasta,
        reserve_len: usize,
        sparse: bool,
    ) -> anyhow::Result<Self> {
        assert!(fasta.aa);
        let base = Seq::from_multifasta(fasta, reserve_len, sparse, true)?;
        let this = Peptide {
            name: base.name,
            seq: base.seq,
            sparse: base.sparse,
            pseudo: false,
        };
        this.qc()?;
        Ok(this)
    }

    pub fn get_seq_alphabet(&self) -> &'static str {
        if self.sparse {
            EXT_SPARSE_TERM_PEPTIDE_ALPHABET
        } else {
            EXT_TERM_PEPTIDE_ALPHABET
        }
    }

    pub fn qc(&self) -> anyhow::Result<()> {
        let base = Seq {
            name: self.name.clone(),
            seq: self.seq.clone(),
            sparse: self.sparse,
        };
        base.qc_name()?;
        base.qc_alphabet(self.get_seq_alphabet())?;
        anyhow::ensure!(
            !self.has_inside_stop() || self.pseudo,
            "Peptide has an internal stop codon"
        );
        Ok(())
    }

    pub fn is_ambiguous(c: char) -> bool {
        c == 'X'
    }

    pub fn has_inside_stop(&self) -> bool {
        let pos = self.seq.find('*');
        pos.is_some_and(|pos| pos != self.seq.len() - 1)
    }

    pub fn trim_stop(&mut self) -> bool {
        if self.seq.ends_with('*') {
            self.seq.pop();
            true
        } else {
            false
        }
    }

    pub fn is_description_partial(&self) -> bool {
        let desc = Seq {
            name: self.name.clone(),
            seq: self.seq.clone(),
            sparse: self.sparse,
        }
        .get_description(false)
        .to_ascii_lowercase();
        desc.contains("fragmented")
            || desc.contains("fragmentary")
            || desc.contains("partial")
            || desc.contains("truncat")
    }

    pub fn is_start_aa(aa: char) -> bool {
        aa != '*' && aa.is_ascii_lowercase()
    }

    pub fn get_peptide_orfs(
        &self,
        translation_start: usize,
        strand: Strand,
        include_initial: bool,
        longest_only: bool,
        len_min: usize,
    ) -> Vec<PeptideOrf> {
        let mut orfs = Vec::new();

        let mut prev = PeptideOrf::default();
        for i in 0..self.seq.len() {
            if Peptide::is_start_aa(self.seq.as_bytes()[i] as char) {
                let orf = PeptideOrf::new(translation_start, strand, self, i);
                assert!(!orf.empty());
                orf.qc();
                if !prev.empty() {
                    assert!(prev.start < orf.start);
                    if longest_only && prev.stop == orf.stop {
                        continue;
                    }
                }
                if orf.size() >= len_min {
                    orfs.push(orf.clone());
                }
                prev = orf;
            }
        }

        if include_initial && (orfs.is_empty() || orfs[0].start > 0) {
            let orf = PeptideOrf::new(translation_start, strand, self, 0);
            assert!(!orf.empty());
            orf.qc();
            if orf.size() >= len_min {
                orfs.push(orf);
            }
        }

        orfs
    }

    pub fn get_complexity_int(&self, start: usize, mut end: usize) -> f64 {
        if self.seq.is_empty() {
            return 0.0;
        }

        const MAX_AA_NUM: usize = 21;
        let mut di_freq = [[0usize; MAX_AA_NUM]; MAX_AA_NUM];
        let mut n = 0usize;
        end = end.min(self.seq.len());
        let bytes = self.seq.as_bytes();
        for i in start + 1..end {
            let a = aa2num(bytes[i - 1] as char);
            let b = aa2num(bytes[i] as char);
            if a < MAX_AA_NUM && b < MAX_AA_NUM {
                di_freq[a][b] += 1;
                n += 1;
            }
        }

        let mut s = 0.0;
        for row in di_freq {
            for count in row {
                if count != 0 {
                    let p = count as f64 / n as f64;
                    s -= p * p.ln();
                }
            }
        }
        assert!(s >= 0.0);
        s
    }

    pub fn get_complexity(&self) -> f64 {
        self.get_complexity_int(0, self.seq.len())
    }

    pub fn get_self_similarity(&self, mat: &SubstMat, start: usize, mut stop: usize) -> f64 {
        assert!(start <= stop);
        assert!(stop <= self.seq.len());

        if start == 0 && stop == 0 {
            stop = self.seq.len();
        }

        let mut s = 0.0;
        for c in self.seq[start..stop].chars() {
            if c != '-' {
                s += mat.sim[c as usize][c as usize];
            }
        }
        s
    }

    pub fn get_similarity(
        &self,
        other: &Peptide,
        mat: &SubstMat,
        gap_open_cost: f64,
        gap_cost: f64,
    ) -> f64 {
        assert!(gap_open_cost >= 0.0);
        assert!(gap_cost >= 0.0);
        assert!(self.sparse);
        assert!(other.sparse);
        assert_eq!(self.seq.len(), other.seq.len());

        let mut s = 0.0;
        let mut c1_prev = '\0';
        let mut c2_prev = '\0';
        for (c1, c2) in self.seq.chars().zip(other.seq.chars()) {
            if c1 == '-' && c2 == '-' {
            } else if c1 == '-' || c2 == '-' {
                s -= gap_cost;
                if (c1 == '-' && c1_prev != '-') || (c2 == '-' && c2_prev != '-') {
                    s -= gap_open_cost;
                }
            } else {
                s += mat.sim[c1 as usize][c2 as usize];
            }
            c1_prev = c1;
            c2_prev = c2;
        }

        s
    }

    pub fn ambig2x(&mut self) -> usize {
        let mut n = 0usize;
        let mut seq = String::with_capacity(self.seq.len());
        for c in self.seq.chars() {
            if "BZJUO".contains(c) {
                seq.push('X');
                n += 1;
            } else {
                seq.push(c);
            }
        }
        self.seq = seq;
        n
    }

    pub fn to_gbmr4(&mut self) {
        let mut seq = String::with_capacity(self.seq.len());
        for c in self.seq.chars() {
            seq.push(match c {
                'A' | 'D' | 'K' | 'E' | 'R' | 'N' | 'T' | 'S' | 'Q' => 'A',
                'Y' | 'F' | 'L' | 'I' | 'V' | 'M' | 'C' | 'W' | 'H' => 'Y',
                'G' | 'P' | 'X' | '*' => c,
                _ => panic!("Unknown aa {c}"),
            });
        }
        self.seq = seq;
    }
}

// --- CDS constants ---

/// Minimum peptide size for a CDS
pub const CDS_PEPTIDE_SIZE_MIN: usize = 20;

// --- Interval ---

/// Represents a genomic interval (0-based, half-open)
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Interval {
    pub start: usize, // NO_INDEX if empty
    pub stop: usize,  // position after segment
    pub strand: Strand,
}

impl Interval {
    pub fn new(start: usize, stop: usize, strand: Strand) -> Self {
        Interval {
            start,
            stop,
            strand,
        }
    }

    pub fn empty(&self) -> bool {
        self.start == NO_INDEX
    }

    pub fn valid(&self) -> bool {
        self.start <= self.stop
    }

    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.stop - self.start
    }

    pub fn frame(&self) -> Frame {
        self.strand * ((self.start % 3) as i8 + 1)
    }

    pub fn inside_eq(&self, other: &Interval, slack: usize) -> bool {
        self.strand == other.strand
            && self.start + slack >= other.start
            && self.stop <= other.stop + slack
    }

    pub fn contains(&self, other: &Interval) -> bool {
        other.inside_eq(self, 0)
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        self.strand == other.strand && self.stop > other.start && self.start < other.stop
    }

    pub fn contains_strongly(&self, other: &Interval) -> bool {
        self.contains(other) && self.overlaps(other)
    }

    pub fn rest(&self, seq_len: usize, upstream: bool) -> usize {
        if (self.strand == -1) != upstream {
            self.start
        } else {
            seq_len - self.stop
        }
    }

    /// Format as 1-based for display
    pub fn format_display(&self) -> String {
        let mut text = if self.strand == -1 {
            format!("{}-{}", self.stop, self.start + 1)
        } else {
            format!("{}-{}", self.start + 1, self.stop)
        };
        if self.strand != 0 {
            text.push('/');
            text.push(strand2char(self.strand));
        }
        text
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.strand
            .cmp(&other.strand)
            .then(self.start.cmp(&other.start))
            .then(self.stop.cmp(&other.stop))
    }
}

// --- Disruption ---

/// Type of gene disruption
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DisruptionType {
    None,
    Smooth,
    Frameshift,
    Deletion, // or replacement
    Insertion,
}

impl DisruptionType {
    pub const NAMES: &[&str] = &["none", "smooth", "fs", "del", "ins"];

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "none" => Some(DisruptionType::None),
            "smooth" => Some(DisruptionType::Smooth),
            "fs" => Some(DisruptionType::Frameshift),
            "del" => Some(DisruptionType::Deletion),
            "ins" => Some(DisruptionType::Insertion),
            _ => None,
        }
    }
}

pub const STOP_SUFFIX: &str = "_STOP";

/// Gene disruption (between HSPs in a blastx merge)
#[derive(Debug, Clone)]
pub struct Disruption {
    pub prev_hsp_idx: Option<usize>, // index into Hsp vector
    pub next_hsp_idx: Option<usize>,
    pub prev_start: usize, // position in prev qseq/sseq
    pub next_stop: usize,  // position in next qseq/sseq
    pub prev_slen: usize,
    pub intron: bool,
    pub s_stop_codon: bool,
    // Cached intervals (computed from HSP data)
    pub q_interval: Interval,
    pub s_interval: Interval,
}

impl Disruption {
    pub fn empty(&self) -> bool {
        self.prev_hsp_idx.is_none()
    }

    pub fn disruption_type(&self) -> DisruptionType {
        if self.empty() {
            return DisruptionType::None;
        }
        if self.q_interval.len() == 0 && self.s_interval.len() == 0 {
            return DisruptionType::Smooth;
        }
        if !self.s_interval.len().is_multiple_of(3) {
            return DisruptionType::Frameshift;
        }
        if self.q_interval.len() > 0 {
            return DisruptionType::Deletion;
        }
        DisruptionType::Insertion
    }

    pub fn get_len(&self) -> usize {
        std::cmp::max(self.q_interval.len(), self.s_interval.len() / 3)
    }

    pub fn reportable(&self) -> bool {
        !matches!(
            self.disruption_type(),
            DisruptionType::None | DisruptionType::Smooth
        )
    }

    pub fn qc(&self) {
        if self.empty() {
            assert!(self.prev_hsp_idx.is_none());
            assert!(self.next_hsp_idx.is_none());
            assert_eq!(self.prev_start, 0);
            assert_eq!(self.next_stop, 0);
            assert_eq!(self.prev_slen, NO_INDEX);
            assert!(!self.intron);
            assert!(!self.s_stop_codon);
            return;
        }

        assert!(self.prev_hsp_idx.is_some());
        assert!(self.next_hsp_idx.is_some());
        assert!(self.q_interval.valid());
        assert!(self.s_interval.valid());
        match self.disruption_type() {
            DisruptionType::Frameshift => {
                assert!(self.intron);
            }
            DisruptionType::Insertion => {
                assert!(self.prev_start > 0);
                assert!(self.q_interval.start > 0);
            }
            _ => {}
        }
        if self.prev_hsp_idx == self.next_hsp_idx && self.prev_start == self.next_stop {
            assert_eq!(self.disruption_type(), DisruptionType::Smooth);
        }
    }

    pub fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        write!(os, "Disruption:")?;
        if self.empty() {
            write!(os, "empty")?;
        } else {
            write!(
                os,
                "{}:{}",
                self.q_interval.format_display(),
                self.s_interval.format_display()
            )?;
            if self.s_stop_codon() {
                write!(os, "/*")?;
            }
            if self.intron {
                write!(os, "/intron")?;
            }
        }
        Ok(())
    }

    pub fn s_stop_codon(&self) -> bool {
        assert!(!self.empty());
        self.prev_hsp_idx == self.next_hsp_idx && !self.intron && self.s_stop_codon
    }

    pub fn genesymbol_raw(&self) -> String {
        assert!(!self.empty());

        let dtype = self.disruption_type();
        assert_ne!(dtype, DisruptionType::Smooth);

        let mut q_interval = self.q_interval.clone();
        let mut s_interval = self.s_interval.clone();
        match dtype {
            DisruptionType::Frameshift => {
                q_interval.stop = q_interval.start + 1;
                if s_interval.strand == 1 {
                    s_interval.stop = self.prev_slen;
                } else {
                    s_interval.start = 0;
                }
            }
            DisruptionType::Insertion => {
                assert_eq!(q_interval.len(), 0);
                q_interval.start -= 1;
                if s_interval.strand == -1 {
                    s_interval.stop += 3;
                } else {
                    s_interval.start -= 3;
                }
            }
            _ => {}
        }

        let mut s = format!(
            "{}_{}_{}_{}_{}_{}",
            DisruptionType::NAMES[dtype as usize],
            q_interval.start,
            q_interval.stop,
            s_interval.start,
            s_interval.stop,
            if s_interval.strand == 1 { 1 } else { 0 }
        );
        if self.s_stop_codon() {
            s.push_str(STOP_SUFFIX);
        }
        s
    }
}

impl PartialEq for Disruption {
    fn eq(&self, other: &Self) -> bool {
        self.q_interval == other.q_interval && self.s_interval == other.s_interval
    }
}

impl Eq for Disruption {}

impl PartialOrd for Disruption {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Disruption {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        assert!(!self.empty());
        assert!(!other.empty());
        assert_eq!(self.s_interval.strand, other.s_interval.strand);
        self.q_interval
            .cmp(&other.q_interval)
            .then(self.s_interval.cmp(&other.s_interval))
    }
}

// --- Hsp ---

/// BLAST High Scoring Pair
///
/// Field naming convention (matching C++ AMRFinderPlus):
/// - **q (query)** = the reference protein from the AMR database
/// - **s (subject)** = the user's input sequence (protein or DNA)
///
/// This is the reverse of BLAST's own convention when using `-outfmt '6 sseqid qseqid ...'`
/// (reverse format). The BLAST output columns are swapped so that when parsed left-to-right,
/// field 1 becomes `qseqid` (reference) and field 2 becomes `sseqid` (user input).
///
/// - `qseqid` contains the reference accession (pipe-delimited metadata)
/// - `sseqid` contains the user's protein/contig identifier
/// - `qlen` = reference sequence length
/// - `slen` = user's sequence length
/// - `q_int` = alignment interval on the reference
/// - `s_int` = alignment interval on the user's sequence
#[derive(Debug, Clone)]
pub struct Hsp {
    pub merged: bool,

    // Protein flags
    pub q_prot: bool, // reference is protein
    pub s_prot: bool, // user input is protein (false for blastx/tblastn DNA)
    pub a_prot: bool, // alignment is in protein space

    // Unit conversion factors (1 or 3 for protein-to-DNA)
    pub a2q: usize,
    pub a2s: usize,

    // BLAST fields — see struct doc for q/s convention
    pub qseqid: String,  // reference identifier (pipe-delimited metadata)
    pub sseqid: String,  // user's sequence identifier
    pub q_int: Interval, // alignment interval on reference
    pub s_int: Interval, // alignment interval on user's sequence
    pub qlen: usize,     // reference sequence length
    pub slen: usize,     // user's sequence length
    pub qseq: String,    // aligned reference sequence (with gaps)
    pub sseq: String,    // aligned user sequence (with gaps)

    // Computed by finish_hsp()
    pub length: usize,
    pub nident: usize,
    pub qgap: usize,
    pub sgap: usize,
    pub qx: usize,
    pub sx: usize,

    pub(crate) pos2q_: Vec<usize>,
    pub(crate) pos2s_: Vec<usize>,

    pub sframe: Frame,
    pub c_complete: Option<bool>, // None = unknown, Some(true) = detected, Some(false) = missing
    pub s_internal_stop: bool,

    pub disrs: Vec<Disruption>,
}

impl Hsp {
    pub const BLASTP_FAST: &str = "  -comp_based_stats 0  -seg no  -max_target_seqs 10000  -dbsize 10000  -evalue 1e-10  -word_size 5";
    pub const BLASTP_SLOW: &str = "  -comp_based_stats 0  -seg no  -max_target_seqs 10000  -dbsize 10000  -evalue 1      -word_size 3";

    pub const FORMAT: [&str; 2] = [
        "sseqid qseqid sstart send slen qstart qend qlen sseq qseq",
        "qseqid sseqid qstart qend qlen sstart send slen qseq sseq",
    ];

    pub fn format_par(forward: bool) -> String {
        format!("  -outfmt '6 {}'", Self::FORMAT[forward as usize])
    }

    pub fn empty(&self) -> bool {
        self.sseqid.is_empty()
    }

    pub fn blastx(&self) -> bool {
        self.q_prot && !self.s_prot
    }

    pub fn q_abs_coverage(&self) -> usize {
        self.q_int.len()
    }

    pub fn q_effective_len(&self) -> usize {
        self.qlen
    }

    pub fn s_abs_coverage(&self) -> usize {
        self.s_int.len()
    }

    pub fn q_len(&self) -> usize {
        assert!(self.disrs.is_empty());
        self.q_abs_coverage() / self.a2q
    }

    pub fn s_len(&self) -> usize {
        assert!(self.disrs.is_empty());
        self.s_abs_coverage() / self.a2s
    }

    pub fn rel_identity(&self) -> f64 {
        self.nident as f64 / self.length as f64
    }

    pub fn q_rel_coverage(&self) -> f64 {
        self.q_abs_coverage() as f64 / self.q_effective_len() as f64
    }

    pub fn s_rel_coverage(&self) -> f64 {
        self.s_abs_coverage() as f64 / self.slen as f64
    }

    pub fn q_complete(&self) -> bool {
        self.q_abs_coverage() == self.q_effective_len() && self.c_complete != Some(false)
    }

    pub fn perfect(&self) -> bool {
        self.q_complete() && self.disrs.is_empty()
    }

    pub fn s_truncated(&self) -> bool {
        (self.s_int.start < self.a2s
            && ((self.s_int.strand == 1 && self.q_int.start > 0)
                || (self.s_int.strand == -1 && self.q_int.stop < self.qlen)))
            || (self.slen - self.s_int.stop < self.a2s
                && ((self.s_int.strand == 1 && self.q_int.stop < self.qlen)
                    || (self.s_int.strand == -1 && self.q_int.start > 0)))
    }

    pub fn s_inside_eq(&self, other: &Hsp, slack: usize) -> bool {
        self.s_int.inside_eq(&other.s_int, slack)
    }

    pub fn char_match(&self, pos: usize) -> bool {
        if self.a_prot {
            aa_match(
                self.qseq.as_bytes()[pos] as char,
                self.sseq.as_bytes()[pos] as char,
            )
        } else {
            nucleotide_match(
                self.qseq.as_bytes()[pos] as char,
                self.sseq.as_bytes()[pos] as char,
            )
        }
    }

    pub fn find_disruption(&self, dtype: DisruptionType) -> Option<&Disruption> {
        self.disrs.iter().find(|d| d.disruption_type() == dtype)
    }

    pub fn has_reportable_disruption(&self) -> bool {
        self.disrs.iter().any(Disruption::reportable)
    }

    pub fn reportable_disruption_count(&self) -> usize {
        self.disrs
            .iter()
            .filter(|disruption| disruption.reportable())
            .count()
    }

    pub fn has_long_disruption(&self, len_min: usize) -> bool {
        self.disrs.iter().any(|d| d.get_len() >= len_min)
    }

    pub fn contains_hsp(&self, other: &Hsp) -> bool {
        self.q_prot == other.q_prot
            && self.s_prot == other.s_prot
            && self.a_prot == other.a_prot
            && self.qseqid == other.qseqid
            && self.sseqid == other.sseqid
            && self.q_int.contains(&other.q_int)
            && self.s_int.contains(&other.s_int)
    }

    pub fn can_merge_blastx_with(&self, next: &Hsp) -> bool {
        Exon::new(0, 0, self.length).arcable(&Exon::new(1, 0, next.length), &[self, next], true)
    }

    pub fn blastx_arcable(prev: &Hsp, next: &Hsp, bacteria: bool) -> bool {
        const INTRON_MAX: usize = 5000;

        if !bacteria
            || !prev.blastx()
            || !next.blastx()
            || prev.qseqid != next.qseqid
            || prev.sseqid != next.sseqid
            || prev.qlen != next.qlen
            || prev.slen != next.slen
            || prev.s_int.strand != next.s_int.strand
        {
            return false;
        }
        let q_center = (prev.q_int.start + prev.q_int.stop) / 2;
        let next_q_center = (next.q_int.start + next.q_int.stop) / 2;
        if q_center >= next_q_center {
            return false;
        }
        if prev.q_int.stop >= next.q_int.start.saturating_add(50) {
            return false;
        }

        let s_center = (prev.s_int.start + prev.s_int.stop) / 2;
        let next_s_center = (next.s_int.start + next.s_int.stop) / 2;
        if prev.s_int.strand == 1 {
            s_center < next_s_center
                && prev.s_int.stop.saturating_add(INTRON_MAX) >= next.s_int.start
        } else {
            next_s_center < s_center
                && next.s_int.start.saturating_add(INTRON_MAX) >= prev.s_int.stop
        }
    }

    pub fn merge_blastx_chain(hsps: &[Hsp]) -> Option<Hsp> {
        Self::merge_blastx_chain_with_sm(hsps, None)
    }

    fn merge_blastx_chain_with_sm(hsps: &[Hsp], sm: Option<&SubstMat>) -> Option<Hsp> {
        let input_disrs: Vec<Vec<Disruption>> = hsps.iter().map(|hsp| hsp.disrs.clone()).collect();
        let clean_hsps;
        let hsps = if input_disrs.iter().any(|disrs| !disrs.is_empty()) {
            clean_hsps = hsps
                .iter()
                .cloned()
                .map(|mut hsp| {
                    hsp.disrs.clear();
                    hsp
                })
                .collect::<Vec<_>>();
            clean_hsps.as_slice()
        } else {
            hsps
        };

        let mut first = hsps.first()?.clone();
        if hsps.len() == 1 {
            for disr in &input_disrs[0] {
                if let Some(mut clipped) =
                    Self::clip_same_hsp_disruption(&first, disr, 0, first.length)
                {
                    clipped.prev_hsp_idx = Some(0);
                    clipped.next_hsp_idx = Some(0);
                    first.disrs.push(clipped);
                }
            }
            return Some(first);
        }
        if hsps
            .windows(2)
            .any(|pair| !pair[0].can_merge_blastx_with(&pair[1]))
        {
            return None;
        }

        let mut starts = vec![0; hsps.len()];
        let mut stops: Vec<usize> = hsps.iter().map(|hsp| hsp.length).collect();

        for (idx, pair) in hsps.windows(2).enumerate() {
            let prev = &pair[0];
            let next = &pair[1];

            if prev.q_int.stop <= next.q_int.start {
                stops[idx] = prev.length;
                starts[idx + 1] = 0;
                continue;
            }

            let q_start = next.q_int.start;
            let q_stop = prev.q_int.stop;
            let q_len = q_stop - q_start;

            let mut prev_scores = Vec::with_capacity(q_len + 1);
            for _ in q_start..prev.q_int.start {
                prev_scores.push(0.0);
            }
            let mut q_pos = prev.q_int.start;
            for i in 0..prev.length {
                if prev.qseq.as_bytes()[i] != b'-' {
                    if q_start <= q_pos && q_pos < q_stop {
                        prev_scores.push(
                            SubstMat::char2score_static(
                                sm,
                                prev.qseq.as_bytes()[i] as char,
                                prev.sseq.as_bytes()[i] as char,
                            )
                            .expect("amino acid score"),
                        );
                    }
                    q_pos += prev.a2q;
                }
            }
            prev_scores.push(0.0);

            let mut next_scores = Vec::with_capacity(q_len + 1);
            next_scores.push(0.0);
            let mut q_pos = next.q_int.start;
            for i in 0..next.length {
                if next.qseq.as_bytes()[i] != b'-' {
                    if q_start <= q_pos && q_pos < q_stop {
                        next_scores.push(
                            SubstMat::char2score_static(
                                sm,
                                next.qseq.as_bytes()[i] as char,
                                next.sseq.as_bytes()[i] as char,
                            )
                            .expect("amino acid score"),
                        );
                    }
                    q_pos += next.a2q;
                }
            }
            for _ in next.q_int.stop..q_stop {
                next_scores.push(0.0);
            }

            if prev_scores.len() != q_len + 1 || next_scores.len() != q_len + 1 {
                return None;
            }

            for i in (0..q_len).rev() {
                prev_scores[i] += prev_scores[i + 1];
            }
            for i in 1..=q_len {
                next_scores[i] += next_scores[i - 1];
            }

            let mut best_score = f64::INFINITY;
            let mut best_split = None;
            for i in (0..=q_len).rev() {
                let score = prev_scores[i] + next_scores[i];
                if score < best_score {
                    best_score = score;
                    best_split = Some(i);
                }
            }
            let q_split = q_start + best_split?;
            let prev_q_center = (prev.q_int.start + prev.q_int.stop) / 2;
            let next_q_center = (next.q_int.start + next.q_int.stop) / 2;
            if q_split < prev_q_center || q_split >= next_q_center || q_split == 0 {
                return None;
            }

            let prev_start = prev.q2pos(q_split - 1, false) + 1;
            let next_stop = next.q2pos(q_split, true);
            if prev_start == 0 || prev_start <= starts[idx] || next_stop >= stops[idx + 1] {
                return None;
            }
            stops[idx] = prev_start;
            starts[idx + 1] = next_stop;
        }

        let last = hsps.last().expect("non-empty HSP chain");
        let mut merged = first.clone();
        merged.merged = true;
        merged.q_int.start = first.pos2q(starts[0], true);
        merged.q_int.stop = last.pos2q(stops[hsps.len() - 1], true);
        merged.s_int = if first.s_int.strand == 1 {
            Interval::new(
                first.pos2s(starts[0], true),
                last.pos2s(stops[hsps.len() - 1], true),
                first.s_int.strand,
            )
        } else {
            Interval::new(
                last.pos2s(stops[hsps.len() - 1], true),
                first.pos2s(starts[0], true),
                first.s_int.strand,
            )
        };
        merged.qseq = hsps
            .iter()
            .zip(starts.iter().zip(stops.iter()))
            .map(|(hsp, (start, stop))| &hsp.qseq[*start..*stop])
            .collect();
        merged.sseq = hsps
            .iter()
            .zip(starts.iter().zip(stops.iter()))
            .map(|(hsp, (start, stop))| &hsp.sseq[*start..*stop])
            .collect();
        merged.disrs.clear();

        for (idx, hsp) in hsps.iter().enumerate() {
            for disr in &input_disrs[idx] {
                if let Some(mut clipped) =
                    Self::clip_same_hsp_disruption(hsp, disr, starts[idx], stops[idx])
                {
                    clipped.prev_hsp_idx = Some(idx);
                    clipped.next_hsp_idx = Some(idx);
                    merged.disrs.push(clipped);
                }
            }
        }

        for idx in 0..hsps.len() - 1 {
            let prev = &hsps[idx];
            let next = &hsps[idx + 1];
            let prev_start = stops[idx];
            let next_stop = starts[idx + 1];
            let q_interval =
                Interval::new(prev.pos2q(prev_start, true), next.pos2q(next_stop, true), 0);
            let s_start = prev.pos2s(prev_start, true);
            let s_stop = next.pos2s(next_stop, true);
            let s_interval = if prev.s_int.strand == -1 {
                Interval::new(s_stop, s_start, prev.s_int.strand)
            } else {
                Interval::new(s_start, s_stop, prev.s_int.strand)
            };
            let disruption = Disruption {
                prev_hsp_idx: Some(idx),
                next_hsp_idx: Some(idx + 1),
                prev_start,
                next_stop,
                prev_slen: prev.slen,
                intron: true,
                s_stop_codon: false,
                q_interval,
                s_interval,
            };
            if !disruption.q_interval.valid() || !disruption.s_interval.valid() {
                return None;
            }
            if disruption.reportable() {
                merged.disrs.push(disruption);
            }
        }

        merged.finish_hsp(false, false);
        Some(merged)
    }

    fn clip_same_hsp_disruption(
        hsp: &Hsp,
        disr: &Disruption,
        start: usize,
        stop: usize,
    ) -> Option<Disruption> {
        if disr.prev_hsp_idx != disr.next_hsp_idx || disr.intron {
            return None;
        }
        assert!(start <= stop);

        let prev_start = disr.prev_start.clamp(start, stop);
        let next_stop = disr.next_stop.clamp(start, stop);
        assert!(prev_start <= next_stop);

        let q_interval = Interval::new(
            hsp.pos2q_[hsp.pos2real_q(prev_start, true)],
            hsp.pos2q_[hsp.pos2real_q(next_stop, true)],
            0,
        );
        let s_start = hsp.pos2s_[hsp.pos2real_s(prev_start, true)];
        let s_stop = hsp.pos2s_[hsp.pos2real_s(next_stop, true)];
        let s_interval = if hsp.s_int.strand == -1 {
            Interval::new(s_stop, s_start, hsp.s_int.strand)
        } else {
            Interval::new(s_start, s_stop, hsp.s_int.strand)
        };
        let clipped = Disruption {
            prev_hsp_idx: disr.prev_hsp_idx,
            next_hsp_idx: disr.next_hsp_idx,
            prev_start,
            next_stop,
            prev_slen: hsp.slen,
            intron: false,
            s_stop_codon: hsp.sseq[prev_start..next_stop].contains('*'),
            q_interval,
            s_interval,
        };

        if clipped.reportable() {
            Some(clipped)
        } else {
            None
        }
    }

    /// Sort comparator: sseqid, strand, qseqid, sInt.start, sInt.stop
    pub fn less(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
        a.sseqid
            .cmp(&b.sseqid)
            .then(a.s_int.strand.cmp(&b.s_int.strand))
            .then(a.qseqid.cmp(&b.qseqid))
            .then(a.s_int.start.cmp(&b.s_int.start))
            .then(a.s_int.stop.cmp(&b.s_int.stop))
    }

    pub fn clean_qseq(&self) -> String {
        self.qseq.replace('-', "")
    }

    pub fn clean_sseq(&self) -> String {
        self.sseq.replace('-', "")
    }

    pub fn q_len_real(&self) -> usize {
        self.length - self.qgap
    }

    pub fn s_len_real(&self) -> usize {
        self.length - self.sgap
    }

    pub fn qc(&self) {
        assert!(self.merged || self.disrs.is_empty() || self.blastx());
        assert!(!(self.q_prot || self.s_prot) || self.a_prot);
        assert!(!self.s_prot || self.q_prot);

        if self.empty() {
            assert!(self.qseqid.is_empty());
            assert!(self.qseq.is_empty());
            assert!(self.q_int.empty());
            assert_eq!(self.qlen, NO_INDEX);
            assert!(self.sseqid.is_empty());
            assert!(self.sseq.is_empty());
            assert!(self.s_int.empty());
            assert_eq!(self.slen, NO_INDEX);
            assert_eq!(self.nident, NO_INDEX);
            assert_eq!(self.qgap, NO_INDEX);
            assert_eq!(self.sgap, NO_INDEX);
            assert_eq!(self.qx, NO_INDEX);
            assert_eq!(self.sx, NO_INDEX);
            return;
        }

        assert!(!self.qseqid.is_empty());
        assert!(!self.sseqid.is_empty());
        assert!(self.q_int.valid());
        assert!(self.s_int.valid());
        assert!(self.s_int.stop <= self.slen);
        assert!(self.q_int.stop <= self.qlen);
        assert_eq!(self.qseq.len(), self.sseq.len());
        assert!(!self.qseq.is_empty());
        assert_eq!(self.s_int.strand != 0, !self.s_prot);
        if self.disrs.is_empty() {
            assert_eq!(self.q_int.len() % self.a2q, 0);
            assert_eq!(self.s_int.len() % self.a2s, 0);
        }

        for (q, s) in self.qseq.bytes().zip(self.sseq.bytes()) {
            assert!(q != b'-' || s != b'-');
        }
        if self.s_truncated() {
            assert!(!self.q_complete());
        }

        assert!(self.length > 0);
        if !self.merged {
            assert!(self.nident > 0);
        }
        assert!(self.qx <= self.length);
        assert!(self.sx <= self.length);
        assert!(self.nident + self.qgap <= self.length);
        assert!(self.nident + self.sgap <= self.length);
        assert_eq!(self.length, self.qseq.len());
        if self.disrs.is_empty() {
            if !self.merged {
                assert_ne!(self.qseq.as_bytes().first(), Some(&b'-'));
                assert_ne!(self.sseq.as_bytes().first(), Some(&b'-'));
                assert_ne!(self.qseq.as_bytes().last(), Some(&b'-'));
                assert_ne!(self.sseq.as_bytes().last(), Some(&b'-'));
            }
            assert!(self.nident <= self.q_len());
            assert!(self.nident <= self.s_len());
            assert!(self.q_len() <= self.length);
            assert!(self.s_len() <= self.length);
            if self.a_prot && !self.q_prot {
                assert_eq!(self.q_abs_coverage() % 3, 0);
            }
            if self.a_prot && !self.s_prot {
                assert_eq!(self.s_abs_coverage() % 3, 0);
            }
        } else {
            assert!(self.merged);
            for disr in &self.disrs {
                disr.qc();
                assert!(!matches!(
                    disr.disruption_type(),
                    DisruptionType::None | DisruptionType::Smooth
                ));
                assert_eq!(disr.s_interval.strand, self.s_int.strand);
                assert!(self.blastx());
            }
        }

        assert_eq!(self.sframe != 0, self.a_prot && !self.s_prot);
        if self.a_prot && !self.s_prot {
            assert_eq!(self.s_int.strand == -1, self.sframe < 0);
            assert!(is_frame(self.sframe));
        }
        assert!(self.c_complete.is_none() || self.a_prot);
        assert!(!self.s_internal_stop || self.a_prot);
    }

    pub fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        write!(
            os,
            "Hsp: merged={} {}({}) {} {}({}) {}",
            self.merged as i32,
            self.qseqid,
            self.qlen,
            self.q_int.format_display(),
            self.sseqid,
            self.slen,
            self.s_int.format_display()
        )?;
        if self.disrs.is_empty() {
            write!(os, " qLen={} sLen={}", self.q_len(), self.s_len())?;
        }
        write!(
            os,
            " length={} nident={} qgap={} sgap={} qx={} sx={} sframe={} sInternalStop={} #disrs={}",
            self.length,
            self.nident,
            self.qgap,
            self.sgap,
            self.qx,
            self.sx,
            self.sframe,
            self.s_internal_stop as i32,
            self.disrs.len()
        )?;
        for disr in &self.disrs {
            write!(os, " ")?;
            disr.save_text(os)?;
        }
        Ok(())
    }

    /// Compute alignment statistics from qseq/sseq.
    /// Must be called after construction.
    pub fn finish_hsp(&mut self, q_stop_codon: bool, bacterial_start_codon: bool) {
        self.a2q = if self.a_prot && !self.q_prot { 3 } else { 1 };
        self.a2s = if self.a_prot && !self.s_prot { 3 } else { 1 };

        if !self.merged {
            self.move_dashes_right();
            self.c_complete = None;

            // Handle stop codon at query end
            if q_stop_codon && self.a_prot && !self.qseq.is_empty() {
                let q_stop = self.q_int.stop;
                let q_len = self.qlen;
                if q_stop == q_len {
                    assert_eq!(
                        self.qseq.as_bytes().last().copied(),
                        Some(b'*'),
                        "ending stop codon is expected"
                    );
                    self.c_complete = Some(self.sseq.as_bytes().last().copied() == Some(b'*'));
                    self.erase_qseq_back();
                    self.erase_sseq_back();
                } else if q_stop == q_len.saturating_sub(1) && self.s_tail(false) >= self.a2s {
                    self.c_complete = Some(false);
                }
                self.qlen = self.qlen.saturating_sub(1);
            }

            if bacterial_start_codon
                && self.q_prot
                && self.q_int.start == 0
                && self.qseq.starts_with('M')
                && matches!(self.sseq.as_bytes().first(), Some(b'L' | b'I' | b'V'))
            {
                self.sseq.replace_range(0..1, "M");
            }
        }

        self.length = self.qseq.len();
        assert_eq!(self.length, self.sseq.len());

        self.nident = 0;
        self.qgap = 0;
        self.sgap = 0;
        self.qx = 0;
        self.sx = 0;

        let q_bytes = self.qseq.as_bytes();
        let s_bytes = self.sseq.as_bytes();

        for i in 0..self.length {
            let qc = q_bytes[i] as char;
            let sc = s_bytes[i] as char;

            if qc == '-' {
                self.qgap += 1;
            } else if sc == '-' {
                self.sgap += 1;
            } else if self.char_match(i) {
                self.nident += 1;
            }

            if is_ambig(qc, self.a_prot) {
                self.qx += 1;
            }
            if is_ambig(sc, self.a_prot) {
                self.sx += 1;
            }
        }
        self.sframe = 0;
        if self.a_prot && !self.s_prot {
            self.sframe = self.s_int.frame();
        }
        self.s_internal_stop = self.a_prot && self.sseq.contains('*');

        self.pos2q_.clear();
        self.pos2s_.clear();

        if !self.disrs.is_empty() {
            return;
        }

        self.pos2q_.reserve(self.length + 1);
        let mut q_pos = self.q_int.start;
        for i in 0..=self.length {
            assert!(q_pos >= self.q_int.start);
            assert!(q_pos <= self.q_int.stop);
            self.pos2q_.push(q_pos);
            if i < self.length && q_bytes[i] != b'-' {
                q_pos += self.a2q;
            }
        }

        self.pos2s_.reserve(self.length + 1);
        let mut s_pos = if self.s_int.strand == -1 {
            self.s_int.stop
        } else {
            self.s_int.start
        };
        for i in 0..=self.length {
            assert!(s_pos >= self.s_int.start);
            assert!(s_pos <= self.s_int.stop);
            self.pos2s_.push(s_pos);
            if i < self.length && s_bytes[i] != b'-' {
                if self.s_int.strand == -1 {
                    s_pos = s_pos.saturating_sub(self.a2s);
                } else {
                    s_pos += self.a2s;
                }
            }
        }
    }

    pub fn pos2real_q(&self, mut pos: usize, forward: bool) -> usize {
        while pos < self.length && self.qseq.as_bytes()[pos] == b'-' {
            if forward {
                pos += 1;
            } else {
                assert!(pos > 0);
                pos -= 1;
            }
        }
        pos
    }

    pub fn pos2real_s(&self, mut pos: usize, forward: bool) -> usize {
        while pos < self.length && self.sseq.as_bytes()[pos] == b'-' {
            if forward {
                pos += 1;
            } else {
                assert!(pos > 0);
                pos -= 1;
            }
        }
        pos
    }

    pub fn pos2q(&self, pos: usize, forward: bool) -> usize {
        assert!(self.disrs.is_empty());
        self.pos2q_[self.pos2real_q(pos, forward)]
    }

    pub fn pos2s(&self, pos: usize, forward: bool) -> usize {
        assert!(self.disrs.is_empty());
        self.pos2s_[self.pos2real_s(pos, forward)]
    }

    pub fn q2pos(&self, q_pos: usize, forward: bool) -> usize {
        assert!(self.disrs.is_empty());
        let pos = self.pos2q_.partition_point(|pos| *pos < q_pos);
        assert!(pos < self.pos2q_.len());
        self.pos2real_q(pos, forward)
    }

    pub fn s2pos(&self, s_pos: usize, forward: bool) -> usize {
        assert!(self.disrs.is_empty());
        let pos = if self.s_int.strand == -1 {
            self.pos2s_.partition_point(|pos| *pos > s_pos)
        } else {
            self.pos2s_.partition_point(|pos| *pos < s_pos)
        };
        assert!(pos < self.pos2s_.len());
        self.pos2real_s(pos, forward)
    }

    pub fn q_map(&self, len: usize) -> String {
        assert!(self.qlen <= len);

        let mut map = String::with_capacity(len);
        map.extend(std::iter::repeat_n('-', self.q_int.start));
        for (q, s) in self.qseq.bytes().zip(self.sseq.bytes()) {
            if q != b'-' {
                map.push(s as char);
            }
        }
        map.extend(std::iter::repeat_n('-', len - self.q_int.stop));
        map
    }

    pub fn q_better_eq(&self, other: &Hsp) -> bool {
        assert_eq!(self.a_prot, other.a_prot);
        assert_eq!(self.sseqid, other.sseqid);
        assert_eq!(self.s_int.strand, other.s_int.strand);

        if !other.s_inside_eq(self, 0) {
            return false;
        }
        if other.rel_identity() < self.rel_identity() {
            return true;
        }
        if self.rel_identity() < other.rel_identity() {
            return false;
        }
        if other.nident < self.nident {
            return true;
        }
        if self.nident < other.nident {
            return false;
        }
        if self.qseqid < other.qseqid {
            return true;
        }
        if other.qseqid < self.qseqid {
            return false;
        }
        true
    }

    pub fn s_start_global(&self) -> isize {
        if self.s_int.strand == 1 {
            self.s_int.start as isize - 3 * self.q_int.start as isize
        } else {
            self.s_int.stop as isize + 3 * self.q_int.start as isize
        }
    }

    fn move_dashes_right(&mut self) {
        loop {
            let b1 = Self::move_dashes_right_one(&self.qseq, &mut self.sseq);
            let b2 = Self::move_dashes_right_one(&self.sseq, &mut self.qseq);
            if !b1 && !b2 {
                break;
            }
        }

        while !self.qseq.is_empty()
            && (self.qseq.as_bytes()[0] == b'-' || self.sseq.as_bytes()[0] == b'-')
        {
            self.erase_qseq_front();
            self.erase_sseq_front();
        }

        while !self.qseq.is_empty()
            && (self.qseq.as_bytes().last() == Some(&b'-')
                || self.sseq.as_bytes().last() == Some(&b'-'))
        {
            self.erase_qseq_back();
            self.erase_sseq_back();
        }
    }

    fn move_dashes_right_one(seq1: &str, seq2: &mut String) -> bool {
        assert_eq!(seq1.len(), seq2.len());
        let seq1 = seq1.as_bytes();
        let mut seq2_bytes = seq2.as_bytes().to_vec();
        let mut changed = false;
        let mut start = None;

        for i in 0..seq1.len() {
            if let Some(start_idx) = start {
                if seq2_bytes[i] != b'-' {
                    if seq1[i] == seq2_bytes[i] && seq2_bytes[i] == seq1[start_idx] {
                        seq2_bytes.swap(i, start_idx);
                        changed = true;
                    } else {
                        start = None;
                    }
                }
            } else if seq2_bytes[i] == b'-' {
                start = Some(i);
            }
        }

        if changed {
            *seq2 = String::from_utf8(seq2_bytes).expect("BLAST sequence is valid UTF-8");
        }
        changed
    }

    fn erase_qseq_front(&mut self) {
        if self.qseq.as_bytes().first() != Some(&b'-') {
            self.q_int.start += self.a2q;
        }
        self.qseq.remove(0);
    }

    fn erase_sseq_front(&mut self) {
        if self.sseq.as_bytes().first() != Some(&b'-') {
            if self.s_int.strand == -1 {
                self.s_int.stop = self.s_int.stop.saturating_sub(self.a2s);
            } else {
                self.s_int.start += self.a2s;
            }
        }
        self.sseq.remove(0);
    }

    fn erase_qseq_back(&mut self) {
        if self.qseq.as_bytes().last() != Some(&b'-') {
            self.q_int.stop = self.q_int.stop.saturating_sub(self.a2q);
        }
        self.qseq.pop();
    }

    fn erase_sseq_back(&mut self) {
        if self.sseq.as_bytes().last() != Some(&b'-') {
            if self.s_int.strand == -1 {
                self.s_int.start += self.a2s;
            } else {
                self.s_int.stop = self.s_int.stop.saturating_sub(self.a2s);
            }
        }
        self.sseq.pop();
    }

    pub fn s_tail(&self, left: bool) -> usize {
        if self.s_int.empty() {
            return 0;
        }
        if (left && self.s_int.strand == 1) || (!left && self.s_int.strand == -1) {
            self.s_int.start
        } else {
            self.slen.saturating_sub(self.s_int.stop)
        }
    }

    /// Parse from BLAST tabular output line
    pub fn from_blast_line(
        line: &str,
        q_prot: bool,
        s_prot: bool,
        a_prot: bool,
        q_stop_codon: bool,
        bacterial_start_codon: bool,
    ) -> anyhow::Result<Self> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 10 {
            anyhow::bail!("BLAST line needs at least 10 fields");
        }

        // Format: qseqid sseqid qstart qend qlen sstart send slen qseq sseq
        let qseqid = fields[0].to_string();
        let sseqid = fields[1].to_string();
        let mut qstart: usize = fields[2].parse()?;
        let mut qend: usize = fields[3].parse()?;
        let qlen: usize = fields[4].parse()?;
        let mut sstart: usize = fields[5].parse()?;
        let mut send: usize = fields[6].parse()?;
        let slen: usize = fields[7].parse()?;
        let mut qseq = if q_prot || a_prot {
            fields[8].to_uppercase()
        } else {
            fields[8].to_lowercase()
        };
        let mut sseq = if s_prot || a_prot {
            fields[9].to_uppercase()
        } else {
            fields[9].to_lowercase()
        };

        let mut s_strand: Strand = 0;
        if !s_prot {
            s_strand = 1;
            if q_prot {
                if sstart > send {
                    s_strand = -1;
                    std::mem::swap(&mut sstart, &mut send);
                }
            } else if qstart > qend {
                std::mem::swap(&mut qstart, &mut qend);
                if !a_prot {
                    reverse_dna(&mut qseq);
                    reverse_dna(&mut sseq);
                }
                s_strand = -1;
            }
        }

        // Convert from 1-based to 0-based
        let qstart = qstart.saturating_sub(1);
        let sstart = sstart.saturating_sub(1);

        let a2q = if !q_prot && a_prot { 3 } else { 1 };
        let a2s = if !s_prot && a_prot { 3 } else { 1 };

        let mut hsp = Hsp {
            merged: false,
            q_prot,
            s_prot,
            a_prot,
            a2q,
            a2s,
            qseqid,
            sseqid,
            q_int: Interval::new(qstart, qend, 0),
            s_int: Interval::new(sstart, send, s_strand),
            qlen,
            slen,
            qseq,
            sseq,
            length: 0,
            nident: 0,
            qgap: 0,
            sgap: 0,
            qx: 0,
            sx: 0,
            pos2q_: Vec::new(),
            pos2s_: Vec::new(),
            sframe: 0,
            c_complete: None,
            s_internal_stop: false,
            disrs: Vec::new(),
        };

        hsp.finish_hsp(q_stop_codon, bacterial_start_codon);
        Ok(hsp)
    }
}

pub struct HspMerge<'a> {
    orig_hsps: Vec<Hsp>,
    sm: Option<&'a SubstMat>,
    intron_score: AlignScore,
    bacteria: bool,
    used_exons: BTreeSet<(usize, usize, usize)>,
}

impl<'a> HspMerge<'a> {
    pub fn new(
        orig_hsps: Vec<Hsp>,
        sm: Option<&'a SubstMat>,
        intron_score: AlignScore,
        bacteria: bool,
    ) -> Self {
        let mut seen = BTreeSet::new();
        for hsp in &orig_hsps {
            assert!(!hsp.merged);
            let key = (
                hsp.qseqid.clone(),
                hsp.sseqid.clone(),
                hsp.q_int.start,
                hsp.q_int.stop,
                hsp.s_int.start,
                hsp.s_int.stop,
                hsp.qseq.clone(),
                hsp.sseq.clone(),
            );
            assert!(seen.insert(key), "Duplicate HSP");
        }

        HspMerge {
            orig_hsps,
            sm,
            intron_score,
            bacteria,
            used_exons: BTreeSet::new(),
        }
    }

    pub fn get(&mut self) -> (Hsp, Option<usize>, AlignScore) {
        let mut best_score = -f64::INFINITY;
        let mut best_start = None;
        let mut best_exons = Vec::new();
        let mut best_hsp = Hsp::default();
        let mut exons = Vec::new();
        for (idx, hsp) in self.orig_hsps.iter().enumerate() {
            let mut hsp_exons: Vec<Exon> = hsp2exons(idx, hsp, self.sm);
            hsp_exons.retain(|exon| {
                !self
                    .used_exons
                    .contains(&(exon.hsp_idx, exon.start, exon.stop))
            });
            exons.extend(hsp_exons);
        }
        exons.sort_by_key(|exon| (exon.hsp_idx, exon.start, exon.stop));
        let hsp_refs: Vec<&Hsp> = self.orig_hsps.iter().collect();
        let mut introns = Vec::new();
        for prev in &exons {
            for next in &exons {
                if !prev.arcable(next, &hsp_refs, self.bacteria) {
                    continue;
                }
                introns.push(Intron::new(
                    prev.clone(),
                    next.clone(),
                    &self.orig_hsps,
                    self.sm,
                ));
            }
        }

        for exon in &exons {
            let mut memo = vec![None; exons.len()];
            let (score, chain_exons) = exon.merge_tail(
                &exons,
                &introns,
                &self.orig_hsps,
                &mut memo,
                self.bacteria,
                self.intron_score,
                self.sm,
            );
            let mut chain = Vec::with_capacity(chain_exons.len());
            for chain_exon in &chain_exons {
                let orig_hsp = &self.orig_hsps[chain_exon.hsp_idx];
                let q_int = chain_exon.q_interval(&self.orig_hsps);
                let s_int = chain_exon.s_interval(&self.orig_hsps);
                let mut hsp = orig_hsp.clone();
                hsp.q_int = q_int;
                hsp.s_int = s_int;
                hsp.qseq = orig_hsp.qseq[chain_exon.start..chain_exon.stop].to_string();
                hsp.sseq = orig_hsp.sseq[chain_exon.start..chain_exon.stop].to_string();
                hsp.disrs.clear();
                hsp.merged = false;
                hsp.finish_hsp(false, false);
                hsp.disrs = chain_exon
                    .disrs
                    .iter()
                    .filter_map(|disr| {
                        let mut local_disr = disr.clone();
                        local_disr.prev_start = local_disr
                            .prev_start
                            .clamp(chain_exon.start, chain_exon.stop)
                            - chain_exon.start;
                        local_disr.next_stop = local_disr
                            .next_stop
                            .clamp(chain_exon.start, chain_exon.stop)
                            - chain_exon.start;
                        local_disr.prev_hsp_idx = Some(0);
                        local_disr.next_hsp_idx = Some(0);
                        local_disr.prev_slen = hsp.slen;
                        Hsp::clip_same_hsp_disruption(&hsp, &local_disr, 0, hsp.length)
                    })
                    .collect();
                chain.push(hsp);
            }

            let Some(mut hsp_new) = Hsp::merge_blastx_chain_with_sm(&chain, self.sm) else {
                continue;
            };
            if !hsp_new.merged {
                hsp_new.merged = true;
                hsp_new.finish_hsp(false, false);
            }
            if hsp_new.nident == 0 {
                continue;
            }

            if score > best_score {
                best_score = score;
                best_start = Some(exon.hsp_idx);
                best_exons = chain_exons;
                best_hsp = hsp_new;
            }
        }

        let Some(start) = best_start else {
            return (Hsp::default(), None, -f64::INFINITY);
        };
        for exon in best_exons {
            self.used_exons
                .insert((exon.hsp_idx, exon.start, exon.stop));
        }

        best_hsp.disrs.sort_by_key(|disr| {
            (
                disr.q_interval.start,
                disr.q_interval.stop,
                disr.s_interval.start,
                disr.s_interval.stop,
                disr.intron,
            )
        });
        best_hsp.qc();
        (best_hsp, Some(start), best_score)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Exon {
    pub is_insertion: bool,
    pub hsp_idx: usize,
    pub start: usize,
    pub stop: usize,
    pub disrs: Vec<Disruption>,
}

impl Exon {
    pub fn new(hsp_idx: usize, start: usize, stop: usize) -> Self {
        Self {
            is_insertion: false,
            hsp_idx,
            start,
            stop,
            disrs: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.stop - self.start
    }

    pub fn q_interval(&self, hsps: &[Hsp]) -> Interval {
        let hsp = &hsps[self.hsp_idx];
        Interval::new(hsp.pos2q(self.start, true), hsp.pos2q(self.stop, true), 0)
    }

    pub fn s_interval(&self, hsps: &[Hsp]) -> Interval {
        let hsp = &hsps[self.hsp_idx];
        let start = hsp.pos2s(self.start, true);
        let stop = hsp.pos2s(self.stop, true);
        if hsp.s_int.strand == -1 {
            Interval::new(stop, start, hsp.s_int.strand)
        } else {
            Interval::new(start, stop, hsp.s_int.strand)
        }
    }

    pub fn arcable(&self, next: &Exon, hsps: &[&Hsp], bacteria: bool) -> bool {
        let hsp = hsps[self.hsp_idx];
        let next_hsp = hsps[next.hsp_idx];

        if hsp.qseqid != next_hsp.qseqid || hsp.sseqid != next_hsp.sseqid {
            return false;
        }

        if self == next {
            return false;
        }
        if self.hsp_idx == next.hsp_idx {
            if bacteria {
                return self.stop == next.start;
            }
        } else if next.is_insertion {
            return false;
        }

        let intron_max = if bacteria { 5000 } else { 30000 };

        if hsp.s_int.strand != next_hsp.s_int.strand {
            return false;
        }

        let q_start = hsp.pos2q(self.start, true);
        let q_stop = hsp.pos2q(self.stop, true);
        let next_q_start = next_hsp.pos2q(next.start, true);
        let next_q_stop = next_hsp.pos2q(next.stop, true);
        let q_center = (q_start + q_stop) / 2;
        let next_q_center = (next_q_start + next_q_stop) / 2;
        if q_center >= next_q_center {
            return false;
        }

        let s_start = hsp.pos2s(self.start, true);
        let s_stop = hsp.pos2s(self.stop, true);
        let next_s_start = next_hsp.pos2s(next.start, true);
        let next_s_stop = next_hsp.pos2s(next.stop, true);
        let s_center = (s_start + s_stop) / 2;
        let next_s_center = (next_s_start + next_s_stop) / 2;
        if hsp.s_int.strand == 1 {
            if s_center >= next_s_center {
                return false;
            }
        } else if next_s_center >= s_center {
            return false;
        }

        if bacteria && q_stop >= next_q_start.saturating_add(50) {
            return false;
        }

        if hsp.s_int.strand == 1 && s_stop.saturating_add(intron_max) < next_s_start {
            return false;
        }
        if hsp.s_int.strand == -1 && next_s_start.saturating_add(intron_max) < s_stop {
            return false;
        }

        true
    }

    pub fn set_best_intron(
        &self,
        exons: &[Exon],
        introns: &[Intron],
        hsps: &[Hsp],
        memo: &mut [Option<(AlignScore, Vec<Exon>)>],
        bacteria: bool,
        intron_score: AlignScore,
        sm: Option<&SubstMat>,
    ) -> Option<(Intron, AlignScore, Vec<Exon>)> {
        let mut best: Option<(Intron, AlignScore, Vec<Exon>)> = None;

        for intron in introns.iter().filter(|intron| intron.prev == *self) {
            let next = &intron.next;
            if !exons.iter().any(|exon| exon == next) {
                continue;
            }
            let score =
                intron.get_total_score(exons, introns, hsps, memo, bacteria, intron_score, sm);
            if best
                .as_ref()
                .is_none_or(|(_, best_score, _)| score > *best_score)
            {
                let (_, tail_indices) =
                    next.merge_tail(exons, introns, hsps, memo, bacteria, intron_score, sm);
                best = Some((intron.clone(), score, tail_indices));
            }
        }

        best
    }

    pub fn merge_tail(
        &self,
        exons: &[Exon],
        introns: &[Intron],
        hsps: &[Hsp],
        memo: &mut [Option<(AlignScore, Vec<Exon>)>],
        bacteria: bool,
        intron_score: AlignScore,
        sm: Option<&SubstMat>,
    ) -> (AlignScore, Vec<Exon>) {
        let self_pos = exons
            .iter()
            .position(|exon| exon == self)
            .expect("Exon must be present in merge graph");
        if let Some(cached) = &memo[self_pos] {
            return cached.clone();
        }

        let hsp = &hsps[self.hsp_idx];
        let mut score = 0.0;
        for i in self.start..self.stop {
            score += SubstMat::char2score_static(
                sm,
                hsp.qseq.as_bytes()[i] as char,
                hsp.sseq.as_bytes()[i] as char,
            )
            .expect("amino acid score");
        }
        let mut best = (score, vec![self.clone()]);

        if let Some((_, tail_score, tail_indices)) =
            self.set_best_intron(exons, introns, hsps, memo, bacteria, intron_score, sm)
        {
            let total_score = score + tail_score;
            if total_score > best.0 {
                let mut indices = Vec::with_capacity(tail_indices.len() + 1);
                indices.push(self.clone());
                indices.extend(tail_indices);
                best = (total_score, indices);
            }
        }

        memo[self_pos] = Some(best.clone());
        best
    }
}

#[derive(Debug, Clone)]
struct DensityState {
    weight_local: [f64; 2],
    weight_global: [f64; 2],
    prev_global_hi_dens: [bool; 2],
}

impl DensityState {
    const LO_DENS_PROB: f64 = 0.10;
    const HI_DENS_PROB: f64 = 0.90;
    const DENS_CHANGE_PROB: f64 = 0.005;

    fn new(match_: bool) -> Self {
        let lo_dens_weight = if match_ {
            -Self::LO_DENS_PROB.ln()
        } else {
            -(1.0 - Self::LO_DENS_PROB).ln()
        };
        let hi_dens_weight = if match_ {
            -Self::HI_DENS_PROB.ln()
        } else {
            -(1.0 - Self::HI_DENS_PROB).ln()
        };
        Self {
            weight_local: [lo_dens_weight, hi_dens_weight],
            weight_global: [0.0, 0.0],
            prev_global_hi_dens: [false, true],
        }
    }

    #[allow(dead_code)]
    fn save_text(&self, os: &mut dyn Write) -> anyhow::Result<()> {
        write!(
            os,
            "{:.2}\t{:.2}\t{:.2}\t{:.2}\t{}\t{}",
            self.weight_local[0],
            self.weight_local[1],
            self.weight_global[0],
            self.weight_global[1],
            self.prev_global_hi_dens[0] as i32,
            self.prev_global_hi_dens[1] as i32
        )?;
        Ok(())
    }
}

pub fn hsp2exons(hsp_idx: usize, hsp: &Hsp, _sm: Option<&SubstMat>) -> Vec<Exon> {
    assert!(hsp.length > 0);

    let mut dss: Vec<DensityState> = Vec::with_capacity(hsp.length);
    for i in 0..hsp.length {
        let mut ds = DensityState::new(hsp.char_match(i));
        if i > 0 {
            for hi_dens in [false, true] {
                let hi_idx = hi_dens as usize;
                ds.weight_global[hi_idx] = f64::INFINITY;
                for prev_hi_dens in [false, true] {
                    let prev_idx = prev_hi_dens as usize;
                    let dens_change_weight = if hi_dens != prev_hi_dens {
                        -DensityState::DENS_CHANGE_PROB.ln()
                    } else {
                        -(1.0 - DensityState::DENS_CHANGE_PROB).ln()
                    };
                    let weight_new = ds.weight_local[hi_idx]
                        + dss[i - 1].weight_global[prev_idx]
                        + dens_change_weight;
                    if weight_new < ds.weight_global[hi_idx] {
                        ds.weight_global[hi_idx] = weight_new;
                        ds.prev_global_hi_dens[hi_idx] = prev_hi_dens;
                    }
                }
            }
        } else {
            ds.weight_global = ds.weight_local;
        }
        dss.push(ds);
    }

    let mut exons = Vec::new();
    let mut hi_dens = dss[hsp.length - 1].weight_global[true as usize]
        <= dss[hsp.length - 1].weight_global[false as usize];
    let mut stop = hsp.length;
    let mut start = hsp.length - 1;
    while start > 0 {
        while start > 0 && dss[start].prev_global_hi_dens[hi_dens as usize] == hi_dens {
            start -= 1;
        }
        let mut exon = Exon::new(hsp_idx, start, stop);
        exon.is_insertion = !hi_dens;
        for i in start..stop {
            if hsp.sseq.as_bytes()[i] != b'*' {
                continue;
            }
            let mut i_prev = hsp.pos2real_q(i, false);
            let mut i_next = hsp.pos2real_q(i, true);
            if i_prev == i_next {
                i_next += 1;
            } else {
                i_prev += 1;
            }
            let q_interval = Interval::new(hsp.pos2q(i_prev, true), hsp.pos2q(i_next, true), 0);
            let s_start = hsp.pos2s(i_prev, true);
            let s_stop = hsp.pos2s(i_next, true);
            let s_interval = if hsp.s_int.strand == -1 {
                Interval::new(s_stop, s_start, hsp.s_int.strand)
            } else {
                Interval::new(s_start, s_stop, hsp.s_int.strand)
            };
            let disruption = Disruption {
                prev_hsp_idx: Some(hsp_idx),
                next_hsp_idx: Some(hsp_idx),
                prev_start: i_prev,
                next_stop: i_next,
                prev_slen: hsp.slen,
                intron: false,
                s_stop_codon: hsp.sseq[i_prev..i_next].contains('*'),
                q_interval,
                s_interval,
            };
            if disruption.reportable() {
                exon.disrs.push(disruption);
            }
        }
        exons.push(exon);
        stop = start;
        hi_dens = !hi_dens;
    }

    exons
}

#[derive(Debug, Clone, PartialEq)]
pub struct Intron {
    pub prev: Exon,
    pub next: Exon,
    pub score: AlignScore,
    pub prev_start: usize,
    pub next_stop: usize,
    pub disr: Disruption,
}

impl Intron {
    pub fn new(prev: Exon, next: Exon, hsps: &[Hsp], sm: Option<&SubstMat>) -> Self {
        let mut intron = Self {
            prev,
            next,
            score: f64::MAX,
            prev_start: NO_INDEX,
            next_stop: NO_INDEX,
            disr: Disruption {
                prev_hsp_idx: None,
                next_hsp_idx: None,
                prev_start: 0,
                next_stop: 0,
                prev_slen: NO_INDEX,
                intron: false,
                s_stop_codon: false,
                q_interval: Interval::default(),
                s_interval: Interval::default(),
            },
        };

        let prev_hsp = &hsps[intron.prev.hsp_idx];
        let next_hsp = &hsps[intron.next.hsp_idx];
        let q_start = next_hsp.pos2q(intron.next.start, true);
        let q_stop = prev_hsp.pos2q(intron.prev.stop, true);
        if q_stop <= q_start {
            intron.score = 0.0;
            intron.prev_start = intron.prev.stop;
            intron.next_stop = intron.next.start;
        } else {
            let q_len = q_stop - q_start;

            let mut prev_scores = Vec::with_capacity(q_len + 1);
            for _ in q_start..prev_hsp.pos2q(intron.prev.start, true) {
                prev_scores.push(0.0);
            }
            let mut q_pos = prev_hsp.pos2q(intron.prev.start, true);
            for i in intron.prev.start..intron.prev.stop {
                if prev_hsp.qseq.as_bytes()[i] != b'-' {
                    if q_start <= q_pos && q_pos < q_stop {
                        prev_scores.push(
                            SubstMat::char2score_static(
                                sm,
                                prev_hsp.qseq.as_bytes()[i] as char,
                                prev_hsp.sseq.as_bytes()[i] as char,
                            )
                            .expect("amino acid score"),
                        );
                    }
                    q_pos += prev_hsp.a2q;
                }
            }
            prev_scores.push(0.0);
            assert_eq!(prev_scores.len(), q_len + 1);

            let mut next_scores = Vec::with_capacity(q_len + 1);
            next_scores.push(0.0);
            let mut q_pos = next_hsp.pos2q(intron.next.start, true);
            for i in intron.next.start..intron.next.stop {
                if next_hsp.qseq.as_bytes()[i] != b'-' {
                    if q_start <= q_pos && q_pos < q_stop {
                        next_scores.push(
                            SubstMat::char2score_static(
                                sm,
                                next_hsp.qseq.as_bytes()[i] as char,
                                next_hsp.sseq.as_bytes()[i] as char,
                            )
                            .expect("amino acid score"),
                        );
                    }
                    q_pos += next_hsp.a2q;
                }
            }
            for _ in next_hsp.pos2q(intron.next.stop, true)..q_stop {
                next_scores.push(0.0);
            }
            assert_eq!(next_scores.len(), q_len + 1);

            for i in (0..q_len).rev() {
                prev_scores[i] += prev_scores[i + 1];
            }
            for i in 1..=q_len {
                next_scores[i] += next_scores[i - 1];
            }

            let mut best_split = NO_INDEX;
            for i in (0..=q_len).rev() {
                let score = prev_scores[i] + next_scores[i];
                if score < intron.score {
                    intron.score = score;
                    best_split = i;
                }
            }
            assert_ne!(best_split, NO_INDEX);

            let q_split = q_start + best_split;
            let prev_q_center = (prev_hsp.pos2q(intron.prev.start, true)
                + prev_hsp.pos2q(intron.prev.stop, true))
                / 2;
            let next_q_center = (next_hsp.pos2q(intron.next.start, true)
                + next_hsp.pos2q(intron.next.stop, true))
                / 2;
            if q_split < prev_q_center || q_split >= next_q_center || q_split == 0 {
                intron.score = f64::MAX;
                return intron;
            }

            intron.prev_start = prev_hsp.q2pos(q_split - 1, false) + 1;
            intron.next_stop = next_hsp.q2pos(q_split, true);
        }

        if intron.prev_start <= intron.prev.start || intron.next_stop >= intron.next.stop {
            intron.score = f64::MAX;
            return intron;
        }

        loop {
            if intron.next_stop == intron.next.stop {
                intron.disr = Disruption {
                    prev_hsp_idx: None,
                    next_hsp_idx: None,
                    prev_start: 0,
                    next_stop: 0,
                    prev_slen: NO_INDEX,
                    intron: false,
                    s_stop_codon: false,
                    q_interval: Interval::default(),
                    s_interval: Interval::default(),
                };
                intron.score = f64::MAX;
                return intron;
            }
            intron.disr = intron.disruption(hsps);
            if intron.disr.s_interval.valid() {
                break;
            }
            intron.next_stop += 1;
        }
        intron.disr.qc();
        intron
    }

    pub fn q_interval(&self, hsps: &[Hsp]) -> Interval {
        Interval::new(
            hsps[self.prev.hsp_idx].pos2q(self.prev_start, true),
            hsps[self.next.hsp_idx].pos2q(self.next_stop, true),
            0,
        )
    }

    pub fn s_interval(&self, hsps: &[Hsp]) -> Interval {
        let prev = &hsps[self.prev.hsp_idx];
        let next = &hsps[self.next.hsp_idx];
        let start = prev.pos2s(self.prev_start, true);
        let stop = next.pos2s(self.next_stop, true);
        if prev.s_int.strand == -1 {
            Interval::new(stop, start, prev.s_int.strand)
        } else {
            Interval::new(start, stop, prev.s_int.strand)
        }
    }

    pub fn disruption(&self, hsps: &[Hsp]) -> Disruption {
        Disruption {
            prev_hsp_idx: Some(self.prev.hsp_idx),
            next_hsp_idx: Some(self.next.hsp_idx),
            prev_start: self.prev_start,
            next_stop: self.next_stop,
            prev_slen: hsps[self.prev.hsp_idx].slen,
            intron: true,
            s_stop_codon: false,
            q_interval: self.q_interval(hsps),
            s_interval: self.s_interval(hsps),
        }
    }

    pub fn get_total_score(
        &self,
        exons: &[Exon],
        introns: &[Intron],
        hsps: &[Hsp],
        memo: &mut [Option<(AlignScore, Vec<Exon>)>],
        bacteria: bool,
        intron_score: AlignScore,
        sm: Option<&SubstMat>,
    ) -> AlignScore {
        assert!(intron_score >= 0.0);

        if self.score == f64::MAX {
            return -f64::MAX;
        }

        let (next_total_score, _) =
            self.next
                .merge_tail(exons, introns, hsps, memo, bacteria, intron_score, sm);
        next_total_score - self.score - self.disr.get_len().min(intron_score as usize) as AlignScore
    }
}

pub fn reverse_dna(seq: &mut String) -> &mut String {
    if seq.is_empty() {
        return seq;
    }

    let mut chars: Vec<char> = seq.chars().collect();
    let len = chars.len();
    for i in 0..len / 2 {
        let j = len - 1 - i;
        if chars[i] != '-' {
            chars[i] = complementary_nucleotide(chars[i]);
        }
        if chars[j] != '-' {
            chars[j] = complementary_nucleotide(chars[j]);
        }
        chars.swap(i, j);
    }

    if len % 2 == 1 {
        let i = len / 2;
        if chars[i] != '-' {
            chars[i] = complementary_nucleotide(chars[i]);
        }
    }

    *seq = chars.into_iter().collect();
    seq
}

impl Default for Hsp {
    fn default() -> Self {
        Hsp {
            merged: false,
            q_prot: false,
            s_prot: false,
            a_prot: false,
            a2q: 1,
            a2s: 1,
            qseqid: String::new(),
            sseqid: String::new(),
            q_int: Interval::default(),
            s_int: Interval::default(),
            qlen: NO_INDEX,
            slen: NO_INDEX,
            qseq: String::new(),
            sseq: String::new(),
            length: NO_INDEX,
            nident: NO_INDEX,
            qgap: NO_INDEX,
            sgap: NO_INDEX,
            qx: NO_INDEX,
            sx: NO_INDEX,
            pos2q_: Vec::new(),
            pos2s_: Vec::new(),
            sframe: 0,
            c_complete: None,
            s_internal_stop: false,
            disrs: Vec::new(),
        }
    }
}
