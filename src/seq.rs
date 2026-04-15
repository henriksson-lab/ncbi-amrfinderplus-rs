// Sequence analysis types — port of seq.hpp
// Focus on types used by amr_report, amrfinder, dna_mutation

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

/// Check if two nucleotides match (considering IUPAC ambiguity codes)
pub fn nucleotide_match(a: char, b: char) -> bool {
    if a == '-' || b == '-' {
        return false;
    }
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    if a == b {
        return true;
    }
    // Check IUPAC ambiguity
    let a_set = iupac_nucleotide_set(a);
    let b_set = iupac_nucleotide_set(b);
    // Match if sets overlap
    a_set & b_set != 0
}

fn iupac_nucleotide_set(c: char) -> u8 {
    // Bit flags: A=1, C=2, G=4, T=8
    match c {
        'A' => 0b0001,
        'C' => 0b0010,
        'G' => 0b0100,
        'T' | 'U' => 0b1000,
        'R' => 0b0101, // A or G
        'Y' => 0b1010, // C or T
        'S' => 0b0110, // G or C
        'W' => 0b1001, // A or T
        'K' => 0b1100, // G or T
        'M' => 0b0011, // A or C
        'B' => 0b1110, // C, G, or T
        'D' => 0b1101, // A, G, or T
        'H' => 0b1011, // A, C, or T
        'V' => 0b0111, // A, C, or G
        'N' => 0b1111, // any
        _ => 0,
    }
}

/// Check if two amino acids match
pub fn aa_match(a: char, b: char) -> bool {
    if a == '-' || b == '-' {
        return false;
    }
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    if a == b {
        return true;
    }
    // X matches anything
    if a == 'X' || b == 'X' {
        return true;
    }
    // B = D or N
    if (a == 'B' && (b == 'D' || b == 'N')) || (b == 'B' && (a == 'D' || a == 'N')) {
        return true;
    }
    // Z = E or Q
    if (a == 'Z' && (b == 'E' || b == 'Q')) || (b == 'Z' && (a == 'E' || a == 'Q')) {
        return true;
    }
    // J = I or L
    if (a == 'J' && (b == 'I' || b == 'L')) || (b == 'J' && (a == 'I' || a == 'L')) {
        return true;
    }
    false
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
        Interval { start, stop, strand }
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
        if (self.strand == -1) == upstream {
            self.start
        } else {
            seq_len - self.stop
        }
    }

    /// Format as 1-based for display
    pub fn format_display(&self) -> String {
        if self.strand == -1 {
            format!("{}-{}", self.stop, self.start + 1)
        } else {
            format!("{}-{}", self.start + 1, self.stop)
        }
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.start
            .cmp(&other.start)
            .then(self.stop.cmp(&other.stop))
            .then(self.strand.cmp(&other.strand))
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
    pub const NAMES: &[&str] = &["NONE", "SMOOTH", "FRAMESHIFT", "DELETION", "INSERTION"];

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "NONE" => Some(DisruptionType::None),
            "SMOOTH" => Some(DisruptionType::Smooth),
            "FRAMESHIFT" => Some(DisruptionType::Frameshift),
            "DELETION" => Some(DisruptionType::Deletion),
            "INSERTION" => Some(DisruptionType::Insertion),
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
    pub prev_start: usize,           // position in prev qseq/sseq
    pub next_stop: usize,            // position in next qseq/sseq
    pub intron: bool,
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
    pub q_prot: bool,  // reference is protein
    pub s_prot: bool,  // user input is protein (false for blastx/tblastn DNA)
    pub a_prot: bool,  // alignment is in protein space

    // Unit conversion factors (1 or 3 for protein-to-DNA)
    pub a2q: usize,
    pub a2s: usize,

    // BLAST fields — see struct doc for q/s convention
    pub qseqid: String,  // reference identifier (pipe-delimited metadata)
    pub sseqid: String,  // user's sequence identifier
    pub q_int: Interval,  // alignment interval on reference
    pub s_int: Interval,  // alignment interval on user's sequence
    pub qlen: usize,      // reference sequence length
    pub slen: usize,      // user's sequence length
    pub qseq: String,     // aligned reference sequence (with gaps)
    pub sseq: String,     // aligned user sequence (with gaps)

    // Computed by finish_hsp()
    pub length: usize,
    pub nident: usize,
    pub qgap: usize,
    pub sgap: usize,
    pub qx: usize,
    pub sx: usize,

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

    pub fn s_abs_coverage(&self) -> usize {
        self.s_int.len()
    }

    pub fn rel_identity(&self) -> f64 {
        self.nident as f64 / self.length as f64
    }

    pub fn q_rel_coverage(&self) -> f64 {
        self.q_abs_coverage() as f64 / self.qlen as f64
    }

    pub fn s_rel_coverage(&self) -> f64 {
        self.s_abs_coverage() as f64 / self.slen as f64
    }

    pub fn q_complete(&self) -> bool {
        self.q_abs_coverage() == self.qlen && self.c_complete != Some(false)
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

    /// Compute alignment statistics from qseq/sseq.
    /// Must be called after construction.
    pub fn finish_hsp(&mut self, q_stop_codon: bool, _bacterial_start_codon: bool) {
        // Handle stop codon at query end
        if q_stop_codon && self.a_prot && !self.qseq.is_empty() {
            let q_bytes = self.qseq.as_bytes();
            let s_bytes = self.sseq.as_bytes();
            let last = q_bytes.len() - 1;
            if q_bytes[last] == b'*'
                || (s_bytes[last] == b'*' && q_bytes[last] == b'-')
            {
                self.qseq.pop();
                self.sseq.pop();
                self.qlen = self.qlen.saturating_sub(1);
                self.q_int.stop = self.q_int.stop.saturating_sub(self.a2q);
                self.c_complete = Some(true);
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
            }
            if sc == '-' {
                self.sgap += 1;
            }

            if self.char_match(i) {
                self.nident += 1;
            }

            // Count ambiguous characters
            if self.a_prot {
                if "XBZJUO".contains(qc.to_ascii_uppercase()) {
                    self.qx += 1;
                }
                if "XBZJUO".contains(sc.to_ascii_uppercase()) {
                    self.sx += 1;
                }
            } else {
                if "BDHKMNRSVWY".contains(qc.to_ascii_uppercase()) {
                    self.qx += 1;
                }
                if "BDHKMNRSVWY".contains(sc.to_ascii_uppercase()) {
                    self.sx += 1;
                }
            }
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
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            anyhow::bail!("BLAST line needs at least 10 fields");
        }

        // Format: qseqid sseqid qstart qend qlen sstart send slen qseq sseq
        let qseqid = fields[0].to_string();
        let sseqid = fields[1].to_string();
        let qstart: usize = fields[2].parse()?;
        let qend: usize = fields[3].parse()?;
        let qlen: usize = fields[4].parse()?;
        let sstart: usize = fields[5].parse()?;
        let send: usize = fields[6].parse()?;
        let slen: usize = fields[7].parse()?;
        let qseq = fields[8].to_uppercase();
        let sseq = fields[9].to_uppercase();

        let s_strand: Strand = if sstart <= send { 1 } else { -1 };
        let (sstart, send) = if s_strand == -1 {
            (send, sstart)
        } else {
            (sstart, send)
        };

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
            sframe: 0,
            c_complete: None,
            s_internal_stop: false,
            disrs: Vec::new(),
        };

        hsp.finish_hsp(q_stop_codon, bacterial_start_codon);
        Ok(hsp)
    }
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
            sframe: 0,
            c_complete: None,
            s_internal_stop: false,
            disrs: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nucleotide_match() {
        assert!(nucleotide_match('A', 'A'));
        assert!(nucleotide_match('a', 'A'));
        assert!(!nucleotide_match('A', 'C'));
        assert!(nucleotide_match('N', 'A')); // N matches anything
        assert!(nucleotide_match('R', 'A')); // R = A or G
        assert!(nucleotide_match('R', 'G'));
        assert!(!nucleotide_match('R', 'C'));
        assert!(!nucleotide_match('-', 'A'));
    }

    #[test]
    fn test_aa_match() {
        assert!(aa_match('A', 'A'));
        assert!(aa_match('X', 'M')); // X matches anything
        assert!(aa_match('B', 'D')); // B = D or N
        assert!(aa_match('B', 'N'));
        assert!(!aa_match('B', 'A'));
        assert!(aa_match('Z', 'E')); // Z = E or Q
        assert!(aa_match('J', 'I')); // J = I or L
        assert!(!aa_match('-', 'A'));
    }

    #[test]
    fn test_interval() {
        let i = Interval::new(10, 20, 1);
        assert_eq!(i.len(), 10);
        assert!(!i.empty());
        assert!(i.valid());

        let j = Interval::new(12, 18, 1);
        assert!(i.contains(&j));
        assert!(i.overlaps(&j));
        assert!(!j.contains(&i));
    }

    #[test]
    fn test_strand2char() {
        assert_eq!(strand2char(1), '+');
        assert_eq!(strand2char(-1), '-');
        assert_eq!(strand2char(0), '?');
    }
}
