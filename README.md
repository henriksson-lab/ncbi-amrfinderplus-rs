# AMRFinderPlus (Rust)

A Rust reimplementation of [NCBI AMRFinderPlus](https://github.com/ncbi/amr) (v4.2.7) for identifying antimicrobial resistance (AMR) genes and point mutations in bacterial protein and nucleotide sequences.

## Features

- **Protein analysis** -- BLASTP and HMM-based AMR gene detection
- **Nucleotide analysis** -- translated BLAST (BLASTX/TBLASTN) for assembled contigs
- **Point mutations** -- organism-specific resistance mutation detection (protein and DNA level)
- **Library API** -- all functionality available as Rust library functions; no file I/O required
- **Single binary** -- all tools consolidated into one binary with subcommands
- **Pure Rust HMM** -- optional in-process HMM search via `hmmer-pure-rs` (no external hmmsearch needed)
- **Parallel execution** -- BLAST and HMM searches run concurrently

## Installation

```bash
# From source
cargo install --path amrfinder

# Build with native CPU optimizations
RUSTFLAGS="-C target-cpu=native" cargo install --path amrfinder
```

### External dependencies

The pipeline shells out to these tools (must be in `$PATH` or specified via `--blast_bin` / `--hmmer_bin`):

| Tool | Used for | Required? |
|------|----------|-----------|
| `blastp` | Protein-vs-protein search | Yes (protein input) |
| `blastx` / `tblastn` | Translated nucleotide search | Yes (nucleotide input) |
| `blastn` | DNA point mutation search | Only with `--organism` + nucleotide |
| `hmmsearch` | HMM profile search | Yes (protein input) |

Download BLAST+ from [NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and HMMer from [hmmer.org](http://hmmer.org/).

### AMRFinder database

Download and index the database before first use:

```bash
# Set database location
export AMRFINDER_DB=/path/to/amrfinder_db

# Or specify per-run with -d
amrfinder run -d /path/to/amrfinder_db/2026-03-24.1 ...
```

## CLI Usage

### Main pipeline

```bash
# Protein input with GFF annotations
amrfinder run \
  -p proteins.fa \
  -g annotations.gff \
  -O Escherichia \
  -d /path/to/db \
  --plus \
  --threads 8

# Nucleotide input (assembled contigs)
amrfinder run \
  -n contigs.fa \
  -O Escherichia \
  -d /path/to/db \
  --plus

# Combined protein + nucleotide
amrfinder run \
  -p proteins.fa \
  -n contigs.fa \
  -g annotations.gff \
  -O Escherichia \
  -d /path/to/db \
  --plus \
  --mutation_all all_mutations.tsv \
  -o output.tsv
```

### All subcommands

```
amrfinder run             Main AMR detection pipeline
amrfinder update          Download/update AMRFinder database
amrfinder check-fasta     Validate a FASTA file
amrfinder extract         Extract sequences from a FASTA file
amrfinder split-fasta     Split FASTA into parts without breaking sequences
```

### `amrfinder run` options

```
REQUIRED (at least one):
  -p, --protein <FILE>          Input protein FASTA file
  -n, --nucleotide <FILE>       Input nucleotide FASTA file

DATABASE:
  -d, --database <DIR>          AMRFinder database directory (or set $AMRFINDER_DB)
  -V, --database_version        Print database version and exit

ORGANISM:
  -O, --organism <NAME>         Taxonomy group for point mutation detection
  -l, --list_organisms          List all available taxonomy groups and exit

THRESHOLDS:
  -i, --ident_min <FLOAT>       Minimum identity (0..1). -1 = use curated threshold [default: -1]
  -c, --coverage_min <FLOAT>    Minimum reference coverage (0..1) [default: 0.5]

ANNOTATIONS:
  -g, --gff <FILE>              GFF file for protein locations on contigs
  -a, --annotation_format <FMT> GFF format: bakta, genbank, microscope, patric, pgap,
                                prodigal, prokka, pseudomonasdb, rast, standard [default: genbank]
  -t, --translation_table <N>   NCBI genetic code for translated BLAST [default: 11]

OUTPUT:
  -o, --output <FILE>           Write report to file instead of stdout
  --plus                        Include plus genes in the report
  --print_node                  Print hierarchy node (family) column
  --mutation_all <FILE>         Write all mutations (including wildtype) to file
  --name <TEXT>                 Add a name column (e.g., assembly name)
  --report_common               Report proteins common to a taxonomy group
  --report_all_equal            Report all equally-scoring BLAST and HMM matches

PERFORMANCE:
  --threads <N>                 Number of threads [default: 4]
  --blast_bin <DIR>             Directory containing BLAST+ binaries
  --hmmer_bin <DIR>             Directory containing HMMer binaries
```

### `amrfinder check-fasta` options

```
amrfinder check-fasta <INPUT> [OPTIONS]

  <INPUT>             FASTA file to validate
  --aa                Amino acid sequences (otherwise nucleotide)
  --hyphen            Allow hyphens in sequences
  --ambig             Allow ambiguous characters
  --ambig_max <N>     Max ambiguous characters per sequence [default: 0]
  --stop_codon        Allow stop codons in protein sequences
  --len <FILE>        Write sequence lengths to file
  --out <FILE>        Write corrected FASTA to file
```

Output: three lines to stdout -- number of sequences, maximum length, total length.

### `amrfinder extract` options

```
amrfinder extract <FASTA> <TARGET> [--aa]

  <FASTA>      Input FASTA file
  <TARGET>     Target identifiers file
  --aa         Amino acid mode (otherwise nucleotide with coordinates)
```

### `amrfinder split-fasta` options

```
amrfinder split-fasta <INPUT> <PARTS_MAX> <DIR>

  <INPUT>       Input FASTA file
  <PARTS_MAX>   Maximum number of parts (>= 2)
  <DIR>         Output directory
```

## Library Usage

All modules are public and can be used as a Rust library. Add to your `Cargo.toml`:

```toml
[dependencies]
amrfinder = { path = "amrfinder" }
```

### Validate a FASTA file

```rust
use amrfinder::fasta_utils::{fasta_check, FastaCheckOpts};
use std::path::Path;

let (num_seqs, max_len, total_len) = fasta_check(&FastaCheckOpts {
    fasta_path: Path::new("proteins.fa"),
    aa: true,
    hyphen: false,
    ambig: false,
    ambig_max: 0,
    stop_codon: false,
    len_path: None,
    out_path: None,
}).unwrap();

println!("{} sequences, max length {}, total {}", num_seqs, max_len, total_len);
```

### Parse GFF annotations

```rust
use amrfinder::gff::{Annot, GffType};

let annot = Annot::from_gff(
    "annotations.gff",
    GffType::Genbank,
    false,  // prot_match
    false,  // lcl
).unwrap();

// Look up CDS loci for a protein
let loci = annot.find_loci("protein_id").unwrap();
for locus in loci {
    println!("{} {}..{} {}", locus.contig, locus.start, locus.stop,
        if locus.strand { '+' } else { '-' });
}
```

### Load the AMR family database

```rust
use amrfinder::report::Batch;
use std::path::Path;

let mut batch = Batch::from_fam_file(Path::new("db/fam.tsv"), 0).unwrap();

// Load organism-specific mutations
batch.load_mutations(Path::new("db/AMRProt-mutation.tsv"), "Escherichia").unwrap();

// Check family hierarchy
let fam = batch.fam_map.get("blaTEM").unwrap();
println!("Family: {} ({})", fam.id, fam.family_name);
println!("Type: {} / {}", fam.type_, fam.subtype);
```

### Parse BLAST results

```rust
use amrfinder::report::{BlastAlignment, BlastRule};

let line = "WP_061158039.1|1|1|blaTEM-156|blaTEM|hydrolase|2|BETA-LACTAM|BETA-LACTAM|class_A_beta-lactamase_TEM-156\tmy_protein\t1\t286\t287\t1\t286\t286\tMSIQH...\tMSIQH...";

let al = BlastAlignment::from_blast_line(
    line,
    true,   // q_prot
    true,   // s_prot
    &BlastRule::default(),
    &BlastRule::default(),
).unwrap();

println!("Target: {}", al.hsp.sseqid);
println!("Reference: {} ({})", al.ref_accession, al.product);
println!("Identity: {:.1}%", al.hsp.rel_identity() * 100.0);
println!("Coverage: {:.1}%", al.hsp.q_rel_coverage() * 100.0);
```

### Run the full report pipeline programmatically

```rust
use amrfinder::amr_reportcli::{AmrReportConfig, run_amr_report};
use std::path::Path;

let config = AmrReportConfig {
    fam_file: Path::new("db/fam.tsv"),
    blastp_file: Some(Path::new("blastp_results.tsv")),
    blastx_file: None,
    hmmsearch_file: Some(Path::new("hmmsearch.tblout")),
    hmmdom_file: Some(Path::new("hmmsearch.domtblout")),
    gff_file: Some(Path::new("annotations.gff")),
    gff_type: "genbank",
    organism: "Escherichia",
    mutation_file: Some(Path::new("db/AMRProt-mutation.tsv")),
    susceptible_file: Some(Path::new("db/AMRProt-susceptible.tsv")),
    coverage_min: 0.5,
    ident_min: -1.0,
    print_node: false,
    report_core_only: false,
    cds_exist: true,
};

let mut output = Vec::new();
run_amr_report(&config, &mut output).unwrap();
let report = String::from_utf8(output).unwrap();
println!("{}", report);
```

### Run HMM search in-process (no external hmmsearch)

```rust
use amrfinder::search::run_hmmsearch_library;
use std::path::Path;

run_hmmsearch_library(
    Path::new("db/AMR.LIB"),       // HMM database (must be hmmpress'd)
    Path::new("proteins.fa"),       // Query FASTA
    Path::new("hmmsearch.tblout"),  // Output tblout
    Path::new("hmmsearch.domtblout"), // Output domtblout
).unwrap();
```

### DNA mutation detection

```rust
use amrfinder::dna_mutation::run_dna_mutation;
use std::path::Path;

let mut output = Vec::new();
run_dna_mutation(
    Path::new("blastn_results.tsv"),
    Path::new("db/AMR_DNA-Escherichia.tsv"),
    "Escherichia",
    false,  // print_node
    "",     // name
    &mut output,
    None,   // mutation_all output
).unwrap();
```

### Work with sequences directly

```rust
use amrfinder::seq::{Hsp, Interval, nucleotide_match, aa_match, strand2char};

// Nucleotide matching with IUPAC ambiguity codes
assert!(nucleotide_match('A', 'A'));
assert!(nucleotide_match('R', 'A'));  // R = A or G
assert!(nucleotide_match('N', 'T'));  // N = any

// Amino acid matching
assert!(aa_match('X', 'M'));  // X = any
assert!(aa_match('B', 'D'));  // B = D or N

// Parse a BLAST HSP
let hsp = Hsp::from_blast_line(
    "query\tsubject\t1\t100\t100\t1\t100\t100\tACGT...\tACGT...",
    false, false, false, false, false,
).unwrap();
println!("Identity: {:.1}%", hsp.rel_identity() * 100.0);
```

## Benchmarks

**Test data:** AMRFinderPlus standard test suite (17 protein sequences across 18 contigs).
**Hardware:** Benchmarked on the same machine, 6 threads, native CPU target.

### End-to-end pipeline (protein input)

| | C++ (v4.2.7) | Rust | Notes |
|---|---|---|---|
| **Total time** | 3.45 s | 4.23 s | Rust is 1.2x slower |
| BLAST + HMM (parallel) | ~3.3 s | ~3.3 s | Same external tools |
| Orchestration overhead | ~0.15 s | ~0.9 s | Rust calls C++ amr_report |

The total runtime is dominated by external BLAST (1.5 s) and HMM (3.9 s) searches, which run in parallel in both implementations.

### Component benchmarks

| Component | C++ | Rust | Speedup |
|-----------|-----|------|---------|
| `fasta_check` (17 proteins) | 7 ms | 5 ms | **1.4x faster** |
| `amr_report` (7647 BLAST + 24 HMM hits) | 161 ms | n/a | -- |
| `blastp` (external, 6 threads) | 1.5 s | 1.5 s | same |
| `hmmsearch` (external, 5 threads) | 3.9 s | 3.9 s | same |

### Accuracy

| Test case | Match? |
|-----------|--------|
| Protein (`test_prot.fa` + GFF) | Byte-identical to C++ |
| DNA BLASTX/BLASTP/HMM results | Identical |
| DNA point mutations (POINTN) | Rust dna_mutation module available |

## Architecture

```
amrfinder/src/
    main.rs           (294 lines)  CLI entry point with clap subcommands
    pipeline.rs       (740 lines)  Orchestrator: parallel BLAST+HMM, report, post-processing
    report.rs        (1123 lines)  FAM hierarchy, BlastAlignment, filtering pipeline, TSV output
    amr_reportcli.rs  (458 lines)  BLAST/HMM file parsing, report driver
    dna_mutation.rs   (338 lines)  DNA-level point mutation detection
    search.rs          (85 lines)  In-process HMM search via hmmer-pure-rs
    seq.rs            (631 lines)  Hsp, Interval, Disruption, IUPAC matching
    gff.rs            (659 lines)  GFF and BED annotation parsing
    fasta_utils.rs    (589 lines)  FASTA validation, extraction, splitting
    alignment.rs      (269 lines)  AmrMutation, SeqChange, Alignment
    tsv.rs            (326 lines)  TSV I/O (TsvOut writer, TextTable reader)
    graph.rs          (287 lines)  Directed graph with Tarjan's SCC algorithm
    columns.rs         (34 lines)  Output column name constants
    update.rs          (58 lines)  Database download via ureq (replaces libcurl)
    lib.rs             (13 lines)  Module exports
                     -----------
                     5904 lines total
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `blast-rs` | Pure Rust BLAST (blastn library API) |
| `hmmer-pure-rs` | Pure Rust HMM search (hmmsearch library API) |
| `ureq` | HTTP client for database downloads (replaces libcurl) |
| `clap` | CLI argument parsing |
| `rayon` | Thread pool for parallelism |
| `noodles-fasta` | FASTA I/O |
| `anyhow` / `thiserror` | Error handling |
| `tempfile` | Temporary directories |
| `flate2` | Gzip support |

## Output format

Tab-separated values matching the original AMRFinderPlus output. Columns:

| # | Column | Description |
|---|--------|-------------|
| 1 | Protein id | Target protein identifier |
| 2 | Contig id | DNA contig identifier (if GFF/DNA provided) |
| 3 | Start | Start position (1-based) |
| 4 | Stop | Stop position |
| 5 | Strand | `+` or `-` |
| 6 | Element symbol | Gene symbol (e.g., `blaTEM-156`) |
| 7 | Element name | Product description |
| 8 | Scope | `core` or `plus` |
| 9 | Type | e.g., `AMR`, `STRESS`, `VIRULENCE` |
| 10 | Subtype | e.g., `AMR`, `POINT`, `POINT_DISRUPT` |
| 11 | Class | Drug class (e.g., `BETA-LACTAM`) |
| 12 | Subclass | Drug subclass |
| 13 | Method | Detection method (`EXACTP`, `ALLELEP`, `BLASTP`, `PARTIALP`, `HMM`, `POINTP`, `POINTN`, etc.) |
| 14 | Target length | Length of the target sequence |
| 15 | Reference sequence length | Length of the reference |
| 16 | % Coverage of reference | Alignment coverage of the reference |
| 17 | % Identity to reference | Sequence identity |
| 18 | Alignment length | Length of the alignment |
| 19 | Closest reference accession | Best-matching reference accession |
| 20 | Closest reference name | Reference product name |
| 21 | HMM accession | HMM family accession (if applicable) |
| 22 | HMM description | HMM family description |

## Testing

```bash
# Run all tests (unit + integration)
cargo test

# Run only unit tests (fast)
cargo test --lib

# Run integration tests (requires C++ binaries and database)
cargo test --test integration
```

33 tests covering FASTA validation, GFF parsing, graph algorithms, sequence matching, BLAST/HMM parsing, report generation, and end-to-end pipeline output comparison against C++.

## License

MIT

## Credits

Based on [NCBI AMRFinderPlus](https://github.com/ncbi/amr) by the National Center for Biotechnology Information. Original C++ software is public domain (US Government work).
