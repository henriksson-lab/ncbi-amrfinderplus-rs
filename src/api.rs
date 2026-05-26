use std::path::PathBuf;

use anyhow::Result;

use crate::pipeline::{run_pipeline, PipelineConfig};

/// High-level AMRFinder pipeline entry point.
///
/// Use [`AmrFinder::builder`] for a stable, chainable crate API instead of
/// constructing [`PipelineConfig`] directly.
pub struct AmrFinder;

impl AmrFinder {
    pub fn builder() -> AmrFinderBuilder {
        AmrFinderBuilder::new()
    }
}

/// Result returned by the high-level pipeline API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmrFinderRun {
    /// Main AMRFinder TSV report.
    pub report: String,
    /// Mutation-all side-output path, when requested.
    pub mutation_all: Option<PathBuf>,
    /// Main report output path, when requested.
    pub output: Option<PathBuf>,
}

/// Builder for running the AMRFinder pipeline from library code.
#[derive(Debug, Clone)]
pub struct AmrFinderBuilder {
    config: PipelineConfig,
}

impl Default for AmrFinderBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl AmrFinderBuilder {
    pub fn new() -> Self {
        Self {
            config: PipelineConfig::default(),
        }
    }

    pub fn protein(mut self, path: impl Into<PathBuf>) -> Self {
        self.config.protein = Some(path.into());
        self
    }

    pub fn nucleotide(mut self, path: impl Into<PathBuf>) -> Self {
        self.config.nucleotide = Some(path.into());
        self
    }

    pub fn gff(mut self, path: impl Into<PathBuf>) -> Self {
        self.config.gff = Some(path.into());
        self
    }

    pub fn database(mut self, path: impl Into<PathBuf>) -> Self {
        self.config.database = path.into();
        self
    }

    pub fn organism(mut self, organism: impl Into<String>) -> Self {
        self.config.organism = organism.into();
        self
    }

    pub fn ident_min(mut self, ident_min: f64) -> Self {
        self.config.ident_min = ident_min;
        self
    }

    pub fn coverage_min(mut self, coverage_min: f64) -> Self {
        self.config.coverage_min = coverage_min;
        self
    }

    pub fn threads(mut self, threads: usize) -> Self {
        self.config.threads = threads;
        self
    }

    pub fn plus(mut self, plus: bool) -> Self {
        self.config.plus = plus;
        self
    }

    pub fn report_common(mut self, report_common: bool) -> Self {
        self.config.report_common = report_common;
        self
    }

    pub fn report_all_equal(mut self, report_all_equal: bool) -> Self {
        self.config.report_all_equal = report_all_equal;
        self
    }

    pub fn print_node(mut self, print_node: bool) -> Self {
        self.config.print_node = print_node;
        self
    }

    pub fn mutation_all(mut self, path: impl Into<PathBuf>) -> Self {
        self.config.mutation_all = Some(path.into());
        self
    }

    pub fn annotation_format(mut self, annotation_format: impl Into<String>) -> Self {
        self.config.annotation_format = annotation_format.into();
        self
    }

    pub fn translation_table(mut self, translation_table: u32) -> Self {
        self.config.translation_table = translation_table;
        self
    }

    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.config.name = name.into();
        self
    }

    pub fn blast_bin(mut self, dir: impl Into<String>) -> Self {
        self.config.blast_bin = dir.into();
        self
    }

    pub fn hmmer_bin(mut self, dir: impl Into<String>) -> Self {
        self.config.hmmer_bin = dir.into();
        self
    }

    pub fn output(mut self, path: impl Into<PathBuf>) -> Self {
        self.config.output = Some(path.into());
        self
    }

    pub fn config(&self) -> &PipelineConfig {
        &self.config
    }

    pub fn into_config(self) -> PipelineConfig {
        self.config
    }

    pub fn run(self) -> Result<AmrFinderRun> {
        let mutation_all = self.config.mutation_all.clone();
        let output = self.config.output.clone();
        let report = run_pipeline(&self.config)?;

        Ok(AmrFinderRun {
            report,
            mutation_all,
            output,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builder_sets_pipeline_config_options() {
        let config = AmrFinder::builder()
            .protein("proteins.fa")
            .nucleotide("contigs.fa")
            .gff("annot.gff")
            .database("db")
            .organism("Escherichia")
            .ident_min(0.95)
            .coverage_min(0.8)
            .threads(8)
            .plus(true)
            .report_common(true)
            .print_node(true)
            .mutation_all("mutation_all.tsv")
            .annotation_format("gff3")
            .translation_table(4)
            .name("sample")
            .blast_bin("/opt/blast/bin")
            .hmmer_bin("/opt/hmmer/bin")
            .output("report.tsv")
            .into_config();

        assert_eq!(config.protein, Some(PathBuf::from("proteins.fa")));
        assert_eq!(config.nucleotide, Some(PathBuf::from("contigs.fa")));
        assert_eq!(config.gff, Some(PathBuf::from("annot.gff")));
        assert_eq!(config.database, PathBuf::from("db"));
        assert_eq!(config.organism, "Escherichia");
        assert_eq!(config.ident_min, 0.95);
        assert_eq!(config.coverage_min, 0.8);
        assert_eq!(config.threads, 8);
        assert!(config.plus);
        assert!(config.report_common);
        assert!(config.print_node);
        assert_eq!(config.mutation_all, Some(PathBuf::from("mutation_all.tsv")));
        assert_eq!(config.annotation_format, "gff3");
        assert_eq!(config.translation_table, 4);
        assert_eq!(config.name, "sample");
        assert_eq!(config.blast_bin, "/opt/blast/bin");
        assert_eq!(config.hmmer_bin, "/opt/hmmer/bin");
        assert_eq!(config.output, Some(PathBuf::from("report.tsv")));
    }
}
