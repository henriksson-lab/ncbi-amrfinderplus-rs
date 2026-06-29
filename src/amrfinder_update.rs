// Database download and update skeleton, shaped after upstream amrfinder_update.cpp.

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::common::{
    exec, shell_quote, which, Application, DataVersion, Key, Run, SoftwareVersion,
};
use crate::curl_easy::Curl;

const NCBI_FTP_BASE: &str =
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database";
const LATEST_DIR: &str = "latest";
const VERSION_FILE: &str = "version.txt";
const CORE_DB_FILES: &[&str] = &[
    "AMR.LIB",
    "AMRProt.fa",
    "AMRProt-mutation.tsv",
    "AMRProt-suppress.tsv",
    "AMRProt-susceptible.fa",
    "AMRProt-susceptible.tsv",
    "AMR_CDS.fa",
    "database_format_version.txt",
    "fam.tsv",
    "taxgroup.tsv",
    VERSION_FILE,
];
const CHANGES_FILE: &str = "changes.txt";

#[derive(Debug, Clone)]
pub struct UpdateOptions {
    pub database_dir: PathBuf,
    pub force: bool,
    pub source_base: String,
    pub current_minor: String,
    pub amrfinder_index: Option<AmrFinderIndexOptions>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmrFinderIndexOptions {
    pub program: PathBuf,
    pub blast_bin: Option<PathBuf>,
    pub hmmer_bin: Option<PathBuf>,
    pub quiet: bool,
    pub debug: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmrFinderIndexInvocation {
    pub version_dir: PathBuf,
    pub program: PathBuf,
    pub blast_bin: Option<PathBuf>,
    pub hmmer_bin: Option<PathBuf>,
    pub quiet: bool,
    pub debug: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UpdateStatus {
    AlreadyCurrent,
    Updated,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UpdateOutcome {
    pub status: UpdateStatus,
    pub latest_minor: String,
    pub latest_data_version: String,
    pub fetched_files: Vec<PathBuf>,
    pub warnings: Vec<String>,
}

pub struct ThisApplication {
    application: Application,
    cur_minor: String,
    amrfinder_index: PathBuf,
    source_base: String,
}

impl ThisApplication {
    pub fn new() -> Self {
        Self {
            application: Application {
                description: concat!(
                    "Update the database for AMRFinder from ",
                    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/",
                    "\nRequirement: the database directory contains subdirectories named by database versions."
                ),
                usage: "amrfinder_update [-database DATABASE_DIR] [-blast_bin BLAST_DIR] [-hmmer_bin HMMER_DIR] [-force_update]",
                positionals: 0,
                needs_arg: false,
                threads_used: false,
                key_args: vec![
                    Key {
                        name: "database",
                        flag: false,
                        default_value: "$BASE/data",
                        acronym: Some('d'),
                    },
                    Key {
                        name: "blast_bin",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "hmmer_bin",
                        flag: false,
                        default_value: "",
                        acronym: None,
                    },
                    Key {
                        name: "force_update",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                    Key {
                        name: "quiet",
                        flag: true,
                        default_value: "false",
                        acronym: Some('q'),
                    },
                    Key {
                        name: "debug",
                        flag: true,
                        default_value: "false",
                        acronym: None,
                    },
                ],
            },
            cur_minor: current_software_minor(),
            amrfinder_index: PathBuf::from("amrfinder_index"),
            source_base: NCBI_FTP_BASE.to_string(),
        }
    }

    pub fn create_latest_link(&self, main_dir: &Path, latest_dir: &str) -> Result<()> {
        create_latest_link(main_dir, latest_dir)
    }

    pub fn shell_body(&self, run: &Run) -> Result<UpdateOutcome> {
        let database_dir = PathBuf::from(run.key_args.get("database").expect("database defaulted"));
        let blast_bin = run.key_args.get("blast_bin").expect("blast_bin defaulted");
        let blast_bin = (!blast_bin.is_empty()).then(|| PathBuf::from(blast_bin));
        let hmmer_bin = run.key_args.get("hmmer_bin").expect("hmmer_bin defaulted");
        let hmmer_bin = (!hmmer_bin.is_empty()).then(|| PathBuf::from(hmmer_bin));
        let force_update = run
            .key_args
            .get("force_update")
            .is_some_and(|value| value == "true");
        let quiet = run
            .key_args
            .get("quiet")
            .is_some_and(|value| value == "true");
        let debug = run
            .key_args
            .get("debug")
            .or_else(|| run.key_args.get("qc"))
            .is_some_and(|value| value == "true");

        shell_body(UpdateOptions {
            database_dir,
            force: force_update,
            source_base: self.source_base.clone(),
            current_minor: self.cur_minor.clone(),
            amrfinder_index: Some(AmrFinderIndexOptions {
                program: self.amrfinder_index.clone(),
                blast_bin,
                hmmer_bin,
                quiet,
                debug,
            }),
        })
    }
}

pub fn main(mut argv: Vec<OsString>) -> Result<i32> {
    let mut app = ThisApplication::new();
    if let Some(program_path) = argv.first().cloned().and_then(|program| {
        let path = Path::new(&program);
        path.canonicalize().ok().or_else(|| {
            if path
                .parent()
                .is_some_and(|parent| !parent.as_os_str().is_empty())
            {
                Some(path.to_path_buf())
            } else {
                let dir = which(&program.to_string_lossy());
                (!dir.is_empty()).then(|| {
                    let found = PathBuf::from(dir).join(path);
                    found.canonicalize().unwrap_or(found)
                })
            }
        })
    }) {
        if Path::new(&argv[0])
            .parent()
            .is_none_or(|parent| parent.as_os_str().is_empty())
        {
            argv[0] = OsString::from(&program_path);
        }
        let Some(exec_dir) = program_path.parent().map(Path::to_path_buf) else {
            return Ok(1);
        };
        app.amrfinder_index = exec_dir.join("amrfinder_index");
    }
    let run = match app.application.run(argv) {
        Ok(run) => run,
        Err(code) => return Ok(code),
    };
    app.shell_body(&run)?;
    Ok(0)
}

fn read_url(curl: &mut Curl, url: &str) -> Result<String> {
    curl.read(url)
}

/// Upstream getLatestMinor analogue: inspect the database root and choose the
/// greatest major.minor directory advertised there.
pub fn get_latest_minor(curl: &mut Curl) -> Result<String> {
    get_latest_minor_from_base(curl, NCBI_FTP_BASE)
}

/// Upstream getLatestDataVersion analogue: inspect a minor-version directory
/// and choose the greatest published data-version directory.
pub fn get_latest_data_version(curl: &mut Curl, minor: &str) -> Result<String> {
    get_latest_data_version_from_base(curl, NCBI_FTP_BASE, minor)
}

/// Upstream fetchAMRFile analogue: copy/download one database artifact.
pub fn fetch_amr_file(
    curl: &mut Curl,
    url_dir: &str,
    local_dir: &Path,
    f_name: &str,
) -> Result<()> {
    assert!(url_dir.ends_with('/'), "urlDir must be a directory name");
    assert!(
        local_dir.display().to_string().ends_with('/'),
        "localDir must be a directory name"
    );
    assert!(!f_name.is_empty());

    let url = format!("{url_dir}{f_name}");
    let dest = format!("{}{f_name}", local_dir.display());

    curl.download(&url, &dest)
}

/// Upstream shellBody analogue: orchestrate version discovery, freshness check,
/// artifact fetch, validation, and optional index building.
pub fn shell_body(options: UpdateOptions) -> Result<UpdateOutcome> {
    let mut curl = Curl::new()?;
    eprintln!(
        "Looking up the published databases at {}/",
        options.source_base.trim_end_matches('/')
    );
    let target = plan_update(&mut curl, &options.source_base, &options.current_minor)?;
    for warning in &target.warnings {
        eprintln!();
        eprintln!("Warning: {warning}");
    }

    let database_dir = dir_name_path(&options.database_dir);
    fs::create_dir_all(&database_dir)
        .with_context(|| format!("failed to create {}", database_dir.display()))?;
    let version_dir = database_dir.join(&target.data_version);
    let version_dir_name = dir_name_path(&version_dir);
    eprintln!(
        "Looking for the target directory: {}",
        version_dir_name.display()
    );

    if version_dir.is_dir() {
        if options.force {
            eprintln!(
                "{}: already exists, overwriting what was there",
                version_dir_name.display()
            );
        } else if local_version_matches_remote(
            &mut curl,
            &options.source_base,
            &target.minor,
            &target.data_version,
            &version_dir,
        )? {
            let version_path = version_dir.join(VERSION_FILE);
            let version_old = fs::read_to_string(&version_path)
                .with_context(|| format!("failed to read {}", version_path.display()))?;
            let version_old = normalize_version(&version_old)?;
            eprintln!(
                "Warning: {}: contains the latest version {version_old}",
                version_dir_name.display()
            );
            eprintln!("Skipping update");
            eprintln!("Use amrfinder --force_update to overwrite the existing database");
            create_latest_link(&database_dir, &target.data_version)?;
            return Ok(UpdateOutcome {
                status: UpdateStatus::AlreadyCurrent,
                latest_minor: target.minor,
                latest_data_version: target.data_version,
                fetched_files: Vec::new(),
                warnings: target.warnings,
            });
        }
    } else {
        fs::create_dir_all(&version_dir)
            .with_context(|| format!("failed to create {}", version_dir.display()))?;
    }

    eprintln!(
        "Downloading AMRFinder database version {} into: {}",
        target.data_version,
        version_dir_name.display()
    );

    let mut fetched_files = Vec::new();
    for file_name in CORE_DB_FILES {
        fetched_files.push(fetch_db_file(
            &mut curl,
            &options.source_base,
            &target.minor,
            &target.data_version,
            &version_dir,
            file_name,
        )?);
    }

    for taxgroup in dna_point_mutation_taxgroups(&version_dir.join("taxgroup.tsv"))? {
        fetched_files.push(fetch_db_file(
            &mut curl,
            &options.source_base,
            &target.minor,
            &target.data_version,
            &version_dir,
            &format!("AMR_DNA-{taxgroup}.fa"),
        )?);
        fetched_files.push(fetch_db_file(
            &mut curl,
            &options.source_base,
            &target.minor,
            &target.data_version,
            &version_dir,
            &format!("AMR_DNA-{taxgroup}.tsv"),
        )?);
    }

    fetched_files.push(fetch_db_file(
        &mut curl,
        &options.source_base,
        &target.minor,
        &target.data_version,
        &version_dir,
        CHANGES_FILE,
    )?);
    create_latest_link(&database_dir, &target.data_version)?;

    if let Some(index_options) = options.amrfinder_index {
        run_amrfinder_index(&AmrFinderIndexInvocation {
            version_dir: dir_name_path(&version_dir),
            program: index_options.program,
            blast_bin: index_options.blast_bin,
            hmmer_bin: index_options.hmmer_bin,
            quiet: index_options.quiet,
            debug: index_options.debug,
        })?;
    }

    Ok(UpdateOutcome {
        status: UpdateStatus::Updated,
        latest_minor: target.minor,
        latest_data_version: target.data_version,
        fetched_files,
        warnings: target.warnings,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PlannedUpdate {
    minor: String,
    data_version: String,
    warnings: Vec<String>,
}

fn plan_update(curl: &mut Curl, source_base: &str, current_minor: &str) -> Result<PlannedUpdate> {
    let published_minor = match get_latest_minor_from_base(curl, source_base) {
        Ok(published_minor) => published_minor,
        Err(err) if format!("{err:#}").contains("no major.minor directory found") => {
            bail!("Cannot get the software minor version of the latest published database version")
        }
        Err(err) => return Err(err),
    };
    let published_data_version = match get_latest_data_version_from_base(
        curl,
        source_base,
        &published_minor,
    ) {
        Ok(published_data_version) => published_data_version,
        Err(err) if format!("{err:#}").contains("no data version directory found") => {
            bail!(
                "Cannot get the latest published database version for the software minor version {published_minor}"
            )
        }
        Err(err) => return Err(err),
    };
    let current_data_version =
        latest_data_version_for_current_minor(curl, source_base, current_minor)?;

    let mut warnings = Vec::new();
    let (minor, data_version) = if let Some(current_data_version) = current_data_version {
        if current_data_version != published_data_version {
            let current = DataVersion::from_text(&current_data_version)
                .expect("current data version was parsed earlier");
            let published = DataVersion::from_text(&published_data_version)
                .expect("published data version was parsed earlier");
            assert!(
                current < published,
                "current data version must be older than published data version"
            );
            warnings.push(format!(
                "A newer version of the database exists ({published_data_version}), but it requires a newer version of the software ({published_minor}) to install.\nSee https://github.com/ncbi/amr/wiki/Upgrading for more information."
            ));
        }
        (current_minor.to_string(), current_data_version)
    } else {
        warnings.push(format!(
            "Cannot get the latest published database version for the current software minor version {current_minor}.\nThe latest published database version {published_data_version} for the latest published software minor version {published_minor} will be used instead"
        ));
        (published_minor, published_data_version)
    };

    Ok(PlannedUpdate {
        minor,
        data_version,
        warnings,
    })
}

fn current_software_minor() -> String {
    minor_from_software_version(include_str!("../amr/version.txt").trim())
}

fn minor_from_software_version(version: &str) -> String {
    let mut parts = version.split('.');
    match (parts.next(), parts.next()) {
        (Some(major), Some(minor)) if !major.is_empty() && !minor.is_empty() => {
            format!("{major}.{minor}")
        }
        _ => version.to_string(),
    }
}

fn get_latest_minor_from_base(curl: &mut Curl, source_base: &str) -> Result<String> {
    let url = latest_minor_index_url(source_base);
    let listing = read_url(curl, &url)?;
    parse_latest_minor(&listing).with_context(|| format!("failed to parse latest minor from {url}"))
}

fn get_latest_data_version_from_base(
    curl: &mut Curl,
    source_base: &str,
    minor: &str,
) -> Result<String> {
    let url = latest_data_version_index_url(source_base, minor);
    let listing = read_url(curl, &url)?;
    parse_latest_data_version(&listing)
        .with_context(|| format!("failed to parse latest data version from {url}"))
}

fn latest_data_version_for_current_minor(
    curl: &mut Curl,
    source_base: &str,
    current_minor: &str,
) -> Result<Option<String>> {
    let url = latest_data_version_index_url(source_base, current_minor);
    let listing = read_url(curl, &url)?;
    match parse_latest_data_version(&listing) {
        Ok(version) => Ok(Some(version)),
        Err(_) => Ok(None),
    }
}

fn local_version_matches_remote(
    curl: &mut Curl,
    source_base: &str,
    minor: &str,
    data_version: &str,
    version_dir: &Path,
) -> Result<bool> {
    let path = version_dir.join(VERSION_FILE);
    let tmp = tempfile::TempDir::new().context("failed to create temporary directory")?;
    let remote_path = tmp.path().join("curl");
    curl.download(
        &db_file_url(source_base, minor, data_version, VERSION_FILE),
        &remote_path.display().to_string(),
    )?;
    let local = fs::read_to_string(&path)
        .with_context(|| format!("Reading file {}", shell_quote(&path.to_string_lossy())))?;
    let remote = fs::read_to_string(&remote_path).with_context(|| {
        format!(
            "Reading file {}",
            shell_quote(&remote_path.to_string_lossy())
        )
    })?;
    match (normalize_version(&local), normalize_version(&remote)) {
        (Ok(local), Ok(remote)) => Ok(local == remote),
        _ => Ok(false),
    }
}

fn latest_minor_index_url(source_base: &str) -> String {
    dir_url(&join_url(source_base, &[]))
}

fn latest_data_version_index_url(source_base: &str, minor: &str) -> String {
    dir_url(&join_url(source_base, &[minor]))
}

fn db_file_url(source_base: &str, minor: &str, data_version: &str, file_name: &str) -> String {
    join_url(source_base, &[minor, data_version, file_name])
}

fn fetch_db_file(
    curl: &mut Curl,
    source_base: &str,
    minor: &str,
    data_version: &str,
    version_dir: &Path,
    file_name: &str,
) -> Result<PathBuf> {
    let dest = version_dir.join(file_name);
    fetch_amr_file(
        curl,
        &dir_url(&join_url(source_base, &[minor, data_version])),
        &dir_name_path(version_dir),
        file_name,
    )?;
    Ok(dest)
}

fn run_amrfinder_index(invocation: &AmrFinderIndexInvocation) -> Result<()> {
    let tmp = tempfile::TempDir::new().context("failed to create temporary directory")?;
    let err = tmp.path().join("amrfinder_index.err");
    let command = build_amrfinder_index_command(invocation, &err);
    exec(&command, &err.to_string_lossy()).map_err(anyhow::Error::msg)
}

fn build_amrfinder_index_command(invocation: &AmrFinderIndexInvocation, err: &Path) -> String {
    let mut command = format!(
        "{} {}",
        shell_quote(&invocation.program.to_string_lossy()),
        shell_quote(&dir_name_path(&invocation.version_dir).to_string_lossy())
    );
    if let Some(blast_bin) = &invocation.blast_bin {
        command.push_str(" -blast_bin ");
        command.push_str(&shell_quote(&dir_name_path(blast_bin).to_string_lossy()));
    }
    if let Some(hmmer_bin) = &invocation.hmmer_bin {
        command.push_str(" -hmmer_bin ");
        command.push_str(&shell_quote(&dir_name_path(hmmer_bin).to_string_lossy()));
    }
    if invocation.quiet {
        command.push_str(" -q");
    }
    if invocation.debug {
        command.push_str(" --debug");
    }
    command.push_str(" > ");
    command.push_str(&shell_quote(&err.to_string_lossy()));
    command
}

fn dir_name_path(path: &Path) -> PathBuf {
    let path = path.display().to_string();
    assert!(!path.is_empty());

    let mut items = Vec::new();
    let mut rest = path.as_str();
    while !rest.is_empty() {
        if let Some(pos) = rest.find('/') {
            items.push(rest[..pos].to_string());
            rest = &rest[pos + 1..];
        } else {
            items.push(rest.to_string());
            rest = "";
        }
    }

    let mut i = 0;
    while i < items.len() {
        if items[i].is_empty() && i != 0 {
            items.remove(i);
        } else {
            i += 1;
        }
    }

    let mut i = 0;
    while i < items.len() {
        if items[i] == "." {
            items.remove(i);
        } else {
            i += 1;
        }
    }

    let mut i = 0;
    while i < items.len() {
        if items[i] == ".." && i != 0 {
            items.remove(i);
            let prev = i - 1;
            if items[prev].is_empty() {
                assert_eq!(prev, 0);
                items[prev] = "..".to_string();
                i = prev;
            } else {
                items.remove(prev);
                i = prev;
            }
        } else {
            i += 1;
        }
    }

    let path = if items.is_empty() {
        ".".to_string()
    } else {
        let joined = items.join("/");
        if joined.is_empty() {
            "/".to_string()
        } else {
            joined
        }
    };
    if path.ends_with('/') {
        PathBuf::from(path)
    } else {
        PathBuf::from(format!("{path}/"))
    }
}

fn join_url(base: &str, parts: &[&str]) -> String {
    let mut url = base.trim_end_matches('/').to_string();
    for part in parts {
        url.push('/');
        url.push_str(part.trim_matches('/'));
    }
    url
}

fn dir_url(url: &str) -> String {
    if url.ends_with('/') {
        url.to_string()
    } else {
        format!("{url}/")
    }
}

fn normalize_version(text: &str) -> Result<String> {
    let Some(version) = text.lines().next().map(str::trim) else {
        bail!("empty version");
    };
    Ok(version.to_string())
}

fn parse_latest_minor(listing: &str) -> Result<String> {
    let mut best: Option<SoftwareVersion> = None;
    for token in version_tokens(listing) {
        if let Ok(version) = SoftwareVersion::from_text(&token, true) {
            if best.is_none_or(|current| version > current) {
                best = Some(version);
            }
        }
    }
    best.map(|version| version.get_minor())
        .context("no major.minor directory found")
}

fn parse_latest_data_version(listing: &str) -> Result<String> {
    let mut best: Option<DataVersion> = None;
    for token in data_version_tokens(listing) {
        if let Ok(version) = DataVersion::from_text(&token) {
            if best.is_none_or(|current| version > current) {
                best = Some(version);
            }
        }
    }
    best.map(|version| version.to_string())
        .context("no data version directory found")
}

fn version_tokens(listing: &str) -> Vec<String> {
    anchor_directory_tokens(listing)
}

fn data_version_tokens(listing: &str) -> Vec<String> {
    anchor_directory_tokens(listing)
}

fn anchor_directory_tokens(listing: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    for line in listing.lines() {
        let line = line.trim();
        if !line.starts_with("<a href=") {
            continue;
        }
        let Some(pos1) = line.find('>') else {
            continue;
        };
        let body = &line[pos1 + 1..];
        let Some(pos2) = body.find("/<") else {
            continue;
        };
        tokens.push(body[..pos2].to_string());
    }
    tokens
}

fn dna_point_mutation_taxgroups(taxgroup_path: &Path) -> Result<Vec<String>> {
    let text = fs::read_to_string(taxgroup_path)
        .with_context(|| format!("failed to read {}", taxgroup_path.display()))?;
    let mut taxgroups = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }

        let mut fields = line.split_whitespace();
        let taxgroup = fields
            .next()
            .with_context(|| format!("Bad {}", taxgroup_path.display()))?;
        let _gpipe = fields.next();
        let mutation_count_s = fields
            .next()
            .with_context(|| format!("Bad {}", taxgroup_path.display()))?;
        let end = mutation_count_s
            .char_indices()
            .take_while(|(idx, c)| c.is_ascii_digit() || (*idx == 0 && matches!(*c, '+' | '-')))
            .map(|(idx, c)| idx + c.len_utf8())
            .last()
            .unwrap_or(0);
        let mutation_count = mutation_count_s[..end]
            .parse::<i64>()
            .with_context(|| format!("Bad {}", taxgroup_path.display()))?;
        if mutation_count < 0 {
            bail!("Bad {}", taxgroup_path.display());
        }
        if mutation_count > 0 {
            taxgroups.push(taxgroup.to_string());
        }
    }
    Ok(taxgroups)
}

fn create_latest_link(main_dir: &Path, data_version: &str) -> Result<()> {
    let latest = main_dir.join(LATEST_DIR);
    if let Ok(metadata) = latest.symlink_metadata() {
        if metadata.is_dir() && !metadata.file_type().is_symlink() {
            let _ = fs::remove_dir(&latest);
        } else {
            let _ = fs::remove_file(&latest);
        }
    }
    if Path::new(data_version).is_absolute() {
        bail!(
            "Cannot make a relative symlink for {:?} as {:?} because {:?} is absolute",
            data_version,
            latest.display().to_string(),
            data_version
        );
    }
    let target = main_dir.join(data_version);
    if target.symlink_metadata().is_err() {
        bail!(
            "Cannot make a relative symlink for {:?} as {:?} because {:?} does not exist",
            data_version,
            latest.display().to_string(),
            target.display().to_string()
        );
    }

    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(data_version, &latest)
            .with_context(|| format!("failed to create {}", latest.display()))?;
    }

    #[cfg(not(unix))]
    {
        fs::write(&latest, data_version)
            .with_context(|| format!("failed to write {}", latest.display()))?;
    }

    Ok(())
}
