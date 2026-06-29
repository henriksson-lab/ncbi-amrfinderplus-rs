use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Once;
use std::time::SystemTime;

pub const DEFAULT_VERSION: &str = "0.0.0";

#[cfg(unix)]
const SIGPIPE_NUM: i32 = 13;
#[cfg(unix)]
const SIG_DFL_HANDLER: usize = 0;

#[cfg(unix)]
static INIT_COMMON: Once = Once::new();
#[cfg(unix)]
static SIGPIPE_EXIT_SUCCESS: AtomicBool = AtomicBool::new(false);

#[cfg(unix)]
unsafe extern "C" {
    fn signal(signum: i32, handler: usize) -> usize;
}

thread_local! {
    static PROGRAM_ARGS: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
    static PROGRAM_NAME: RefCell<String> = const { RefCell::new(String::new()) };
    static EXEC_DIR: RefCell<String> = const { RefCell::new(String::new()) };
    static PROG2DIR: RefCell<HashMap<String, String>> = RefCell::new(HashMap::new());
    static QC_ON: RefCell<bool> = const { RefCell::new(false) };
    static QUIET_ENABLED: RefCell<bool> = const { RefCell::new(false) };
    static PROGRESS_ENABLED: RefCell<bool> = const { RefCell::new(true) };
    static PROFILE_ENABLED: RefCell<bool> = const { RefCell::new(false) };
    static SIGPIPE: RefCell<bool> = const { RefCell::new(false) };
    static THREADS_MAX: RefCell<usize> = const { RefCell::new(1) };
    static REQUIRED_GROUPS: RefCell<HashMap<usize, HashMap<&'static str, String>>> = RefCell::new(HashMap::new());
}

#[cfg(unix)]
extern "C" fn sigpipe_handler(_sig_num: i32) {
    unsafe {
        signal(SIGPIPE_NUM, SIG_DFL_HANDLER);
    }
    std::process::exit(if SIGPIPE_EXIT_SUCCESS.load(Ordering::Relaxed) {
        0
    } else {
        1
    });
}

#[cfg(unix)]
pub fn init_common() -> bool {
    INIT_COMMON.call_once(|| unsafe {
        signal(SIGPIPE_NUM, sigpipe_handler as usize);
    });
    true
}

#[cfg(not(unix))]
pub fn init_common() -> bool {
    true
}

pub fn str_null(s: &str) -> bool {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return true;
    }
    matches!(trimmed.to_ascii_uppercase().as_str(), "NULL" | "NA" | "N/A")
}

pub fn get_scientific(number_s: &str) -> (bool, bool, usize) {
    let s = number_s.to_ascii_uppercase();
    let e_pos = s.find('E').unwrap_or(s.len());
    let point_pos = s[..e_pos].find('.');
    let has_point = point_pos.is_some();
    let decimals = point_pos.map_or(0, |point| e_pos - point - 1);
    (e_pos < s.len(), has_point, decimals)
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct SoftwareVersion {
    pub major: u32,
    pub minor: u32,
    pub patch: u32,
}

impl SoftwareVersion {
    pub fn new(major: u32, minor: u32, patch: u32) -> Self {
        Self {
            major,
            minor,
            patch,
        }
    }

    pub fn from_file(f_name: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(f_name).map_err(|_| {
            format!(
                "Cannot read software version. One line is expected in the file: '{}'",
                f_name.display()
            )
        })?;
        let mut lines = text.lines();
        let Some(line) = lines.next() else {
            return Err(format!(
                "Cannot read software version. One line is expected in the file: '{}'",
                f_name.display()
            ));
        };
        if lines.next().is_some() {
            return Err(format!(
                "Cannot read software version. One line is expected in the file: '{}'",
                f_name.display()
            ));
        }
        Self::from_text(line, false)
    }

    pub fn from_text(s: &str, minor_only: bool) -> Result<Self, String> {
        let s = s
            .split_whitespace()
            .next()
            .ok_or_else(|| "Cannot read software version".to_string())?;
        let mut parts = s.split('.');
        let major = parts
            .next()
            .ok_or_else(|| "Cannot read software version".to_string())?
            .parse::<u32>()
            .map_err(|_| "Cannot read software version".to_string())?;
        if minor_only {
            let minor = parts
                .next()
                .ok_or_else(|| "Cannot read software version".to_string())?
                .parse::<u32>()
                .map_err(|_| "Cannot read software version".to_string())?;
            if parts.next().is_some() {
                return Err("Cannot read software version".to_string());
            }
            Ok(Self::new(major, minor, 0))
        } else {
            let minor = parts
                .next()
                .ok_or_else(|| "Cannot read software version".to_string())?
                .parse::<u32>()
                .map_err(|_| "Cannot read software version".to_string())?;
            let patch = parts
                .next()
                .ok_or_else(|| "Cannot read software version".to_string())?
                .parse::<u32>()
                .map_err(|_| "Cannot read software version".to_string())?;
            if parts.next().is_some() {
                return Err("Cannot read software version".to_string());
            }
            Ok(Self::new(major, minor, patch))
        }
    }

    pub fn save_text(&self, os: &mut dyn Write) -> std::io::Result<()> {
        write!(os, "{}.{}.{}", self.major, self.minor, self.patch)
    }

    pub fn get_minor(&self) -> String {
        format!("{}.{}", self.major, self.minor)
    }
}

impl std::fmt::Display for SoftwareVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}.{}", self.major, self.minor, self.patch)
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct DataVersion {
    pub year: u32,
    pub month: u32,
    pub day: u32,
    pub num: u32,
}

impl DataVersion {
    pub fn new(year: u32, month: u32, day: u32, num: u32) -> Self {
        Self {
            year,
            month,
            day,
            num,
        }
    }

    pub fn from_file(f_name: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(f_name).map_err(|_| {
            format!(
                "Cannot read data version. One line is expected in the file: '{}'",
                f_name.display()
            )
        })?;
        let mut lines = text.lines();
        let Some(line) = lines.next() else {
            return Err(format!(
                "Cannot read data version. One line is expected in the file: '{}'",
                f_name.display()
            ));
        };
        if lines.next().is_some() {
            return Err(format!(
                "Cannot read data version. One line is expected in the file: '{}'",
                f_name.display()
            ));
        }
        Self::from_text(line)
    }

    pub fn from_text(s: &str) -> Result<Self, String> {
        let s = s
            .split_whitespace()
            .next()
            .ok_or_else(|| "Cannot read data version".to_string())?;
        let (year, rest) = s
            .split_once('-')
            .ok_or_else(|| "Cannot read data version".to_string())?;
        let (month, rest) = rest
            .split_once('-')
            .ok_or_else(|| "Cannot read data version".to_string())?;
        let (day, num) = rest
            .split_once('.')
            .ok_or_else(|| "Cannot read data version".to_string())?;
        Ok(Self::new(
            year.parse()
                .map_err(|_| "Cannot read data version".to_string())?,
            month
                .parse()
                .map_err(|_| "Cannot read data version".to_string())?,
            day.parse()
                .map_err(|_| "Cannot read data version".to_string())?,
            num.parse()
                .map_err(|_| "Cannot read data version".to_string())?,
        ))
    }

    pub fn save_text(&self, os: &mut dyn Write) -> std::io::Result<()> {
        write!(
            os,
            "{}-{:02}-{:02}.{}",
            self.year, self.month, self.day, self.num
        )
    }
}

impl std::fmt::Display for DataVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}-{:02}-{:02}.{}",
            self.year, self.month, self.day, self.num
        )
    }
}

pub struct Key {
    pub name: &'static str,
    pub flag: bool,
    pub default_value: &'static str,
    pub acronym: Option<char>,
}

pub struct Application {
    pub description: &'static str,
    pub usage: &'static str,
    pub positionals: usize,
    pub needs_arg: bool,
    pub threads_used: bool,
    pub key_args: Vec<Key>,
}

pub struct Run {
    pub positional_args: Vec<OsString>,
    pub key_args: HashMap<&'static str, String>,
}

pub fn get_command_line() -> String {
    PROGRAM_ARGS.with(|program_args| {
        program_args
            .borrow()
            .iter()
            .map(|arg| {
                if arg.is_empty()
                    || arg.chars().any(|c| {
                        matches!(
                            c,
                            ' ' | '|'
                                | ';'
                                | '#'
                                | '*'
                                | '?'
                                | '$'
                                | '('
                                | ')'
                                | '<'
                                | '>'
                                | '~'
                                | '\''
                                | '"'
                                | '\\'
                        )
                    })
                {
                    shell_quote(arg)
                } else {
                    arg.clone()
                }
            })
            .collect::<Vec<_>>()
            .join(" ")
    })
}

pub fn shell_quote(s: &str) -> String {
    format!("'{}'", s.replace('\'', "'\\''"))
}

fn un_quote(s: &str) -> String {
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        s[1..s.len() - 1].replace("'\\''", "'")
    } else {
        s.to_string()
    }
}

pub fn which(prog_name: &str) -> String {
    if prog_name.is_empty() {
        return String::new();
    }
    let Some(path) = std::env::var_os("PATH") else {
        return String::new();
    };
    for dir in std::env::split_paths(&path) {
        if dir.as_os_str().is_empty() {
            continue;
        }
        if dir.join(prog_name).exists() {
            let mut s = dir.to_string_lossy().into_owned();
            if !s.ends_with(std::path::MAIN_SEPARATOR) {
                s.push(std::path::MAIN_SEPARATOR);
            }
            return s;
        }
    }
    String::new()
}

pub fn qc_on() -> bool {
    QC_ON.with(|qc_on| *qc_on.borrow())
}

pub fn quiet_enabled() -> bool {
    QUIET_ENABLED.with(|quiet_enabled| *quiet_enabled.borrow())
}

pub fn progress_enabled() -> bool {
    PROGRESS_ENABLED.with(|progress_enabled| *progress_enabled.borrow())
}

pub fn profile_enabled() -> bool {
    PROFILE_ENABLED.with(|profile_enabled| *profile_enabled.borrow())
}

pub fn sigpipe_enabled() -> bool {
    SIGPIPE.with(|sigpipe| *sigpipe.borrow())
}

pub fn get_threads_max() -> usize {
    THREADS_MAX.with(|threads_max| *threads_max.borrow())
}

pub fn get_now() -> String {
    let output = Command::new("date")
        .arg("+%Y-%m-%d %H:%M:%S")
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .map(|s| s.trim_end().to_string());
    output.unwrap_or_else(|| "0000-00-00 00:00:00".to_string())
}

pub fn make_temp_dir() -> Result<PathBuf, String> {
    let tmp_dir = std::env::var_os("TMPDIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp"));
    let program_name = PROGRAM_NAME.with(|program_name| program_name.borrow().clone());
    let program_name = if program_name.is_empty() {
        "program"
    } else {
        &program_name
    };
    let temp_dir = tempfile::Builder::new()
        .prefix(&format!("{program_name}."))
        .tempdir_in(&tmp_dir)
        .map_err(|_| {
            format!(
                "Error creating a temporary directory in {}",
                tmp_dir.display()
            )
        })?;
    let temp_path = temp_dir.keep();
    let test_path = temp_path.join("test");
    fs::write(&test_path, "abc\n").map_err(|_| {
        format!(
            "{} is full, make space there or use environment variable TMPDIR to change location for temporary files",
            tmp_dir.display()
        )
    })?;
    fs::remove_file(&test_path).map_err(|e| e.to_string())?;
    Ok(temp_path)
}

pub fn remove_directory(dir_name: &Path) -> Result<(), String> {
    fs::remove_dir_all(dir_name).map_err(|e| format!("Cannot remove directory {:?}: {e}", dir_name))
}

pub fn exec(cmd: &str, log_fname: &str) -> Result<(), String> {
    if cmd.is_empty() {
        return Err("empty command".to_string());
    }
    let status = Command::new("sh")
        .arg("-c")
        .arg(cmd)
        .status()
        .map_err(|e| e.to_string())?;
    if status.success() {
        return Ok(());
    }

    let mut err = format!("{cmd}\nstatus = {status}");
    if !log_fname.is_empty() {
        if let Ok(log) = fs::read_to_string(log_fname) {
            let tail = log.lines().take(10).collect::<Vec<_>>().join("\n");
            if !tail.is_empty() {
                err.push('\n');
                err.push_str(&tail);
            }
        }
    }
    Err(err)
}

impl Application {
    pub fn set_required_group(
        &mut self,
        key_name: &'static str,
        required_group: &str,
    ) -> Result<(), String> {
        if required_group.is_empty() {
            return Err("Empty required group".to_string());
        }
        if !self.key_args.iter().any(|key| key.name == key_name) {
            return Err(format!("Unknown key: \"{}\"", key_name));
        }
        REQUIRED_GROUPS.with(|required_groups| {
            required_groups
                .borrow_mut()
                .entry(self as *const Self as usize)
                .or_default()
                .insert(key_name, required_group.to_string());
        });
        Ok(())
    }

    fn qc_args(&self) -> Result<(), String> {
        let mut names = HashSet::new();
        for key in &self.key_args {
            let name = key.name.trim();
            if name.is_empty() {
                return Err("Empty key argument name".to_string());
            }
            if key.name != name || name.contains(' ') || name.contains('=') || name.starts_with('-')
            {
                return Err(format!("Bad key argument name \"{}\"", key.name));
            }
            if name == "help" || name == "version" {
                return Err(format!("Reserved key argument name \"{}\"", key.name));
            }
            if !names.insert(name) {
                return Err(format!("Duplicate key argument name \"{}\"", key.name));
            }
        }
        let mut acronyms = HashSet::new();
        for key in &self.key_args {
            if let Some(acronym) = key.acronym {
                if !acronyms.insert(acronym) {
                    return Err(format!("Duplicate key argument acronym \"{}\"", acronym));
                }
            }
        }
        Ok(())
    }

    pub fn add_default_args(&mut self) {
        const DEFAULT_ARGS: &[Key] = &[
            Key {
                name: "qc",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "verbose",
                flag: false,
                default_value: "0",
                acronym: None,
            },
            Key {
                name: "noprogress",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "profile",
                flag: true,
                default_value: "false",
                acronym: None,
            },
            Key {
                name: "seed",
                flag: false,
                default_value: "1",
                acronym: None,
            },
            Key {
                name: "json",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "log",
                flag: false,
                default_value: "",
                acronym: None,
            },
            Key {
                name: "sigpipe",
                flag: true,
                default_value: "false",
                acronym: None,
            },
        ];
        if self.threads_used && !self.key_args.iter().any(|key| key.name == "threads") {
            self.key_args.push(Key {
                name: "threads",
                flag: false,
                default_value: "1",
                acronym: None,
            });
        }
        for key in DEFAULT_ARGS {
            if !self
                .key_args
                .iter()
                .any(|existing| existing.name == key.name)
            {
                self.key_args.push(Key {
                    name: key.name,
                    flag: key.flag,
                    default_value: key.default_value,
                    acronym: key.acronym,
                });
            }
        }
    }

    pub fn run(&mut self, argv: Vec<OsString>) -> Result<Run, i32> {
        init_common();
        self.add_default_args();
        if let Err(e) = self.qc_args() {
            eprintln!("{e}");
            return Err(1);
        }

        let mut program_args: Vec<OsString> = Vec::new();
        for (i, arg) in argv.into_iter().enumerate() {
            let key = arg.to_string_lossy();
            if i > 0 && key.starts_with('-') && !key.is_empty() {
                if let Some((key, value)) = key.split_once('=') {
                    program_args.push(OsString::from(key));
                    program_args.push(OsString::from(value));
                    continue;
                }
            }
            program_args.push(arg);
        }
        PROGRAM_ARGS.with(|stored_args| {
            *stored_args.borrow_mut() = program_args
                .iter()
                .map(|arg| arg.to_string_lossy().into_owned())
                .collect();
        });
        if let Some(program) = program_args.first() {
            let name = Path::new(program)
                .file_name()
                .unwrap_or(program.as_os_str())
                .to_string_lossy()
                .into_owned();
            PROGRAM_NAME.with(|program_name| *program_name.borrow_mut() = name);
        }
        let program_args_string = program_args
            .iter()
            .map(|arg| arg.to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        let exec_dir = program_args
            .first()
            .and_then(|program| {
                let path = Path::new(program);
                path.canonicalize()
                    .ok()
                    .or_else(|| Some(path.to_path_buf()))
                    .and_then(|path| path.parent().map(Path::to_path_buf))
            })
            .and_then(|path| path.to_str().map(str::to_string));
        if let Some(exec_dir) = &exec_dir {
            let mut dir = exec_dir.clone();
            if !dir.ends_with(std::path::MAIN_SEPARATOR) {
                dir.push(std::path::MAIN_SEPARATOR);
            }
            EXEC_DIR.with(|stored_exec_dir| *stored_exec_dir.borrow_mut() = dir);
        }

        let mut key_args: HashMap<&'static str, String> = HashMap::new();
        let mut keys_read: HashSet<&'static str> = HashSet::new();
        let mut positional_args: Vec<OsString> = Vec::new();
        let mut value_key: Option<&Key> = None;
        for (i, s) in program_args.iter().enumerate() {
            if i == 0 {
                continue;
            }
            if let Some(key) = value_key {
                key_args.insert(key.name, s.to_string_lossy().into_owned());
                value_key = None;
            } else {
                let s_lossy = s.to_string_lossy();
                if !s_lossy.is_empty() && s_lossy.starts_with('-') {
                    let short = !s_lossy.starts_with("--");
                    let name = if let Some(name) = s_lossy.strip_prefix("--") {
                        name
                    } else {
                        &s_lossy[1..]
                    };
                    if name.is_empty() {
                        if s_lossy.starts_with("--") {
                            eprintln!("Dashes with no key\n\n{}", self.usage);
                        } else {
                            eprintln!("Dash with no key\n\n{}", self.usage);
                        }
                        return Err(1);
                    }
                    if name == "help" || (short && name == "h") {
                        if program_args.len() > 2 {
                            eprintln!("Single keyword \"help\" is expected\n\n{}", self.usage);
                            return Err(1);
                        }
                        println!("{}", self.usage);
                        return Err(0);
                    }
                    if name == "version" || (short && name == "v") {
                        if program_args.len() > 2 {
                            eprintln!("Single keyword \"version\" is expected\n\n{}", self.usage);
                            return Err(1);
                        }
                        println!("{}", DEFAULT_VERSION);
                        return Err(0);
                    }
                    let key = if short && name.chars().count() == 1 {
                        let acronym = name.chars().next().unwrap();
                        self.key_args
                            .iter()
                            .find(|key| key.acronym == Some(acronym))
                            .or_else(|| self.key_args.iter().find(|key| key.name == name))
                    } else {
                        self.key_args.iter().find(|key| key.name == name)
                    };
                    if let Some(key) = key {
                        if !keys_read.insert(key.name) {
                            eprintln!("Parameter \"{}\" is used more than once", key.name);
                            return Err(1);
                        }
                        if key.flag {
                            key_args.insert(key.name, "true".to_string());
                        } else {
                            value_key = Some(key);
                        }
                    } else if positional_args.len() == self.positionals {
                        eprintln!("\"{}\" is not a valid option\n\n{}", s_lossy, self.usage);
                        return Err(1);
                    } else {
                        positional_args.push(s.clone());
                    }
                } else if positional_args.len() == self.positionals {
                    eprintln!(
                        "\"{}\" cannot be a positional parameter\n\n{}",
                        s_lossy, self.usage
                    );
                    return Err(1);
                } else {
                    positional_args.push(s.clone());
                }
            }
        }
        if let Some(key) = value_key {
            eprintln!("Key with no value: {}\n\n{}", key.name, self.usage);
            return Err(1);
        }
        if program_args.len() == 1 && (self.positionals > 0 || self.needs_arg) {
            println!("{}", self.usage);
            return Err(1);
        }
        if positional_args.len() != self.positionals {
            eprintln!("Too few positional parameters\n\n{}", self.usage);
            return Err(1);
        }
        for key in &self.key_args {
            if !keys_read.contains(key.name) {
                let mut default_value = key.default_value.to_string();
                if !key.flag {
                    if let Some(exec_dir) = &exec_dir {
                        default_value = default_value.replace("$BASE", exec_dir);
                    }
                }
                key_args.insert(key.name, default_value);
            }
        }
        let required_groups = REQUIRED_GROUPS.with(|required_groups| {
            required_groups
                .borrow()
                .get(&(self as *const Self as usize))
                .cloned()
                .unwrap_or_default()
        });
        if !required_groups.is_empty() {
            let mut required_group2names: HashMap<String, Vec<String>> = HashMap::new();
            let mut required_group2given: HashMap<String, Vec<String>> = HashMap::new();
            for key in &self.key_args {
                let Some(required_group) = required_groups.get(key.name) else {
                    continue;
                };
                let value_given = if key.flag {
                    keys_read.contains(key.name)
                } else {
                    key_args
                        .get(key.name)
                        .is_some_and(|value| !value.is_empty())
                };
                if value_given {
                    let given = required_group2given
                        .entry(required_group.clone())
                        .or_default();
                    if !given.is_empty() {
                        eprintln!(
                            "Parameter \"{}\" is conflicting with {}",
                            key.name,
                            given.join(", ")
                        );
                        return Err(1);
                    }
                    given.push(format!("\"{}\"", key.name));
                }
                required_group2names
                    .entry(required_group.clone())
                    .or_default()
                    .push(format!("\"{}\"", key.name));
            }
            for (required_group, names) in required_group2names {
                if required_group2given
                    .get(&required_group)
                    .is_none_or(Vec::is_empty)
                {
                    eprintln!(
                        "One of required parameters is missing: {}",
                        names.join(", ")
                    );
                    return Err(1);
                }
            }
        }
        if let Some(seed) = key_args.get("seed") {
            match seed.parse::<u64>() {
                Ok(0) => {
                    eprintln!("Seed cannot be 0");
                    return Err(1);
                }
                Ok(_) => {}
                Err(_) => {
                    eprintln!("Cannot convert \"{}\" to number", seed);
                    return Err(1);
                }
            }
        }
        if let Some(threads) = key_args.get("threads") {
            match threads.parse::<usize>() {
                Ok(0) => {
                    eprintln!("Number of threads cannot be 0");
                    return Err(1);
                }
                Ok(threads) => {
                    if threads > 1
                        && key_args
                            .get("profile")
                            .is_some_and(|profile| profile == "true")
                    {
                        eprintln!("Cannot profile with threads");
                        return Err(1);
                    }
                    THREADS_MAX.with(|threads_max| *threads_max.borrow_mut() = threads);
                }
                Err(_) => {
                    eprintln!("Cannot convert \"{}\" to number", threads);
                    return Err(1);
                }
            }
        } else {
            THREADS_MAX.with(|threads_max| *threads_max.borrow_mut() = 1);
        }
        QC_ON.with(|qc_on| {
            *qc_on.borrow_mut() = key_args
                .get("qc")
                .or_else(|| key_args.get("debug"))
                .is_some_and(|value| value == "true")
        });
        QUIET_ENABLED.with(|quiet_enabled| {
            *quiet_enabled.borrow_mut() = key_args.get("quiet").is_some_and(|value| value == "true")
        });
        PROGRESS_ENABLED.with(|progress_enabled| {
            *progress_enabled.borrow_mut() = !key_args
                .get("noprogress")
                .is_some_and(|value| value == "true")
        });
        PROFILE_ENABLED.with(|profile_enabled| {
            *profile_enabled.borrow_mut() =
                key_args.get("profile").is_some_and(|value| value == "true")
        });
        SIGPIPE.with(|sigpipe| {
            let enabled = key_args.get("sigpipe").is_some_and(|value| value == "true");
            *sigpipe.borrow_mut() = enabled;
            #[cfg(unix)]
            SIGPIPE_EXIT_SUCCESS.store(enabled, Ordering::Relaxed);
        });
        if let Some(log_fname) = key_args.get("log") {
            if !log_fname.is_empty() {
                match OpenOptions::new().create(true).append(true).open(log_fname) {
                    Ok(mut log) => {
                        let command_line = program_args_string
                            .iter()
                            .map(|arg| {
                                if arg.is_empty()
                                    || arg.chars().any(|c| {
                                        matches!(
                                            c,
                                            ' ' | '|'
                                                | ';'
                                                | '#'
                                                | '*'
                                                | '?'
                                                | '$'
                                                | '('
                                                | ')'
                                                | '<'
                                                | '>'
                                                | '~'
                                                | '\''
                                                | '"'
                                                | '\\'
                                        )
                                    })
                                {
                                    shell_quote(arg)
                                } else {
                                    arg.clone()
                                }
                            })
                            .collect::<Vec<_>>()
                            .join(" ");
                        if let Err(e) = writeln!(log, "\n{}\n$ {}", get_now(), command_line) {
                            eprintln!("{e}");
                            return Err(1);
                        }
                    }
                    Err(e) => {
                        eprintln!("{e}");
                        return Err(1);
                    }
                }
            }
        }
        Ok(Run {
            positional_args,
            key_args,
        })
    }
}

pub struct ShellApplication {
    pub application: Application,
    pub use_tmp: bool,
    pub tmp: Option<PathBuf>,
    pub stderr_quiet: bool,
    pub start_time: Option<SystemTime>,
}

impl ShellApplication {
    pub fn init_var(&mut self) -> Result<(), String> {
        if self.tmp.is_some() {
            return Ok(());
        }
        if self.use_tmp {
            self.tmp = Some(make_temp_dir()?);
        }
        self.stderr_quiet = quiet_enabled();
        self.start_time = Some(SystemTime::now());
        if !self.stderr_quiet {
            eprintln!("Running: {}", get_command_line());
        }
        Ok(())
    }

    pub fn body<F>(&mut self, shell_body: F) -> Result<(), String>
    where
        F: FnOnce(&Self) -> Result<(), String>,
    {
        if let Some(tmp) = &self.tmp {
            if !quiet_enabled()
                && self
                    .application
                    .key_args
                    .iter()
                    .any(|key| key.name == "log" && !key.default_value.is_empty())
            {
                eprintln!("{}", tmp.display());
            }
        }
        shell_body(self)
    }

    pub fn find_prog(&self, prog_name: &str) -> Result<(), String> {
        if prog_name.is_empty() {
            return Err("Empty program name".to_string());
        }
        if prog_name.contains(std::path::MAIN_SEPARATOR) {
            return Err(format!("Program name contains a slash: {prog_name}"));
        }

        let cached = PROG2DIR.with(|prog2dir| prog2dir.borrow().contains_key(prog_name));
        if cached {
            return Ok(());
        }

        let exec_dir = EXEC_DIR.with(|exec_dir| exec_dir.borrow().clone());
        let dir = if !exec_dir.is_empty() && Path::new(&exec_dir).join(prog_name).exists() {
            exec_dir
        } else {
            which(prog_name)
        };
        if dir.is_empty() {
            let program_name = PROGRAM_NAME.with(|program_name| program_name.borrow().clone());
            return Err(format!(
                "Binary {} is not found.\nPlease make sure that {} is in the same directory as {} or is in your $PATH.",
                shell_quote(prog_name),
                shell_quote(prog_name),
                shell_quote(&program_name)
            ));
        }
        PROG2DIR.with(|prog2dir| {
            prog2dir.borrow_mut().insert(prog_name.to_string(), dir);
        });
        Ok(())
    }

    pub fn full_prog(&self, prog_name: &str) -> Result<String, String> {
        PROG2DIR.with(|prog2dir| {
            prog2dir
                .borrow()
                .get(prog_name)
                .map(|dir| format!("{} ", shell_quote(&format!("{dir}{prog_name}"))))
                .ok_or_else(|| format!("Program {} is not found", shell_quote(prog_name)))
        })
    }

    pub fn exec2str(&self, cmd: &str, tmp_name: &str, log_fname: &str) -> Result<String, String> {
        if tmp_name.contains(' ') {
            return Err(format!("Temporary file name contains a space: {tmp_name}"));
        }
        let tmp = self
            .tmp
            .as_ref()
            .ok_or_else(|| "Temporary directory is not initialized".to_string())?;
        let out = tmp.join(tmp_name);
        exec(
            &format!("{cmd} > {}", shell_quote(&out.to_string_lossy())),
            log_fname,
        )?;
        let text = fs::read_to_string(&out).map_err(|e| e.to_string())?;
        let lines = text.lines().collect::<Vec<_>>();
        if lines.len() != 1 {
            return Err(format!("{cmd}\nOne line is expected"));
        }
        Ok(lines[0].to_string())
    }

    pub fn uncompress(&self, quoted_fname: &str, suffix: &str) -> Result<String, String> {
        let tmp = self
            .tmp
            .as_ref()
            .ok_or_else(|| "Temporary directory is not initialized".to_string())?;
        let res = shell_quote(&tmp.join(suffix).to_string_lossy());
        if quoted_fname == res {
            return Err("Input and output file names are the same".to_string());
        }
        let s = un_quote(quoted_fname);
        if s.ends_with(".gz") {
            exec(&format!("gunzip -c {quoted_fname} > {res}"), "")?;
            return Ok(res);
        }
        Ok(quoted_fname.to_string())
    }

    pub fn get_blast_threads_param(
        &self,
        blast: &str,
        threads_max_max: usize,
    ) -> Result<String, String> {
        let tmp = self
            .tmp
            .as_ref()
            .ok_or_else(|| "Temporary directory is not initialized".to_string())?;
        let t = get_threads_max().min(threads_max_max);
        if t <= 1 {
            return Ok(String::new());
        }

        let blast_help = tmp.join("blast_help");
        exec(
            &format!(
                "{} -help > {}",
                self.full_prog(blast)?,
                shell_quote(&blast_help.to_string_lossy())
            ),
            "",
        )?;
        let help = fs::read_to_string(&blast_help).map_err(|e| e.to_string())?;
        let mut num_threads = false;
        let mut mt_mode = false;
        for line in help.lines() {
            let line = line.trim();
            if line.contains("-num_threads ") {
                num_threads = true;
            }
            if line.contains("-mt_mode ") {
                mt_mode = true;
            }
        }

        if !num_threads {
            return Ok(String::new());
        }
        let mut s = format!("  -num_threads {t}");
        #[cfg(not(target_os = "macos"))]
        let mt_mode_works = true;
        #[cfg(target_os = "macos")]
        let mut mt_mode_works = false;
        #[cfg(target_os = "macos")]
        {
            let blast_version = tmp.join("blast_version");
            exec(
                &format!(
                    "{} -version > {}",
                    self.full_prog(blast)?,
                    shell_quote(&blast_version.to_string_lossy())
                ),
                "",
            )?;
            let version = fs::read_to_string(&blast_version).map_err(|e| e.to_string())?;
            let prefix = format!("{blast}: ");
            if let Some(line) = version.lines().next() {
                let line = line.trim();
                if let Some(version) = line.strip_prefix(&prefix) {
                    let version = version.trim_end_matches('+');
                    let mut parts = version
                        .split('.')
                        .take(3)
                        .map(|part| part.parse::<u64>().unwrap_or(0));
                    let major = parts.next().unwrap_or(0);
                    let minor = parts.next().unwrap_or(0);
                    let patch = parts.next().unwrap_or(0);
                    mt_mode_works = (major, minor, patch) >= (2, 13, 0);
                }
            }
        }
        if mt_mode && mt_mode_works {
            s.push_str("  -mt_mode 1");
        }
        Ok(s)
    }
}

impl Drop for ShellApplication {
    fn drop(&mut self) {
        let keeps_tmp = self
            .application
            .key_args
            .iter()
            .any(|key| key.name == "log" && !key.default_value.is_empty());
        if !keeps_tmp {
            if let Some(tmp) = &self.tmp {
                let _ = remove_directory(tmp);
            }
        }
        if let Some(start_time) = self.start_time {
            if !self.stderr_quiet {
                let elapsed = start_time.elapsed().unwrap_or_default().as_secs();
                let program_name = PROGRAM_NAME.with(|program_name| program_name.borrow().clone());
                eprintln!("{program_name} took {elapsed} seconds to complete");
            }
        }
    }
}
