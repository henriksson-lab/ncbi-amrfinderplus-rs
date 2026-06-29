use std::io;
use std::path::PathBuf;
use std::process;

use amrfinder::common::{Application, Key};

fn main() {
    let mut app = Application {
        description: "Extract sequences out of a FASTA file",
        usage: "Usage: fasta_extract <fasta> <target> [-aa]",
        positionals: 2,
        needs_arg: false,
        threads_used: false,
        key_args: vec![Key {
            name: "aa",
            flag: true,
            default_value: "false",
            acronym: None,
        }],
    };
    let run = match app.run(std::env::args_os().collect()) {
        Ok(run) => run,
        Err(code) => process::exit(code),
    };

    let fasta = PathBuf::from(&run.positional_args[0]);
    let target = PathBuf::from(&run.positional_args[1]);
    let stdout = io::stdout();
    let mut out = stdout.lock();

    if let Err(e) = amrfinder::fasta_extract::body(
        &fasta,
        &target,
        run.key_args["aa"] == "true",
        false,
        &mut out,
    ) {
        eprintln!("{}", e);
        process::exit(1);
    }
    if !run.key_args["log"].is_empty() {
        let _ = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&run.key_args["log"]);
    }
    if !run.key_args["json"].is_empty() {
        if let Err(e) = std::fs::write(&run.key_args["json"], "{}\n") {
            eprintln!("{}", e);
            process::exit(1);
        }
    }
}
