use std::path::PathBuf;
use std::process;

use amrfinder::common::Application;

fn main() {
    let mut app = Application {
        description: "Split FASTA file into approximately equal parts",
        usage: "Usage: fasta2parts <in> <parts_max> <dir>",
        positionals: 3,
        needs_arg: false,
        threads_used: false,
        key_args: Vec::new(),
    };
    let run = match app.run(std::env::args_os().collect()) {
        Ok(run) => run,
        Err(code) => process::exit(code),
    };

    let fasta = PathBuf::from(&run.positional_args[0]);
    let parts_max = match run.positional_args[1].to_string_lossy().parse::<usize>() {
        Ok(parts_max) => parts_max,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(1);
        }
    };
    let dir = PathBuf::from(&run.positional_args[2]);

    if let Err(e) = amrfinder::fasta2parts::body(&fasta, parts_max, &dir) {
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
