use std::process;

fn main() {
    match amrfinder::amrfinder_update::main(std::env::args_os().collect()) {
        Ok(code) => process::exit(code),
        Err(e) => {
            eprintln!("{}", e);
            process::exit(1);
        }
    }
}
