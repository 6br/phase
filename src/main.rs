#![feature(box_syntax)]
#[macro_use] extern crate lazy_static;
#[macro_use] extern crate log;
extern crate env_logger;
extern crate rust_htslib;
extern crate gcollections;
extern crate interval;

extern crate getopts;
extern crate regex;

mod lib;

fn main() {
    env_logger::init().unwrap();
    let args = lib::option_parser();
    match args {
        Some (mut a) => match a.flag_thread {
            true => {lib::run()},
            false => {lib::run_sequencial(&mut a)}
        },
        None => return
    }
}

