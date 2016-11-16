#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rust_htslib;
extern crate gcollections;
extern crate interval;

extern crate getopts;
extern crate regex;

use getopts::Options;
use std::env;
use std::cmp::Ord;
use std::path::Path;
use regex::Regex;

//use interval::Interval;
use std::thread;
use interval::interval_set::ToIntervalSet;
use interval::ops::*;
use gcollections::ops::*;

//use rust_htslib::htslib::vcf;
//use std::collections::BTreeMap;
use rust_htslib::bcf;
use std::cmp::Ordering;
use rust_htslib::bcf::Writer;
use std::collections::BinaryHeap;

struct Vcf {
    idx: i32,
    pq: i32,
    range: interval::interval_set::IntervalSet<u32>,
    raw_data: rust_htslib::bcf::Record
}

struct VcfIndex {
    idx: i32,
    raw_data: rust_htslib::bcf::Record
}

impl PartialEq for Vcf {
    fn eq(&self, other: &Self) -> bool { self.idx == other.idx }
}

impl Eq for Vcf {}

impl PartialOrd for Vcf {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { self.pq.partial_cmp(&other.pq) }
}

impl Ord for Vcf {
    fn cmp(&self, other: &Self) -> Ordering { self.pq.cmp(&other.pq) }
}


impl PartialEq for VcfIndex {
    fn eq(&self, other: &Self) -> bool { self.idx == other.idx }
}

impl Eq for VcfIndex {}

impl PartialOrd for VcfIndex {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { self.idx.partial_cmp(&other.idx) }
}

impl Ord for VcfIndex {
    fn cmp(&self, other: &Self) -> Ordering { self.idx.cmp(&other.idx) }
}


fn genotypes(record:&mut rust_htslib::bcf::Record) -> Option<rust_htslib::bcf::record::Genotypes<>> {
    return record.genotypes().ok();
}


fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
}

pub struct Args {
    flag_thread: bool,
    cmd_output: bool,
    arg_output: String,
    arg_input: String,
    cmd_chr: bool,
    arg_chr: String,
    //flag_help: bool,
    //flag_version: bool
}

struct WriteArgs<'a> {
    path: &'a Path,
    uncompressed: bool,
    vcf: bool,
}

fn option_parser() -> Option<Args> {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "", "set output file name", "FILE");
    opts.optopt("c", "chr", "filiter by chromosome id", "CHRID");
    opts.optflag("t", "thread", "enable multi-thread");
    opts.optflag("h", "help", "print this help menu");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!(f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return None;
    }
    let output = matches.opt_str("o");
    let thread = matches.opt_str("t");
    let chr = matches.opt_str("c");
    let input = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, opts);
        return None;
    };
    return Some (Args {cmd_chr: chr.is_some() , flag_thread: thread.is_some(), cmd_output: output.is_some(), arg_input: input, arg_output: output.unwrap(), arg_chr: chr.unwrap()})
}

fn main() {
    env_logger::init().unwrap();
    // Input should be sorted by vcf-sort, and normalized by bcftools norm.
    let input_file = Path::new("../chr1.recode.norm.vcf.gz");
    if !input_file.exists() {
        println!("Error: No input file.");
        return;
    }
    let output_file_0 = Path::new("../chr1_phased_haploid_0.vcf");
    let output_file_1 = Path::new("../chr1_phased_haploid_1.vcf");

    let output_path = Path::new("../chr1_phased_haploid.vcf");
    let handle1 = thread::spawn(move || {
        //let ext = "_0.".to_string() + path.extension().unwrap().to_str().unwrap();
        solve_chromosome(0, input_file, output_file_0);//path.with_extension(ext).as_path());
    });
    let handle2 = thread::spawn(move || {
        //let ext = "_1.".to_string() + path.extension().unwrap().to_str().unwrap();
        solve_chromosome(1, input_file, output_file_1);//path.with_extension(ext).as_path());
    });
    let result = handle1.join();
    let result2 = handle2.join();
    assert!(!result.is_err());
    assert!(!result2.is_err());
}

fn gen_output_filename(output_file: &Path, suffix: String) -> WriteArgs {
    let extend_regex_vcfgz = Regex::new("(.*).vcf.gz").unwrap();
    let extend_regex_vcf = Regex::new("(.*).vcf").unwrap();
    let extend_regex_bcf = Regex::new("(.*).bcf").unwrap();

    return WriteArgs{path:output_file, vcf: true, uncompressed: true}
}

fn solve_chromosome<'a>(haploid: usize, input_file: &Path, output_file: &Path) {
    println!("#{} => {}", input_file.to_str().unwrap(), output_file.to_str().unwrap());
    let bcf : bcf::Reader<> = bcf::Reader::new(&input_file).ok().expect("Error opening vcf.");
    let mut header = bcf::Header::with_template(&bcf.header);
    let header_version = concat!("##vcf-phenotype-quality-filter=", env!("CARGO_PKG_VERSION"));
    let header2 = header.push_record(header_version.as_bytes());

    // Last 2 Argument means (uncompressed: bool, vcf: bool)
    let mut out = bcf::Writer::new(&output_file, header2, true, true).ok().expect("Error opening vcf.");

    let mut heap = BinaryHeap::new();
    let mut index = 0;
    let mut previous_end = 0;
    //let mut chr = "";

    for r in bcf.records() {
        let mut record = r.ok().expect("Error reading Vcf file.");
        out.translate(&mut record);
        out.subset(&mut record);
        //record.trim_alleles().ok().except("Error trimming alleles");

        debug! ("x={}", record.pos());
        let pos = record.pos();
        //let alt = record.alleles()[chromosome].len();
        let alt = record.inner().rlen;
        //println!("x={}, {}", record.pos(), record.inner().rlen);
        let end = pos + alt as u32;
        let interval = interval::IntervalSet::new(pos + 1, end);

        //println!("{} {}",pos, previous_end);
        if previous_end < pos {
            //println!("{} resdup",pos);
            resolve_duplication(&mut heap, &mut out);
            heap.clear();
            //if pos > 811368 { break; }
        }

        let pq = {
            let pq_opt = record.format(b"PQ").integer();
            match pq_opt {
                Ok(a) => a[0][0],
                Err(_) => 0
            }
        };
        let genotype = genotypes(&mut record).unwrap().get(0);
        //if genotype.fmt() == "1|1" {
        if genotype[haploid].index().unwrap() == 1 {
            heap.push(Vcf {idx: index,raw_data: record, pq: pq, range:interval});
        }
        index += 1;
        previous_end = end;
    }
    resolve_duplication(&mut heap, &mut out);
}

fn resolve_duplication(mut heap : &mut BinaryHeap<Vcf>, out: &mut Writer) {
    if heap.len() == 1 {
        output_from_heap(&mut heap, out)
    } else {
        let mut output_heap = BinaryHeap::new();
        let mut intervalset = vec![].to_interval_set();//interval::IntervalSet::new();

        while let Some(Vcf { pq, idx, raw_data, range }) = heap.pop() {
            //println!("{}, {}", pq, range);
            //println!("{}", intervalset.overlap(&range.upper()));
            if range.upper() - range.lower() == 0 &&  intervalset.overlap(&range.upper()) {
                println!("#Remove SNPs at {}, quality: {}", range, pq);
            } else if intervalset.overlap(&range) {
                println!("#Remove SNPs at {}, quality: {}", range, pq);
            } else {
                debug!("{}, {}", pq, range);
                intervalset = intervalset.union(&range);
                output_heap.push(VcfIndex{idx:idx, raw_data:raw_data});
            }
        }
        output_from_heap(&mut output_heap, out)
    }
}

fn output_from_heap<T : std::cmp::Ord>(output_heap : &mut BinaryHeap<T>, out: &mut Writer){
    while let Some(a) = output_heap.pop() {
        out.write(&a.raw_data).ok().expect("Error Writing Vcf file.");
    }
}
