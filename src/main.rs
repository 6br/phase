extern crate rust_htslib;
extern crate gcollections;
extern crate interval;

//use interval::Interval;
//use interval::ops::*;
//use gcollections::ops::*;

//use rust_htslib::htslib::vcf;
use rust_htslib::bcf;
use std::collections::BTreeMap;
use std::cmp::Ordering;
use rust_htslib::bcf::Record;
//use rust_htslib::htslib::Read;
//o#[derive(Debug)]u
#[derive(Debug)]
struct Vcf {
    //haploid: &Genotype,
    //pq: u32,
    //range: Interval,
    rawData: &rust_htslib::bcf::Record
}
/*
fn haploid(genotype:&GenotypeAllele) u32 {
    match genotype.index() {
        Some id => id
        
    } 
}
*/
fn main() {
    let bcf = bcf::Reader::new(&"../chr1.recode.vcf").ok().expect("Error opening vcf.");
    let header = bcf::Header::with_template(&bcf.header);
    let mut out = bcf::Writer::new(&"../chr1_phased_haploid_1.vcf", &header, true, true).ok().expect("Error opening vcf.");
    let mut out = bcf::Writer::new(&"../chr1_phased_haploid_2.vcf", &header, true, true).ok().expect("Error opening vcf.");

    let mut map_haploid_1 = BTreeMap::new();
    let mut map_haploid_2 = BTreeMap::new();

    for r in bcf.records() {
        //let record = r.ok().except("Error reading vcf file.");
        let record = r.ok().expect("Error reading BAM file.");
        println! ("x={}", record.pos())
        let pq = record.info("PQ").ok().expect("Error")
        println! ("x={}", pq)
        map_haploid_1.insert(Vcf {rawData = record}, pq)
    }
    println!("Hello, world!");

    threed::spawn(move || {
        
        //let result = do_something(data);
        //println!("{:?}", result); // => 2
    });
}
