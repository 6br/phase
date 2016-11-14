#[macro_use]
extern crate log;
extern crate env_logger;
extern crate rust_htslib;
extern crate gcollections;
extern crate interval;

//use interval::Interval;
//use std::thread;
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

fn genotypes(record:&mut rust_htslib::bcf::Record) -> Option<rust_htslib::bcf::record::Genotypes<>> {
    return record.genotypes().ok();
}

fn main() {
    let input_file = "../chr1.recode.norm.vcf.gz";
    //let output_file = |x| ["../chr1_phased_haploid_",x,".vcf.gz"].join("");
    let output_file_0 = "../chr1_phased_haploid_0.vcf.gz";
    let output_file_1 = "../chr1_phased_haploid_1.vcf.gz";
    solve_chromosome(0, input_file, output_file_0);
    solve_chromosome(1, input_file, output_file_1);
}
fn solve_chromosome<'a>(chromosome: usize, input_file:&str, output_file:&str) {
    let bcf : bcf::Reader<> = bcf::Reader::new(&input_file).ok().expect("Error opening vcf.");
    let header = bcf::Header::with_template(&bcf.header);
    let mut out = bcf::Writer::new(&output_file, &header, true, true).ok().expect("Error opening vcf.");
    //let mut out = bcf::Writer::new(&"../chr1_phased_haploid_2.vcf", &header, true, true).ok().expect("Error opening vcf.");

    let mut heap = BinaryHeap::new();
    let mut index = 0;
    let mut previous_end = 0;

    for r in bcf.records() {
        let mut record = r.ok().expect("Error reading Vcf file.");
        out.translate(&mut record);
        out.subset(&mut record);
        //record.trim_alleles().ok().except("Error trimming alleles");

        debug! ("x={}", record.pos());
        let pos = record.pos();
        let alt = record.alleles()[chromosome].len();
        let end = pos + alt as u32;
        let interval = interval::IntervalSet::new(pos + 1, end);
        let pq = {
            let pq_opt = record.format(b"PQ").integer();
            match pq_opt {
                Ok(a) => a[0][0],
                Err(_) => 256
            }
        };
        let genotype = genotypes(&mut record).unwrap().get(0);
        //if genotype.fmt() == "1|1" {
        if genotype[chromosome].index().unwrap() == 1 {
            heap.push(Vcf {idx: index,raw_data: record, pq: pq, range:interval});
        }
        if previous_end < pos {
            resolve_duplication(&mut heap, &mut out);
            heap.clear();
        }
        index += 1;
        previous_end = end;
        //if (genotype[1].index().unwrap() == 1){
        //    map_haploid_2.insert(Vcf {raw_data : &record, pq: pq}, pos);
        //}
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
            if intervalset.overlap(&range) {
                info!("#Remove SNPs at {}, quality: {}", range, pq);
            } else {
                debug!("{}, {}", pq, range);
                intervalset = intervalset.union(&range);
                output_heap.push(Vcf{pq:pq, idx:idx,raw_data:raw_data,range:range});
            }
        }
        output_from_heap(&mut output_heap, out)
    }
}

fn output_from_heap(output_heap : &mut BinaryHeap<Vcf>, out: &mut Writer){
    while let Some(a) = output_heap.pop() {
        out.write(&a.raw_data).ok().expect("Error writing Vcf file.");
    }
}
