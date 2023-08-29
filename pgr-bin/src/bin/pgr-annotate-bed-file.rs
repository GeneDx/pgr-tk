const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use flate2::bufread::MultiGzDecoder;
use iset::IntervalMap;
//use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Align long contigs and identify potential SV regions with respect to the reference fasta file
#[derive(Parser, Debug)]
#[clap(name = "pgr-annotate-bed-file")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the the a bed file
    bed_path: String,
    /// path to the annotation file (gzipped)
    annotation_path: String,
    /// the prefix of the output files
    output_path: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}
fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let mut reader = BufReader::new(File::open(Path::new(&args.annotation_path)).unwrap());

    let annotation_reader = BufReader::new(MultiGzDecoder::new(&mut reader));
    let mut annotation_interval = FxHashMap::<String, IntervalMap<u32, (char, String)>>::default();
    // we support https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz for now
    annotation_reader.lines().for_each(|line| {
        if let Ok(line) = line {
            let err_msg = format!("faile to parse on {}", line);
            let fields = line.split('\t').collect::<Vec<&str>>();
            let chr = fields[0].to_string();
            let f_type = fields[2].to_string();
            let fs = fields[3].parse::<u32>().expect(&err_msg);
            let fe = fields[4].parse::<u32>().expect(&err_msg);
            let strand = fields[6].chars().next().expect(&err_msg);
            let attribute = fields[8].to_string();
            if f_type == "transcript" {
                let e = annotation_interval
                    .entry(chr)
                    .or_insert(IntervalMap::<u32, (char, String)>::default());
                e.insert(fs..fe, (strand, attribute));
            }
        }
    });

    let mut out_bed = BufWriter::new(File::create(Path::new(&args.output_path)).unwrap());

    let bed_reader = BufReader::new(File::open(Path::new(&args.bed_path)).unwrap());
    bed_reader.lines().for_each(|line| {
        if let Ok(line) = line {
            if line.starts_with('#') {
                return;
            };
            let err_msg = format!("faile to parse on {}", line);
            let fields = line.split('\t').collect::<Vec<&str>>();
            let chr = fields[0].to_string();
            let bgn = fields[1].parse::<u32>().expect(&err_msg);
            let end = fields[2].parse::<u32>().expect(&err_msg);
            let annotation = fields[3].to_string();
            if let Some(i_map) = annotation_interval.get(&chr) {
                // TODO, we only pick the first overlap for now
                let mut annotations = FxHashSet::<String>::default();
                for (_strand, attributes) in i_map.values(bgn..end) {
                    // TODO: need a proper parser
                    let attributes = attributes.trim_end_matches(';').to_string();
                    let a_fields = attributes.split(';').collect::<Vec<&str>>();
                    let gn = a_fields.last().unwrap().to_string();
                    let gn = gn.split(' ').collect::<Vec<&str>>();
                    let gn = gn.last().unwrap().to_string();
                    let gn = gn.trim_matches('"');
                    annotations.insert(gn.to_string()); 
                };
                if annotations.is_empty() { return };
                let gn = annotations.into_iter().collect::<Vec<_>>().join("/"); 

                
                writeln!(
                    out_bed,
                    "{}\t{}\t{}\t{}>{}",
                    chr, bgn, end, annotation, gn
                )
                .expect("fail to write the vcf file");
            };
        }
    });

    Ok(())
}
