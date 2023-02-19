const VERSION_STRING: &str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, CommandFactory, Parser};

use pgr_bin::{pair_shmmrs, sequence_to_shmmrs, SeqIndexDB, ShmmrSpec};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

/// Compare SHIMMER pair count in two input sequence files 
#[derive(Parser, Debug)]
#[clap(name = "pgr-make-frgdb")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the path to the file contains the paths to the first set of sequence
    filepath0: String,
    /// the path to the file contains the paths to the seconde set of sequence
    filepath1: String,
    /// output prefix
    prefix: String,
    /// minimizer window size
    #[clap(long, short, default_value_t = 80)]
    w: u32,
    /// minimizer k-mer size
    #[clap(long, short, default_value_t = 56)]
    k: u32,
    /// sparse minimizer (shimmer) reduction factor
    #[clap(long, short, default_value_t = 4)]
    r: u32,
    /// min span for neighboring minimiers
    #[clap(long, short, default_value_t = 64)]
    min_span: u32,
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    // TODO: to log file
    //println!("read data from files in {:?}", args.filepath);
    //println!("output prefix {:?}", args.prefix);
    let shmmr_spec = ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: false,
    };
    let mut sdb0 = SeqIndexDB::new();
    let input_files = BufReader::new(
        File::open(Path::new(&args.filepath0))
            .expect("can't open the input file that contains the paths to the fastx files"),
    );
    input_files
        .lines()
        .into_iter()
        .enumerate()
        .for_each(|(fid, filename)| {
            let filepath = filename
                .expect("can't get fastx file name")
                .trim()
                .to_string();
            if fid == 0 {
                sdb0.load_from_fastx(filepath.clone(), args.w, args.k, args.r, args.min_span)
                    .unwrap_or_else(|_| panic!("fail to read the fastx file: {}", filepath));
            } else {
                sdb0.append_from_fastx(filepath.clone())
                    .unwrap_or_else(|_| panic!("fail to read the fastx file: {}", filepath));
            }
        });

    let mut sdb1 = SeqIndexDB::new();
    let input_files = BufReader::new(
        File::open(Path::new(&args.filepath1))
            .expect("can't open the input file that contains the paths to the fastx files"),
    );
    input_files
        .lines()
        .into_iter()
        .enumerate()
        .for_each(|(fid, filename)| {
            let filepath = filename
                .expect("can't get fastx file name")
                .trim()
                .to_string();
            if fid == 0 {
                sdb1.load_from_fastx(filepath.clone(), args.w, args.k, args.r, args.min_span)
                    .unwrap_or_else(|_| panic!("fail to read the fastx file: {}", filepath));
            } else {
                sdb1.append_from_fastx(filepath.clone())
                    .unwrap_or_else(|_| panic!("fail to read the fastx file: {}", filepath));
            }
        });

    let mut output_file0 = BufWriter::new(
        File::create(Path::new(&args.prefix).with_extension("0.bed"))
            .expect("can't create the output file"),
    );

    let mut output_file1 = BufWriter::new(
        File::create(Path::new(&args.prefix).with_extension("1.bed"))
            .expect("can't create the output file"),
    );

    let frag_map0 = &sdb0.seq_db.as_ref().unwrap().frag_map;
    let frag_map1 = &sdb1.seq_db.as_ref().unwrap().frag_map;
    sdb0.seq_info.as_ref().unwrap().iter().for_each(|(sid, v)| {
        let seq = sdb0.get_sub_seq_by_id(*sid, 0, v.2 as usize).unwrap();
        let shmmrs = sequence_to_shmmrs(*sid, &seq, &shmmr_spec, false);
        let smps = pair_shmmrs(&shmmrs);
        smps.iter().for_each(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.x >> 8;
            let s1 = s1.x >> 8;
            let k = if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            };
            let c0 = if let Some(v) = frag_map0.get(&(k.0, k.1)) {
                v.len()
            } else {
                0
            };
            let c1 = if let Some(v) = frag_map1.get(&(k.0, k.1)) {
                v.len()
            } else {
                0
            };
            assert!(c0> 0); 
            let ctg = v.0.clone();
            let r = c1 as f32 / c0 as f32;
            let _ = writeln!(
                output_file0,
                "{}\t{}\t{}\t{}\t{}\t{}",
                ctg, k.2, k.3, r, c0, c1
            );
        });
    });

    sdb1.seq_info.as_ref().unwrap().iter().for_each(|(sid, v)| {
        let seq = sdb1.get_sub_seq_by_id(*sid, 0, v.2 as usize).unwrap();
        let shmmrs = sequence_to_shmmrs(*sid, &seq, &shmmr_spec, false);
        let smps = pair_shmmrs(&shmmrs);
        smps.iter().for_each(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.x >> 8;
            let s1 = s1.x >> 8;
            let k = if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            };
            let c0 = if let Some(v) = frag_map0.get(&(k.0, k.1)) {
                v.len()
            } else {
                0
            };
            let c1 = if let Some(v) = frag_map1.get(&(k.0, k.1)) {
                v.len()
            } else {
                0
            };
            assert!(c1 > 0);
            let ctg = v.0.clone();
            let r = c0 as f32 / c1 as f32;
            let _ = writeln!(
                output_file1,
                "{}\t{}\t{}\t{}\t{}\t{}",
                ctg, k.2, k.3, r, c1, c0
            );
        });
    });
}
