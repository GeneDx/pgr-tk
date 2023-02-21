const VERSION_STRING: &str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, CommandFactory, Parser};

use pgr_bin::{pair_shmmrs, sequence_to_shmmrs, SeqIndexDB, ShmmrSpec};
use rustc_hash::FxHashSet;
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
    /// option to process data from pre-build AGC backed index
    agc_idx_prefix: Option<String>,
    /// option to process data from pre-build frg backed index
    frg_idx_prefix: Option<String>,
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

fn generate_bed_graph_from_fastx_files(args: &CmdOptions) {

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
        File::create(Path::new(&args.prefix).with_extension("0.bedgraph"))
            .expect("can't create the output file"),
    );

    let mut output_file1 = BufWriter::new(
        File::create(Path::new(&args.prefix).with_extension("1.bedgraph"))
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
            assert!(c0 > 0);
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

fn generate_bed_graph_from_sdb(args: &CmdOptions, input_type: &str) {
    let mut seq_index_db = SeqIndexDB::new();
    if input_type == "AGC" {
        let _ = seq_index_db.load_from_agc_index(args.agc_idx_prefix.as_ref().unwrap().clone());
    } else if input_type == "FRG" {
        let _ = seq_index_db.load_from_frg_index(args.frg_idx_prefix.as_ref().unwrap().clone());
    } else {
        panic!("input type has to be specified  AGC or FRG backends")
    };

    let shmmr_spec = ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: false,
    };

    let mut sample_set0 = FxHashSet::<String>::default();
    let input_files = BufReader::new(
        File::open(Path::new(&args.filepath0))
            .expect("can't open the input file that contains the paths to the fastx files"),
    );
    input_files
        .lines()
        .into_iter()
        .for_each(|filename| {
            let filepath = filename
                .expect("can't get fastx file name")
                .trim()
                .to_string();
            sample_set0.insert(filepath);
        });

    let mut sample_set1 = FxHashSet::<String>::default();
    let input_files = BufReader::new(
        File::open(Path::new(&args.filepath1))
            .expect("can't open the input file that contains the paths to the fastx files"),
    );
    input_files
        .lines()
        .into_iter()
        .for_each(|filename| {
            let filepath = filename
                .expect("can't get fastx file name")
                .trim()
                .to_string();
            sample_set1.insert(filepath);
        });


    let mut sample_id_set0 = FxHashSet::<u32>::default();
    let mut sample_id_set1 = FxHashSet::<u32>::default();
    seq_index_db.seq_info.as_ref().unwrap().iter().for_each(|(sid, (_ctg, src, _))|{
        let src = src.clone().unwrap_or(String::new());
        if sample_set0.contains(&src) {
            sample_id_set0.insert(*sid);
        }
        if sample_set1.contains(&src) {
            sample_id_set1.insert(*sid);
        }
    });

    let frag_map = seq_index_db.get_shmmr_map_internal().unwrap();
    let mut output_file0 = BufWriter::new(
        File::create(Path::new(&args.prefix).with_extension("0.bedgraph"))
            .expect("can't create the output file"),
    );

    let mut output_file1 = BufWriter::new(
        File::create(Path::new(&args.prefix).with_extension("1.bedgraph"))
            .expect("can't create the output file"),
    );
    sample_id_set0.iter().for_each(|sid| {
        let seq_info = seq_index_db.seq_info.as_ref().unwrap().get(sid).unwrap(); 
        let ctg = seq_info.0.clone();
        let seq = seq_index_db.get_seq_by_id(*sid).unwrap();
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
            let (c0, c1) = if let Some(hits) = frag_map.get(&(k.0, k.1)) {
                let mut c0 = 0_u32;
                let mut c1 = 0_u32;
                hits.iter().for_each(|v| {
                    if sample_id_set0.contains(&v.1) {
                        c0 += 1;
                    }

                    if sample_id_set1.contains(&v.1) {
                        c1 += 1;
                    }
                });
                (c0, c1)
            } else {
                (0, 0)
            };
            assert!(c0 > 0);
            let r = c0 as f32 / c1 as f32;
            let _ = writeln!(
                output_file0,
                "{}\t{}\t{}\t{}\t{}\t{}",
                ctg, k.2, k.3, r, c1, c0
            );
        });
    });

    sample_id_set1.iter().for_each(|sid| {
        let seq_info = seq_index_db.seq_info.as_ref().unwrap().get(sid).unwrap(); 
        let ctg = seq_info.0.clone();
        let seq = seq_index_db.get_seq_by_id(*sid).unwrap();
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
            let (c0, c1) = if let Some(hits) = frag_map.get(&(k.0, k.1)) {
                let mut c0 = 0_u32;
                let mut c1 = 0_u32;
                hits.iter().for_each(|v| {
                    if sample_id_set0.contains(&v.1) {
                        c0 += 1;
                    }

                    if sample_id_set1.contains(&v.1) {
                        c1 += 1;
                    }
                });
                (c0, c1)
            } else {
                (0, 0)
            };
            assert!(c1 > 0);
            let r = c1 as f32 / c0 as f32;
            let _ = writeln!(
                output_file1,
                "{}\t{}\t{}\t{}\t{}\t{}",
                ctg, k.2, k.3, r, c1, c0
            );
        });
    });
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    if let Some(_agc_idx_prefix) = Some(args.agc_idx_prefix.clone()) {
        generate_bed_graph_from_sdb(&args, "AGC");
    } else if let Some(_frg_idx_prefix) = Some(args.frg_idx_prefix.clone()) {
        generate_bed_graph_from_sdb(&args, "FRG");
    } else {
        generate_bed_graph_from_fastx_files(&args);
    }
}
