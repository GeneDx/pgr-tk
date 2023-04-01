const VERSION_STRING: &str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, CommandFactory, Parser};

use pgr_bin::{pair_shmmrs, sequence_to_shmmrs, SeqIndexDB};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

/// Compare SHIMMER pair count in two input sequence files
#[derive(Parser, Debug)]
#[clap(name = "pgr-compare-cov2")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// option to process data from pre-build AGC backed index
    #[clap(long, short)]
    agc_idx_prefix: Option<String>,
    /// option to process data from pre-build frg backed index
    #[clap(long, short)]
    frg_idx_prefix: Option<String>,
    /// the path to the file contains the paths to the first set of sequence
    input: String,
    /// coverage threshold
    #[clap(long, short, default_value_t = 2.0)]
    threshold: f32,
}

fn filter_and_group_regions(
    regions: &Vec<(u32, u32, f32, usize, usize)>,
    max_dist: usize,
    min_range: usize,
) -> Vec<(u32, u32, f32, f32, f32)> {
    if regions.is_empty() {
        return vec![];
    }

    let mut chunk = Vec::<_>::new();
    let mut chunks = Vec::<_>::new();
    regions.iter().for_each(|&v| {
        if chunk.is_empty() {
            chunk.push(v);
            return;
        }

        if ((v.0 - chunk[chunk.len() - 1].1) as usize) < max_dist {
            chunk.push(v);
        } else {
            if ((chunk[chunk.len() - 1].1 - chunk[0].0) as usize) > min_range {
                chunks.push(chunk.clone());
            }
            chunk.clear();
        }
    });

    if !chunk.is_empty() && ((chunk[chunk.len() - 1].1 - chunk[0].0) as usize) > min_range {
        chunks.push(chunk);
    }
    chunks
        .into_iter()
        .map(|v| {
            let mut s0 = 0_f32;
            let mut s1 = 0_f32;
            let mut s2 = 0_f32;
            let bgn = v[0].0;
            let end = v[v.len() - 1].1;
            let len = v.len() as f32;
            v.into_iter().for_each(|vv| {
                s0 += vv.2 as f32;
                s1 += vv.3 as f32;
                s2 += vv.4 as f32;
            });
            (bgn, end, s0 / len, s1 / len, s2 / len)
        })
        .collect::<Vec<(u32, u32, f32, f32, f32)>>()
}

fn output_cov_bed(
    out_data: &Vec<(u32, u32, f32, usize, usize)>,
    ctg: String,
    prefix: String,
    threshold: f32,
    output_bed_file0: &mut BufWriter<File>,
) {
    let cov_high = out_data
        .iter()
        .filter(|&v| v.2 > threshold + 0.0001)
        .map(|v| v.clone())
        .collect::<Vec<_>>();

    let cov_high = filter_and_group_regions(&cov_high, 10000, 10000);

    let cov_low = out_data
        .iter()
        .filter(|&v| v.2 < threshold - 0.0001)
        .map(|v| v.clone())
        .collect::<Vec<_>>();

    let cov_low = filter_and_group_regions(&cov_low, 100, 20000);

    let mut cov = Vec::<_>::new();
    cov.extend_from_slice(&cov_high);
    cov.extend_from_slice(&cov_low);

    cov.sort_by(|a, b| a.0.cmp(&b.0));

    cov.into_iter().for_each(|v| {
        let _ = writeln!(
            output_bed_file0,
            "{}\t{}\t{}\t{}:{}\t{}\t{}",
            ctg, v.0, v.1, prefix, v.2, v.3, v.4
        );
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

    let shmmr_spec = seq_index_db.shmmr_spec.clone().unwrap();

    let input_files = BufReader::new(
        File::open(Path::new(&args.input))
            .expect("can't open the input file that contains the paths to the fastx files"),
    );
    input_files
        .lines()
        .into_iter()
        .collect::<Vec<_>>()
        .into_par_iter()
        .for_each(|filename| {
            let mut sample_set0 = FxHashSet::<String>::default();
            let mut sample_set1 = FxHashSet::<String>::default();
            let files = filename.expect("can't get fastx file name");
            let files = files.trim();
            let files = files.to_string();
            let mut files = files.split("\t");
            let prefix = files.next().unwrap().to_string();
            let refrenece = files.next().unwrap().to_string();
            println!("reference: {}", refrenece);
            sample_set0.insert(refrenece.clone());
            files.into_iter().for_each(|s| {
                println!("sample: {}", s);
                sample_set1.insert(s.to_string());
            });
            let mut sample_id_set0 = FxHashSet::<u32>::default();
            let mut sample_id_set1 = FxHashSet::<u32>::default();
            seq_index_db
                .seq_info
                .as_ref()
                .unwrap()
                .iter()
                .for_each(|(sid, (_ctg, src, _))| {
                    let src = src.clone().unwrap_or(String::new());
                    if sample_set0.contains(&src) {
                        sample_id_set0.insert(*sid);
                    }
                    if sample_set1.contains(&src) {
                        sample_id_set1.insert(*sid);
                    }
                });

            let frag_map = seq_index_db.get_shmmr_map_internal().unwrap();
            let mut output_bedgraph_file0 = BufWriter::new(
                File::create(Path::new(&prefix).with_extension("0.bedgraph"))
                    .expect("can't create the output file"),
            );

            let mut output_bed_file0 = BufWriter::new(
                File::create(Path::new(&prefix).with_extension("0.bed"))
                    .expect("can't create the output file"),
            );

            let mut output_bedgraph_file1 = BufWriter::new(
                File::create(Path::new(&prefix).with_extension("1.bedgraph"))
                    .expect("can't create the output file"),
            );

            let mut output_bed_file1 = BufWriter::new(
                File::create(Path::new(&prefix).with_extension("1.bed"))
                    .expect("can't create the output file"),
            );
            sample_id_set0.iter().for_each(|sid| {
                let seq_info = seq_index_db.seq_info.as_ref().unwrap().get(sid).unwrap();
                let ctg = seq_info.0.clone();
                let seq = seq_index_db.get_seq_by_id(*sid).unwrap();
                let shmmrs = sequence_to_shmmrs(*sid, &seq, &shmmr_spec, false);
                let smps = pair_shmmrs(&shmmrs);
                let out_data = smps
                    .iter()
                    .map(|(s0, s1)| {
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
                            let mut c0 = 0_usize;
                            let mut c1 = 0_usize;
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
                        let r = c1 as f32 / c0 as f32;
                        (k.2, k.3, r, c0, c1)
                    })
                    .collect::<Vec<_>>();
                output_cov_bed(
                    &out_data,
                    ctg.clone(),
                    prefix.clone(),
                    args.threshold,
                    &mut output_bed_file0,
                );

                out_data.into_iter().for_each(|v| {
                    let _ = writeln!(
                        output_bedgraph_file0,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        ctg, v.0, v.1, v.2, v.3, v.4
                    );
                });
            });

            sample_id_set1.iter().for_each(|sid| {
                let seq_info = seq_index_db.seq_info.as_ref().unwrap().get(sid).unwrap();
                let ctg = seq_info.0.clone();
                let seq = seq_index_db.get_seq_by_id(*sid).unwrap();
                let shmmrs = sequence_to_shmmrs(*sid, &seq, &shmmr_spec, false);
                let smps = pair_shmmrs(&shmmrs);
                let out_data = smps
                    .iter()
                    .map(|(s0, s1)| {
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
                            let mut c0 = 0_usize;
                            let mut c1 = 0_usize;
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
                        let r = c0 as f32 / c1 as f32;
                        (k.2, k.3, r, c1, c0)
                    })
                    .collect::<Vec<_>>();
                output_cov_bed(
                    &out_data,
                    ctg.clone(),
                    prefix.clone(),
                    1.0 / args.threshold,
                    &mut output_bed_file1,
                );

                out_data.into_iter().for_each(|v| {
                    let _ = writeln!(
                        output_bedgraph_file1,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        ctg, v.0, v.1, v.2, v.3, v.4
                    );
                });
            });
        });
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    if let Some(_agc_idx_prefix) = args.agc_idx_prefix.clone() {
        generate_bed_graph_from_sdb(&args, "AGC");
    } else if let Some(_frg_idx_prefix) = args.frg_idx_prefix.clone() {
        generate_bed_graph_from_sdb(&args, "FRG");
    } else {
        panic!("need a AGC/FRC backed seq index and db")
    }
}
