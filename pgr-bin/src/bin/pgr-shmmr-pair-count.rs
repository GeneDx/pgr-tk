const VERSION_STRING: &'static str = env!("VERSION_STRING");

//use std::path::PathBuf;
use clap::{self, IntoApp, Parser};

use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use pgr_db::seq_db;

#[derive(Parser, Debug)]
#[clap(name = "pgr-shmmr-pair-count")]
#[clap(author, version)]
#[clap(about = "count shimmer pairs in a shimmer database", long_about = None)]
struct CmdOptions {
    prefix: String,
    output_path: String,
    //max_unique_count
    #[clap(long, short, default_value_t = 1)]
    max_unique_count: usize,
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let (_shmmr_spec, shmmr_pair_to_frags) =
        seq_db::read_mdb_file(args.prefix.clone() + ".mdb").unwrap();
    let mut seq_index = HashMap::<(String, Option<String>), (u32, u32)>::new();
    let mut seq_info = HashMap::<u32, (String, Option<String>, u32)>::new();
    let midx_file = BufReader::new(File::open(args.prefix.clone() + ".midx")?);

    let mut sources = HashSet::<String>::new();
    midx_file
        .lines()
        .into_iter()
        .try_for_each(|line| -> Result<(), std::io::Error> {
            let line = line.unwrap();
            let mut line = line.as_str().split("\t");
            let sid = line.next().unwrap().parse::<u32>().unwrap();
            let len = line.next().unwrap().parse::<u32>().unwrap();
            let ctg_name = line.next().unwrap().to_string();
            let source = line.next().unwrap().to_string();
            sources.insert(source.clone());
            seq_index.insert((ctg_name.clone(), Some(source.clone())), (sid, len));
            seq_info.insert(sid, (ctg_name, Some(source), len));
            Ok(())
        })?;

    let source_to_id = sources
        .iter()
        .enumerate()
        .map(|v| (v.1.clone(), v.0 as u32))
        .collect::<HashMap<String, u32>>();

    let mut sid_to_source_id_lookup = vec![0_u32; seq_info.len()];

    seq_info.iter().for_each(|(k, v)| {
        sid_to_source_id_lookup[*k as usize] = *source_to_id.get(v.1.as_ref().unwrap()).unwrap();
    });
    let mut out_file = BufWriter::new(File::create(args.output_path)?);
    let out_vec = shmmr_pair_to_frags
        .par_iter()
        .map(|(k, v)| {
            let mut count = FxHashMap::<u32, usize>::default();
            v.iter().for_each(|v| {
                let sid = (*v).1;
                let source_id = *sid_to_source_id_lookup.get(sid as usize).unwrap();
                *count.entry(source_id).or_insert(0) += 1;
            });
            let v = count
                .into_iter()
                .filter(|(_k, v)| {
                    let muc = args.max_unique_count;
                    if *v > muc {
                        false
                    } else {
                        true
                    }
                })
                .count();
            (k.0, k.1, v)
        })
        .collect::<Vec<(u64, u64, usize)>>();

    out_vec
        .iter()
        .try_for_each(|v| -> Result<(), std::io::Error> {
            writeln!(&mut out_file, "{} {} {}", v.0, v.1, v.2)?;
            Ok(())
        })?;

    Ok(())
}
