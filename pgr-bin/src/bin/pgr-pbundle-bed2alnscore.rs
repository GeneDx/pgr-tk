const VERSION_STRING: &'static str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader};
use std::{fs::File, path};
use svg::node::{self, element, Node};
use svg::Document;

#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-svg")]
#[clap(author, version)]
#[clap(about = "generate alignment scores between contigs using bundle decomposition from a principal bundle bed file", long_about = None)]
struct CmdOptions {
    bed_file_path: String,
    output_prefix: String,
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct BundleSegement {
    bgn: u32,
    end: u32,
    bundle_id: u32,
    bundle_v_count: u32,
    bundle_dir: u32,
    bundle_v_bgn: u32,
    bundle_v_end: u32,
}

enum AlnType {
    Match,
    MisMatch,
    Deletion,
    Insertion,
}

fn align_bundles(q_bundles: &Vec<BundleSegement>, t_bundles: &Vec<BundleSegement>) -> i64 {
    let q_count = q_bundles.len();
    let t_count = t_bundles.len();
    let mut S = FxHashMap::<(usize, usize), i64>::default();
    //let mut T = FxHashMap::<(usize, usize), AlnType>::default();

    let get_aln_direction_with_best_score =
        |q_idx: usize, t_idx: usize, S: &FxHashMap<(usize, usize), i64>| -> (AlnType, i64) {
            let mut best = (AlnType::Match, 0);
            let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
            let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
            let (max_len, min_len) = if q_len > t_len {
                (q_len, t_len)
            } else {
                (t_len, q_len)
            };
            best = match q_bundles[q_idx].bundle_id == t_bundles[t_idx].bundle_id {
                true => (
                    AlnType::Match,
                    2 * min_len + S.get(&(q_idx - 1, t_idx - 1)).unwrap(),
                ),
                false => (
                    AlnType::MisMatch,
                    -2 * max_len + S.get(&(q_idx - 1, t_idx - 1)).unwrap(),
                ),
            };
            let insert_score = -2 * q_len + S.get(&(q_idx, t_idx - 1)).unwrap();
            let delete_score = -2 * t_len + S.get(&(q_idx - 1, t_idx)).unwrap();
            if insert_score > best.1 {
                best = (AlnType::Insertion, insert_score)
            };
            if delete_score > best.1 {
                best = (AlnType::Deletion, delete_score)
            }
            best
        };

    let mut best_score = 0;

    (0..t_count)
        .flat_map(|t_idx| (0..q_count).map(move |q_idx| (q_idx, t_idx)))
        .for_each(|(q_idx, t_idx)| {
            //println!("{} {}", q_idx, t_idx);
            if q_idx == 0 || t_idx == 0 {
                S.insert((q_idx, t_idx), 0);
                return;
            };
            let (_, score) = get_aln_direction_with_best_score(q_idx, t_idx, &S);
            S.insert((q_idx, t_idx), score);
            if score > best_score {
                best_score = score
            }
        });
    best_score
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let bed_file_path = path::Path::new(&args.bed_file_path);
    let bed_file = BufReader::new(File::open(bed_file_path)?);
    let mut ctg_data = FxHashMap::<String, Vec<_>>::default();

    bed_file.lines().into_iter().for_each(|line| {
        let line = line.unwrap();
        let bed_fields = line.split("\t").collect::<Vec<&str>>();
        let ctg: String = bed_fields[0].to_string();
        let bgn: u32 = bed_fields[1].parse().unwrap();
        let end: u32 = bed_fields[2].parse().unwrap();
        let pbundle_fields = bed_fields[3].split(":").collect::<Vec<&str>>();
        let bundle_id: u32 = pbundle_fields[0].parse().unwrap();
        let bundle_v_count: u32 = pbundle_fields[1].parse().unwrap();
        let bundle_dir: u32 = pbundle_fields[2].parse().unwrap();
        let bundle_v_bgn: u32 = pbundle_fields[3].parse().unwrap();
        let bundle_v_end: u32 = pbundle_fields[4].parse().unwrap();

        let e = ctg_data.entry(ctg).or_insert(vec![]);
        let b_seg = BundleSegement {
            bgn,
            end,
            bundle_id,
            bundle_v_count,
            bundle_dir,
            bundle_v_bgn,
            bundle_v_end,
        };
        e.push(b_seg);
    });

    let mut ctg_data = ctg_data
        .into_iter()
        .map(|(k, mut v)| {
            v.sort();
            (k, v)
        })
        .collect::<Vec<_>>();

    ctg_data.sort();
    let n_ctg = ctg_data.len();

    (0..n_ctg)
        .flat_map(|ctg_idx0| (0..n_ctg).map(move |ctg_idx1| (ctg_idx0, ctg_idx1)))
        .for_each(|(ctg_idx0, ctg_idx1)| {
            let (ctg0, bundles0) = &ctg_data[ctg_idx0];
            let (ctg1, bundles1) = &ctg_data[ctg_idx1];
            let score = align_bundles(bundles0, bundles1);
            println!("{} {} {}", ctg0, ctg1, score);
        });
    Ok(())
}
