const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use pgr_bin::SeqIndexDB;
use rustc_hash::{FxHashMap, FxHashSet};
//use std::fs::File;
use std::{
    fs::File,
    io::BufWriter,
    io::{BufRead, BufReader, Write},
    path::Path,
};

#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-decomp")]
#[clap(author, version)]
#[clap(about = "take a fasta file and output the principal bundle decomposition though MAP Graph", long_about = None)]
struct CmdOptions {
    fastx_path: String,
    output_prefix: String,
    #[clap(long, short, default_value = None)]
    include: Option<String>,
    #[clap(short, default_value_t = 48)]
    w: u32,
    #[clap(short, default_value_t = 56)]
    k: u32,
    #[clap(short, default_value_t = 4)]
    r: u32,
    #[clap(long, default_value_t = 12)]
    min_span: u32,
    #[clap(long, default_value_t = 0)]
    min_cov: usize,
    #[clap(long, default_value_t = 8)]
    min_branch_size: usize,
    #[clap(long, default_value_t = 2500)]
    bundle_length_cutoff: usize,
    #[clap(long, default_value_t = 10000)]
    bundle_merge_distance: usize,
}

#[allow(clippy::type_complexity)]
fn group_smps_by_principle_bundle_id(
    smps: &[((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)],
    bundle_length_cutoff: usize,
    bundle_merge_distance: usize,
) -> Vec<Vec<((u64, u64, u32, u32, u8), usize, u32, usize)>> {
    let mut pre_bundle_id: Option<usize> = None;
    let mut pre_direction: Option<u32> = None;
    let mut all_partitions = vec![];
    let mut new_partition = vec![];
    smps.iter().for_each(|&(smp, bundle_info)| {
        if bundle_info.is_none() {
            return;
        };
        let bundle_info = bundle_info.unwrap();
        let d = if smp.4 == bundle_info.1 { 0_u32 } else { 1_u32 };
        let bid = bundle_info.0;
        let bpos = bundle_info.2;
        if pre_bundle_id.is_none() {
            new_partition.clear();
            new_partition.push((smp, bid, d, bpos));
            pre_bundle_id = Some(bid);
            pre_direction = Some(d);
            return;
        };
        if bid != pre_bundle_id.unwrap() || d != pre_direction.unwrap() {
            let l = new_partition.len();
            if new_partition[l - 1].0 .3 as usize - new_partition[0].0 .2 as usize
                > bundle_length_cutoff
            {
                all_partitions.push(new_partition.clone());
                new_partition.clear();
            } else {
                new_partition.clear();
            };
            pre_bundle_id = Some(bid);
            pre_direction = Some(d);
        };
        new_partition.push((smp, bid, d, bpos));
    });
    let l = new_partition.len();
    if l > 0
        && new_partition[l - 1].0 .3 as usize - new_partition[0].0 .2 as usize
            > bundle_length_cutoff
    {
        all_partitions.push(new_partition);
    };

    let mut rtn_partitions = vec![];

    if all_partitions.is_empty() {
        return rtn_partitions;
    }
    let mut partition = all_partitions[0].clone();
    (1..all_partitions.len()).into_iter().for_each(|idx| {
        let p = all_partitions[idx].clone();
        let p_len = partition.len();
        let p_end = partition[p_len - 1].0 .3;
        let p_bid = partition[p_len - 1].1;
        let p_d = partition[p_len - 1].2;
        let np_bgn = p[0].0 .2;
        let np_bid = p[0].1;
        let np_d = p[0].2;
        if p_bid == np_bid
            && p_d == np_d
            && (np_bgn as i64 - p_end as i64).abs() < bundle_merge_distance as i64
        {
            partition.extend(p);
        } else {
            rtn_partitions.push(partition.clone());
            partition = p;
        }
    });
    if !partition.is_empty() {
        rtn_partitions.push(partition);
    }
    rtn_partitions
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let cmd_string = std::env::args().collect::<Vec<String>>().join(" ");
    let mut seq_index_db = SeqIndexDB::new();
    let fastx_path = args.fastx_path.clone();
    seq_index_db
        .load_from_fastx(fastx_path.clone(), args.w, args.k, args.r, args.min_span)
        .unwrap_or_else(|_| panic!("can't read file {}", fastx_path));

    if args.include.is_some() {
        let f = BufReader::new(
            File::open(Path::new(&args.include.unwrap())).expect("can't opne the inlude file"),
        );
        let include_ctgs = f.lines().map(|c| c.unwrap()).collect::<FxHashSet<String>>();
        let seq_list = include_ctgs
            .into_iter()
            .map(|ctg| {
                let seq = seq_index_db
                    .get_seq(fastx_path.clone(), ctg.clone())
                    .expect("fail to fetch sequence");
                (ctg, seq)
            })
            .collect::<Vec<_>>();
        let mut new_seq_index_db = SeqIndexDB::new();
        let _ = new_seq_index_db.load_from_seq_list(
            seq_list,
            Some(fastx_path.as_str()),
            args.w,
            args.k,
            args.r,
            args.min_span,
        );
        seq_index_db = new_seq_index_db;
    };

    let (principal_bundles, sid_smps) =
        seq_index_db.get_principal_bundle_decomposition(args.min_cov, args.min_branch_size, None);

    let bid_to_size = principal_bundles
        .iter()
        .map(|v| (v.0, v.2.len()))
        .collect::<FxHashMap<usize, usize>>();
    let sid_smps: FxHashMap<u32, Vec<_>> = sid_smps.into_iter().collect();

    let output_prefix_path = Path::new(&args.output_prefix);
    seq_index_db.generate_mapg_gfa(
        0,
        output_prefix_path
            .with_extension("mapg.gfa")
            .to_str()
            .unwrap(),
        "from_fragmap",
        None,
    )?;

    seq_index_db.write_mapg_idx(
        output_prefix_path
            .with_extension("mapg.idx")
            .to_str()
            .unwrap(),
    )?;

    seq_index_db.generate_principal_mapg_gfa(
        args.min_cov,
        args.min_branch_size,
        output_prefix_path
            .with_extension("pmapg.gfa")
            .to_str()
            .unwrap(),
        None,
    )?;
    let mut outpu_bed_file =
        BufWriter::new(File::create(output_prefix_path.with_extension("bed"))?);

    writeln!(outpu_bed_file, "# cmd: {}", cmd_string).expect("bed file write error");

    let mut seq_info = seq_index_db
        .seq_info
        .unwrap()
        .into_iter()
        .map(|(k, v)| (k, v))
        .collect::<Vec<_>>();

    seq_info.sort_by_key(|k| k.1 .0.clone());

    seq_info.into_iter().for_each(|(sid, sdata)| {
        let (ctg, _src, _len) = sdata;
        let smps = sid_smps.get(&sid).unwrap();
        let smp_partitions = group_smps_by_principle_bundle_id(
            smps,
            args.bundle_length_cutoff,
            args.bundle_merge_distance,
        );
        smp_partitions.into_iter().for_each(|p| {
            let b = p[0].0 .2;
            let e = p[p.len() - 1].0 .3 + args.k;
            let bid = p[0].1;
            let direction = p[0].2;
            let _ = writeln!(
                outpu_bed_file,
                "{}\t{}\t{}\t{}:{}:{}:{}:{}",
                ctg,
                b,
                e,
                bid,
                bid_to_size[&bid],
                direction,
                p[0].3,
                p[p.len() - 1].3
            );
        });
    });
    Ok(())
}
