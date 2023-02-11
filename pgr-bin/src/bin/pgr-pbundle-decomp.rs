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

/// Generat the principal bundle decomposition though MAP Graph from a fasta file
#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-decomp")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the path to the input fasta file for building the principle bundles
    fastx_path: String,
    /// the prefix of the output files
    output_prefix: String,
    #[clap(long, short, default_value = None)]
    /// the path to the fasta file for principal bundle decomposition. if not specified, using the same one from from "<FASTX_PATH>"
    decomp_fastx_path: Option<String>,
    /// the path to the file that contains a list of contig name to be analyzed
    #[clap(long, short, default_value = None)]
    include: Option<String>,
    /// the SHIMMER parameter w
    #[clap(short, default_value_t = 48)]
    w: u32,
    /// the SHIMMER parameter k
    #[clap(short, default_value_t = 56)]
    k: u32,
    /// the SHIMMER parameter r
    #[clap(short, default_value_t = 4)]
    r: u32,
    /// the SHIMMER parameter minimum span length
    #[clap(long, default_value_t = 12)]
    min_span: u32,
    /// minimum coverage to be included in principal bundles
    #[clap(long, default_value_t = 0)]
    min_cov: usize,
    /// the minimum branch length to be included in the principal bundles
    #[clap(long, default_value_t = 8)]
    min_branch_size: usize,
    /// the minimum local project bundle size to includes
    #[clap(long, default_value_t = 2500)]
    bundle_length_cutoff: usize,
    /// merge two bundles with the same id with the specified length
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

    let decomp_seq_index_db = if let Some(decomp_fastx_path) = args.decomp_fastx_path {
        let decomp_fastx_path = decomp_fastx_path.clone();
        let mut decomp_seq_index_db = SeqIndexDB::new();
        decomp_seq_index_db
            .load_from_fastx(decomp_fastx_path, args.w, args.k, args.r, args.min_span)
            .unwrap_or_else(|_| panic!("can't read file {}", fastx_path));
        decomp_seq_index_db
    } else {
        //The file is read using a Mmap which is not clonable, need to rebuild the database. TODO: fix this.
        let decomp_fastx_path = fastx_path.clone();
        let mut decomp_seq_index_db = SeqIndexDB::new();
        decomp_seq_index_db
            .load_from_fastx(decomp_fastx_path, args.w, args.k, args.r, args.min_span)
            .unwrap_or_else(|_| panic!("can't read file {}", fastx_path));
        decomp_seq_index_db
    };

    let (principal_bundles, sid_smps) = seq_index_db.get_principal_bundle_decomposition(
        args.min_cov,
        args.min_branch_size,
        Some(&decomp_seq_index_db),
        None,
    );

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

    let mut output_ctg_summary_file = BufWriter::new(File::create(
        output_prefix_path.with_extension("ctg.summary.tsv"),
    )?);

    writeln!(outpu_bed_file, "# cmd: {}", cmd_string).expect("bed file write error");

    let mut seq_info = decomp_seq_index_db
        .seq_info
        .unwrap()
        .into_iter()
        .map(|(k, v)| (k, v))
        .collect::<Vec<_>>();

    seq_info.sort_by_key(|k| k.1 .0.clone());

    let mut repeat_count = FxHashMap::<u32, Vec<u32>>::default();
    let mut non_repeat_count = FxHashMap::<u32, Vec<u32>>::default();

    seq_info.iter().for_each(|(sid, sdata)| {
        let (ctg, _src, _len) = sdata;
        let smps = sid_smps.get(&sid).unwrap();
        let smp_partitions = group_smps_by_principle_bundle_id(
            smps,
            args.bundle_length_cutoff,
            args.bundle_merge_distance,
        );
        let mut ctg_bundle_count = FxHashMap::<usize, usize>::default();
        smp_partitions.iter().for_each(|p| {
            let bid = p[0].1;
            *ctg_bundle_count.entry(bid).or_insert_with(|| 0) += 1;
        });
        smp_partitions.into_iter().for_each(|p| {
            let b = p[0].0 .2;
            let e = p[p.len() - 1].0 .3 + args.k;
            let bid = p[0].1;
            let direction = p[0].2;
            let is_repeat;
            if *ctg_bundle_count.get(&bid).unwrap_or(&0) > 1 {
                repeat_count
                    .entry(*sid)
                    .or_insert_with(|| vec![])
                    .push(e - b - args.k);
                is_repeat = "R";
            } else {
                non_repeat_count
                    .entry(*sid)
                    .or_insert_with(|| vec![])
                    .push(e - b - args.k);
                is_repeat = "U";
            }
            let _ = writeln!(
                outpu_bed_file,
                "{}\t{}\t{}\t{}:{}:{}:{}:{}:{}",
                ctg,
                b,
                e,
                bid,
                bid_to_size[&bid],
                direction,
                p[0].3,
                p[p.len() - 1].3,
                is_repeat
            );
        });
    });
    let _ = writeln!(
        output_ctg_summary_file,
        "#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        "ctg",
        "length",
        "repeat_bundle_count",
        "repeat_bundle_sum",
        "repeat_bunlde_percentage",
        "repeat_bundle_mean",
        "repeat_bundle_min",
        "repeat_bundle_max",
        "non_repeat_bundle_count",
        "non_repeat_bundle_sum",
        "non_repeat_bunlde_percentage",
        "non_repeat_bundle_mean",
        "non_repeat_bundle_min",
        "non_repeat_bundle_max",
        "total_bundle_count",
        "total_bundle_coverage_percentage"
    );
    seq_info.into_iter().for_each(|(sid, sdata)| {
        let (ctg, _src, len) = sdata;
        let repeat_bundle_count = repeat_count.get(&sid).unwrap_or(&vec![]).len();
        let non_repeat_bundle_count = non_repeat_count.get(&sid).unwrap_or(&vec![]).len();
        let repeat_sum = repeat_count
            .get(&sid)
            .unwrap_or(&vec![])
            .iter()
            .fold(0, |acc, x| acc + x);
        let repeat_bundle_max = repeat_count
            .get(&sid)
            .unwrap_or(&vec![])
            .into_iter()
            .fold(0, |x, &y| if x > y { x } else { y });
        let repeat_bundle_min = repeat_count
            .get(&sid)
            .unwrap_or(&vec![])
            .into_iter()
            .fold(len, |x, &y| if x < y { x } else { y });
        let non_repeat_sum = non_repeat_count
            .get(&sid)
            .unwrap_or(&vec![])
            .into_iter()
            .fold(0, |acc, x| acc + x);
        let non_repeat_bundle_max = non_repeat_count
            .get(&sid)
            .unwrap_or(&vec![])
            .into_iter()
            .fold(0, |x, &y| if x > y { x } else { y });
        let non_repeat_bundle_min = non_repeat_count
            .get(&sid)
            .unwrap_or(&vec![])
            .into_iter()
            .fold(len, |x, &y| if x < y { x } else { y });
        let repeat_bundle_mean = if repeat_bundle_count > 0 {
            format!("{}", repeat_sum as f32 / repeat_bundle_count as f32)
        } else {
            "NA".to_string()
        };
        let non_repeat_bundle_mean = if non_repeat_bundle_count > 0 {
            format!("{}", non_repeat_sum as f32 / non_repeat_bundle_count as f32)
        } else {
            "NA".to_string()
        };
        let repeat_bundle_min = if repeat_bundle_count > 0 {
            format!("{}", repeat_bundle_min)
        } else {
            "NA".to_string()
        };
        let repeat_bundle_max = if repeat_bundle_count > 0 {
            format!("{}", repeat_bundle_max)
        } else {
            "NA".to_string()
        };
        let non_repeat_bundle_min = if non_repeat_bundle_count > 0 {
            format!("{}", non_repeat_bundle_min)
        } else {
            "NA".to_string()
        };
        let non_repeat_bundle_max = if non_repeat_bundle_count > 0 {
            format!("{}", non_repeat_bundle_max)
        } else {
            "NA".to_string()
        };

        let _ = writeln!(
            output_ctg_summary_file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            ctg,
            len,
            repeat_bundle_count,
            repeat_sum,
            100.0 * repeat_sum as f32 / len as f32,
            repeat_bundle_mean,
            repeat_bundle_min,
            repeat_bundle_max,
            non_repeat_bundle_count,
            non_repeat_sum,
            100.0 * non_repeat_sum as f32 / len as f32,
            non_repeat_bundle_mean,
            non_repeat_bundle_min,
            non_repeat_bundle_max,
            repeat_bundle_count + non_repeat_bundle_count,
            100.0 * (repeat_sum + non_repeat_sum) as f32 / len as f32,
        );
    });
    Ok(())
}
