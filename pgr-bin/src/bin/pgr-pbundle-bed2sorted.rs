const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::{fs::File, path};

/// Generate annotation file with a sorting order from the principal bundle decomposition
#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-bed2sorted")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the path to the pricipal bundle bed file
    bed_file_path: String,
    /// the prefix of the output file
    output_prefix: String,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
struct BundleSegement {
    bgn: u32,
    end: u32,
    bundle_id: u32,
    bundle_v_count: u32,
    bundle_dir: u32,
    bundle_v_bgn: u32,
    bundle_v_end: u32,
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let bed_file_path = path::Path::new(&args.bed_file_path);
    let bed_file = BufReader::new(File::open(bed_file_path)?);
    let mut ctg_data = FxHashMap::<String, Vec<_>>::default();
    let bed_file_parse_err_msg = "bed file parsing error";
    let mut node_length = FxHashMap::<(u32, u32), Vec<_>>::default();
    bed_file.lines().into_iter().for_each(|line| {
        let line = line.unwrap().trim().to_string();
        if line.is_empty() {
            return;
        }
        if &line[0..1] == "#" {
            return;
        }
        let bed_fields = line.split('\t').collect::<Vec<&str>>();
        let ctg: String = bed_fields[0].to_string();
        let bgn: u32 = bed_fields[1].parse().expect(bed_file_parse_err_msg);
        let end: u32 = bed_fields[2].parse().expect(bed_file_parse_err_msg);
        let pbundle_fields = bed_fields[3].split(':').collect::<Vec<&str>>();
        let bundle_id: u32 = pbundle_fields[0].parse().expect(bed_file_parse_err_msg);
        let bundle_v_count: u32 = pbundle_fields[1].parse().expect(bed_file_parse_err_msg);
        let bundle_dir: u32 = pbundle_fields[2].parse().expect(bed_file_parse_err_msg);
        let bundle_v_bgn: u32 = pbundle_fields[3].parse().expect(bed_file_parse_err_msg);
        let bundle_v_end: u32 = pbundle_fields[4].parse().expect(bed_file_parse_err_msg);

        let e = ctg_data.entry(ctg).or_default();
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
        if (bundle_v_bgn as i64 - bundle_v_end as i64).abs() as f32 > (bundle_v_count as f32) * 0.5
        {
            let e = node_length.entry((bundle_id, bundle_dir)).or_default();
            e.push((end as i64 - bgn as i64).unsigned_abs() as u64);
        }
    });

    let mut node_length = node_length
        .into_iter()
        .map(|(n, v)| {
            let c = v.len() as f64;
            let sum = v.into_iter().sum::<u64>() as f64;
            (sum / c, n)
        })
        .collect::<Vec<(f64, (u32, u32))>>();
    node_length.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let mut ctg_data = ctg_data
        .into_iter()
        .map(|(ctg, mut bundle_segs)| {
            bundle_segs.sort();
            let mut node_count = FxHashMap::<(u32, u32), u32>::default();
            bundle_segs.iter().for_each(|vv| {
                let node = (vv.bundle_id, vv.bundle_dir);
                if (vv.bundle_v_bgn as i64 - vv.bundle_v_end as i64).abs() as f32
                    > (vv.bundle_v_count as f32) * 0.5
                {
                    let e = node_count.entry(node).or_insert(0);
                    *e += 1;
                }
            });
            let mut sort_key = vec![];
            node_length.iter().for_each(|&(_, n)| {
                sort_key.push(*node_count.get(&n).unwrap_or(&0));
            });

            (sort_key, ctg, bundle_segs)
        })
        .collect::<Vec<_>>();

    ctg_data.sort();
    ctg_data.reverse();

    let out_path = Path::new(&args.output_prefix).with_extension("ord");
    let mut out_file = BufWriter::new(File::create(out_path)?);

    ctg_data.into_iter().for_each(|(sort_key, ctg, _)| {
        let sort_key = sort_key
            .into_iter()
            .map(|k| format!("{}", k))
            .collect::<Vec<String>>();
        let sort_key = sort_key.join(",");
        writeln!(
            out_file,
            "{}\t{
        }",
            ctg, sort_key
        )
        .expect("writing error");
    });

    Ok(())
}
