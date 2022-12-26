const VERSION_STRING: &'static str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::{fs::File, path};

#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-svg")]
#[clap(author, version)]
#[clap(about = "generate alignment scores between contigs using bundle decomposition from a principal bundle bed file", long_about = None)]
struct CmdOptions {
    bed_file_path: String,
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

#[derive(Clone, Copy, Debug)]
enum AlnType {
    Match,
    Deletion,
    Insertion,
}

fn align_bundles(
    q_bundles: &Vec<BundleSegement>,
    t_bundles: &Vec<BundleSegement>,
) -> (f32, usize, usize) {
    let q_count = q_bundles.len();
    let t_count = t_bundles.len();
    let mut s_map = FxHashMap::<(usize, usize), i64>::default();
    let mut t_map = FxHashMap::<(usize, usize), AlnType>::default();

    let mut get_aln_direction_with_best_score =
        |q_idx: usize, t_idx: usize, s_map: &FxHashMap<(usize, usize), i64>| -> (AlnType, i64) {
            let mut best = (AlnType::Match, i64::MIN);
            let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
            let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
            let min_len = if q_len > t_len { t_len } else { q_len };
            let q_b_seg = q_bundles[q_idx]; 
            let t_b_seg = t_bundles[t_idx]; 
            if q_idx == 0 && t_idx == 0 {
                if (q_b_seg.bundle_id == t_b_seg.bundle_id)
                    && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
                {
                    best = (AlnType::Match, 2 * min_len)
                }
            }
            if q_idx > 0 && t_idx > 0 {
                if (q_b_seg.bundle_id == t_b_seg.bundle_id)
                    && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
                {
                    best = (
                        AlnType::Match,
                        2 * min_len + s_map.get(&(q_idx - 1, t_idx - 1)).unwrap(),
                    )
                };
            }
            if t_idx > 0 {
                let insert_score = -2 * q_len + s_map.get(&(q_idx, t_idx - 1)).unwrap();
                if insert_score > best.1 {
                    best = (AlnType::Insertion, insert_score)
                };
            };
            if q_idx > 0 {
                let delete_score = -2 * t_len + s_map.get(&(q_idx - 1, t_idx)).unwrap();
                if delete_score > best.1 {
                    best = (AlnType::Deletion, delete_score)
                }
            }
            t_map.insert((q_idx, t_idx), best.0.clone());
            best
        };

    //let mut best_score = 0;
    //let mut best_q_idx = 0;
    //let mut best_t_idx = 0;

    (0..t_count)
        .flat_map(|t_idx| (0..q_count).map(move |q_idx| (q_idx, t_idx)))
        .for_each(|(q_idx, t_idx)| {
            //println!("{} {}", q_idx, t_idx);
            let (_, score) = get_aln_direction_with_best_score(q_idx, t_idx, &s_map);
            s_map.insert((q_idx, t_idx), score);
            /*
            if score > best_score {
                best_score = score;
                best_q_idx = q_idx;
                best_t_idx = t_idx;
            }
            */
        });
    let mut q_idx = q_count - 1;
    let mut t_idx = t_count - 1;
    let mut diff_len = 0_usize;
    let mut max_len = 1_usize;
    loop {
        if let Some(aln_type) = t_map.get(&(q_idx, t_idx)) {
            // let qq_idx = q_idx;
            // let tt_idx = t_idx;
            let (diff_len_delta, max_len_delta) = match aln_type {
                AlnType::Match => {
                    let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                    let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                    let diff_len_delta = (q_len - t_len).abs() as usize;
                    let max_len_delata = if q_len > t_len {
                        q_len as usize
                    } else {
                        t_len as usize
                    };
                    q_idx -= 1;
                    t_idx -= 1;
                    (diff_len_delta, max_len_delata)
                }
                AlnType::Deletion => {
                    let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                    q_idx -= 1;
                    (q_len as usize, q_len as usize)
                }
                AlnType::Insertion => {
                    let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                    t_idx -= 1;
                    (t_len as usize, t_len as usize)
                }
            };
            diff_len += diff_len_delta;
            max_len += max_len_delta;
            /*
            println!(
                "{} {} {:?} {:?} {:?} {} {}",
                qq_idx, tt_idx, aln_type, q_bundles[qq_idx].bundle_id, t_bundles[tt_idx].bundle_id, diff_len_delta, max_len_delta
            );
            */
        } else {
            break;
        }
    }
    (diff_len as f32 / max_len as f32, diff_len, max_len)
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let bed_file_path = path::Path::new(&args.bed_file_path);
    let bed_file = BufReader::new(File::open(bed_file_path)?);
    let mut ctg_data = FxHashMap::<String, Vec<_>>::default();
    let bed_file_parse_err_msg = "bed file parsing error";
    bed_file.lines().into_iter().for_each(|line| {
        let line = line.unwrap();
        let bed_fields = line.split("\t").collect::<Vec<&str>>();
        let ctg: String = bed_fields[0].to_string();
        let bgn: u32 = bed_fields[1].parse().expect(bed_file_parse_err_msg);
        let end: u32 = bed_fields[2].parse().expect(bed_file_parse_err_msg);
        let pbundle_fields = bed_fields[3].split(":").collect::<Vec<&str>>();
        let bundle_id: u32 = pbundle_fields[0].parse().expect(bed_file_parse_err_msg);
        let bundle_v_count: u32 = pbundle_fields[1].parse().expect(bed_file_parse_err_msg);
        let bundle_dir: u32 = pbundle_fields[2].parse().expect(bed_file_parse_err_msg);
        let bundle_v_bgn: u32 = pbundle_fields[3].parse().expect(bed_file_parse_err_msg);
        let bundle_v_end: u32 = pbundle_fields[4].parse().expect(bed_file_parse_err_msg);

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

    let out_path = Path::new(&args.output_prefix).with_extension("dist");
    let mut out_file = BufWriter::new(File::create(out_path)?);

    (0..n_ctg)
        .flat_map(|ctg_idx0| (0..n_ctg).map(move |ctg_idx1| (ctg_idx0, ctg_idx1)))
        .for_each(|(ctg_idx0, ctg_idx1)| {
            if ctg_idx1 > ctg_idx0 {
                return;
            };
            let (ctg0, bundles0) = &ctg_data[ctg_idx0];
            let (ctg1, bundles1) = &ctg_data[ctg_idx1];
            let (dist0, diff_len0, max_len0) = align_bundles(bundles0, bundles1);
            let (dist1, diff_len1, max_len1) = align_bundles(bundles1, bundles0);
            let (dist, diff_len, max_len) = if dist0 > dist1 {
                (dist0, diff_len0, max_len0)
            } else {
                (dist1, diff_len1, max_len1)
            };
            writeln!(
                out_file,
                "{} {} {} {} {}",
                ctg0, ctg1, dist, diff_len, max_len
            )
            .expect("writing error");
            if ctg_idx1 != ctg_idx0 {
                writeln!(
                    out_file,
                    "{} {} {} {} {}",
                    ctg1, ctg0, dist, diff_len, max_len
                )
                .expect("writing error");
            }
        });
    Ok(())
}
