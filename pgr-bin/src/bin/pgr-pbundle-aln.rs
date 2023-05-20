const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::{fs::File, path};

/// Generate alignment between sequences using bundle decomposition from a principal bundle bed file
#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-bed2dist")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the path to the principal bundle bed file
    bed_file_path: String,
    /// a file contain two lines of the contig ids that should be aligned to each other
    aln_spec: String,
    /// the prefix of the output file
    output_prefix: String,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
struct BundleSegment {
    bgn: u32,
    end: u32,
    bundle_id: u32,
    bundle_v_count: u32,
    bundle_dir: u32,
    bundle_v_bgn: u32,
    bundle_v_end: u32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum AlnType {
    Match,
    Insertion,
    Deletion,
}

type AlnPathElement = (usize, usize, AlnType, u32, u32, usize, usize);
type AlnPath = Vec<AlnPathElement>;

fn align_bundles(
    q_bundles: &Vec<BundleSegment>,
    t_bundles: &Vec<BundleSegment>,
) -> (f32, usize, usize, AlnPath) {
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
            if q_idx == 0
                && t_idx == 0
                && (q_b_seg.bundle_id == t_b_seg.bundle_id)
                && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
            {
                best = (AlnType::Match, 2 * min_len)
            };
            if q_idx > 0
                && t_idx > 0
                && q_b_seg.bundle_id == t_b_seg.bundle_id
                && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
            {
                best = (
                    AlnType::Match,
                    2 * min_len + s_map.get(&(q_idx - 1, t_idx - 1)).unwrap(),
                )
            };
            if t_idx > 0 {
                let score = -2 * q_len + s_map.get(&(q_idx, t_idx - 1)).unwrap();
                if score > best.1 {
                    best = (AlnType::Deletion, score)
                };
            };
            if q_idx > 0 {
                let score = -2 * t_len + s_map.get(&(q_idx - 1, t_idx)).unwrap();
                if score > best.1 {
                    best = (AlnType::Insertion, score)
                }
            }
            t_map.insert((q_idx, t_idx), best.0);
            best
        };

    //let mut best_score = 0;
    //let mut best_q_idx = 0;
    //let mut best_t_idx = 0;
    let mut aln_path = AlnPath::new();

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
    while let Some(aln_type) = t_map.get(&(q_idx, t_idx)) {
        let qq_idx = q_idx;
        let tt_idx = t_idx;
        let (diff_len_delta, max_len_delta) = match aln_type {
            AlnType::Match => {
                let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                let diff_len_delta = (q_len - t_len).unsigned_abs() as usize;
                let max_len_delata = if q_len > t_len {
                    q_len as usize
                } else {
                    t_len as usize
                };
                q_idx -= 1;
                t_idx -= 1;
                (diff_len_delta, max_len_delata)
            }
            AlnType::Insertion => {
                let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                q_idx -= 1;
                (q_len as usize, q_len as usize)
            }
            AlnType::Deletion => {
                let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                t_idx -= 1;
                (t_len as usize, t_len as usize)
            }
        };
        diff_len += diff_len_delta;
        max_len += max_len_delta;
        aln_path.push((
            qq_idx,
            tt_idx,
            aln_type.clone(),
            q_bundles[qq_idx].bundle_id,
            t_bundles[tt_idx].bundle_id,
            diff_len_delta,
            max_len_delta,
        ));
    }
    aln_path.reverse();
    (
        diff_len as f32 / max_len as f32,
        diff_len,
        max_len,
        aln_path,
    )
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let bed_file_path = path::Path::new(&args.bed_file_path);
    let bed_file = BufReader::new(File::open(bed_file_path).expect("can't open the bed file"));
    let mut ctg_data = FxHashMap::<String, Vec<_>>::default();
    let bed_file_parse_err_msg = "bed file parsing error";
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
        let b_seg = BundleSegment {
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

    let aln_spec = path::Path::new(&args.aln_spec);
    let spec_file = BufReader::new(File::open(aln_spec).expect("can't open the aln_spec file"));
    let mut ctg_of_interests = Vec::<String>::new();
    spec_file.lines().into_iter().for_each(|line| {
        let line = line.unwrap().trim().to_string();
        ctg_of_interests.push(line);
    });

    let mut ctg_data = ctg_of_interests
        .into_iter()
        .map(|k| {
            let v = ctg_data
                .get(&k)
                .expect(format!("ctg name nof found: {}", k).as_str());
            (k, v)
        })
        .collect::<Vec<_>>();

    let n_ctg = ctg_data.len();

    let out_path = Path::new(&args.output_prefix).with_extension("bdl.aln");
    let mut out_file =
        BufWriter::new(File::create(out_path).expect("can't create the bundle alignment file"));

    (1..n_ctg).into_iter().for_each(|ctg_idx1| {
        let ctg_idx0 = 0;

        // the first sequence is the "target"
        let (target_ctg, target_bundles) = &ctg_data[ctg_idx0];
        let (query_ctg, query_bundles) = &ctg_data[ctg_idx1];
        let (dist0, diff_len0, max_len0, aln_path) = align_bundles(query_bundles, target_bundles);

        aln_path.into_iter().for_each(
            |(
                qq_idx,
                tt_idx,
                aln_type,
                q_bundle_id,
                t_bundle_id,
                diff_len_delta,
                max_len_delta,
            )| {
                let target_data = target_bundles.get(tt_idx).unwrap();
                let query_data = query_bundles.get(qq_idx).unwrap();
                if aln_type == AlnType::Match
                    && t_bundle_id == q_bundle_id
                    && target_data.bundle_v_end == query_data.bundle_v_end
                    && target_data.bundle_v_bgn == query_data.bundle_v_bgn
                {
                    println!(
                        "* {} {} {:?} {} {} {:?} {:?} {} {}",
                        tt_idx,
                        qq_idx,
                        aln_type,
                        t_bundle_id,
                        q_bundle_id,
                        target_bundles.get(tt_idx).unwrap(),
                        query_bundles.get(qq_idx).unwrap(),
                        diff_len_delta,
                        max_len_delta
                    )
                } else {
                    println!(
                        "- {} {} {:?} {} {} {:?} {:?} {} {}",
                        tt_idx,
                        qq_idx,
                        aln_type,
                        t_bundle_id,
                        q_bundle_id,
                        target_bundles.get(tt_idx).unwrap(),
                        query_bundles.get(qq_idx).unwrap(),
                        diff_len_delta,
                        max_len_delta
                    )
                };
            },
        );
    });
    Ok(())
}
