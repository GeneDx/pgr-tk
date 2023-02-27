const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::{fs::File, path};

/// Generate alignment scores between sequences using bundle decomposition from a principal bundle bed file
#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-bed2dist")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the path to the pricipal bundle bed file
    bed_file_path: String,
    /// the prefix of the output file
    output_prefix: String,
    /// the path the annotation files
    #[clap(long)]
    ctgs_of_interest: Option<String>,
    /// specify the bundle to anchor on
    #[clap(long, default_value_t = false)]
    alt_anchoring_mode: bool,
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
) -> (
    f32,
    usize,
    usize,
    Vec<(usize, usize, AlnType, u32, u32, i64)>,
) {
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
                } else {
                    best = (AlnType::Match, 0)
                }
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
                let insert_score = -q_len + s_map.get(&(q_idx, t_idx - 1)).unwrap();
                if insert_score > best.1 {
                    best = (AlnType::Insertion, insert_score)
                };
            };
            if q_idx > 0 {
                let delete_score = -t_len + s_map.get(&(q_idx - 1, t_idx)).unwrap();
                if delete_score > best.1 {
                    best = (AlnType::Deletion, delete_score)
                }
            }
            t_map.insert((q_idx, t_idx), best.0);
            best
        };

    //let mut best_score = 0;
    //let mut best_q_idx = 0;
    //let mut best_t_idx = 0;
    let mut aln_path = Vec::<_>::new();
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
        let score = s_map.get(&(qq_idx, tt_idx)).unwrap_or(&0);
        aln_path.push((
            qq_idx,
            tt_idx,
            aln_type.clone(),
            q_bundles[qq_idx].bundle_id,
            t_bundles[tt_idx].bundle_id,
            *score,
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

    let mut ctg_to_annotation = FxHashMap::<String, String>::default();
    let ctg_data_vec = if args.ctgs_of_interest.is_some() {
        let filename = args.ctgs_of_interest.unwrap();
        let path = path::Path::new(&filename);
        let ctg_of_interest_file = BufReader::new(File::open(path)?);
        let ctg_data_vec: Vec<_> = ctg_of_interest_file
            .lines()
            .map(|line| {
                let line = line.unwrap().trim().to_string();
                if line.is_empty() {
                    return None;
                }
                if &line[0..1] == "#" {
                    return None;
                }
                let mut ctg_annotation = line.split('\t');
                let ctg = ctg_annotation
                    .next()
                    .expect("error parsthe ctgs_of_interest file")
                    .to_string();

                // println!("{:?} {:?}", line, ctg);
                let data = ctg_data.get(&ctg).unwrap().to_owned();
                if let Some(annotation) = ctg_annotation.next() {
                    ctg_to_annotation.insert(ctg.clone(), annotation.to_string());
                    //annotations.push( (ctg.clone(), annotation) );
                    Some((ctg, annotation.to_string(), data))
                } else {
                    ctg_to_annotation.insert(ctg.clone(), "".to_string());

                    Some((ctg, "".to_string(), data))
                }
            })
            .flatten()
            .collect();
        ctg_data_vec
    } else {
        let mut ctg_data_vec = ctg_data.iter().map(|(k, v)| (k, v)).collect::<Vec<_>>();
        ctg_data.keys().into_iter().for_each(|ctg| {
            ctg_to_annotation.insert(ctg.clone(), ctg.clone());
        });
        ctg_data_vec.sort();
        ctg_data_vec
            .into_iter()
            .map(|(ctg, data)| (ctg.clone(), ctg.clone(), data.clone()))
            .collect()
    };

    let n_ctg = ctg_data_vec.len();

    let out_path = Path::new(&args.output_prefix).with_extension("offset");
    let mut out_file =
        BufWriter::new(File::create(out_path).expect("can't create the bundle offset file"));

    let (ctg1, _anootation, bundles1) = &ctg_data_vec[0];
    let _ = writeln!(out_file, "{}\t{}", ctg1, 0);
    (1..n_ctg).for_each(|ctg_idx| {
        let (ctg0, _annotation, bundles0) = &ctg_data_vec[ctg_idx];
        let (_dist0, _diff_len0, _max_len0, alns) = align_bundles(bundles0, bundles1);
        let mut best_anchor_point = None;
        let mut best_single_match_anchor_point = None;

        let mut score = 0_i64;
        let mut last_global_score = 0_i64;
        let mut current_score = 0_i64;
        let mut best_score = 0_i64;
        let mut best_single_bundle_score = 0_i64;
        alns.into_iter().for_each(
            |(qq_idx, tt_idx, _aln_type, _q_bid, _t_bid, global_score)| {
                score = global_score - last_global_score;
                if score > best_single_bundle_score {
                    best_single_bundle_score = score;
                    best_single_match_anchor_point = Some((qq_idx, tt_idx));
                };
                current_score += score;
                if current_score < 0 {
                    current_score = 0;
                }
                if current_score > best_score {
                    best_score = current_score;
                    best_anchor_point = Some((qq_idx, tt_idx));
                }
                /*
                println!(
                    "{} {} {} {} {:?} {:?} {:?} {} {}",
                    ctg0, ctg1, qq_idx, tt_idx, aln_type, q_bid, t_bid, global_score, current_score
                );
                */
                last_global_score = global_score;
            },
        );

        let b0 = if args.alt_anchoring_mode {
            if let Some(anchor_point) = best_single_match_anchor_point {
                bundles0[anchor_point.0].bgn
            } else {
                0
            }
        } else {
            if let Some(anchor_point) = best_anchor_point {
                bundles0[anchor_point.0].bgn
            } else {
                0
            }
        };
        let b1 = if args.alt_anchoring_mode {
            if let Some(anchor_point) = best_single_match_anchor_point {
                bundles1[anchor_point.1].bgn
            } else {
                0
            }
        } else {
            if let Some(anchor_point) = best_anchor_point {
                bundles1[anchor_point.1].bgn
            } else {
                0
            }
        };
        let offset = b1 as i64 - b0 as i64;
        //println!("XXX {} {} {} {:?} {:?}", best_score, best_anchor_point.0, best_anchor_point.1, bundles0[best_anchor_point.0], bundles1[best_anchor_point.1]);
        let _ = writeln!(out_file, "{}\t{}", ctg0, offset);
    });

    Ok(())
}
