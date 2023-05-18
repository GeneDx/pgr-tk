const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use kodama::{linkage, Method};
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

#[derive(Clone, Copy, Debug)]
enum AlnType {
    Match,
    Deletion,
    Insertion,
}

fn align_bundles(
    q_bundles: &Vec<BundleSegment>,
    t_bundles: &Vec<BundleSegment>,
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
            t_map.insert((q_idx, t_idx), best.0);
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
    while let Some(aln_type) = t_map.get(&(q_idx, t_idx)) {
        // let qq_idx = q_idx;
        // let tt_idx = t_idx;
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
        /*
        println!(
            "{} {} {:?} {:?} {:?} {} {}",
            qq_idx, tt_idx, aln_type, q_bundles[qq_idx].bundle_id, t_bundles[tt_idx].bundle_id, diff_len_delta, max_len_delta
        );
            */
    }
    (diff_len as f32 / max_len as f32, diff_len, max_len)
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
    let mut out_file = BufWriter::new(File::create(out_path).expect("can't create the dist file"));

    let mut dist_map = FxHashMap::<(usize, usize), f32>::default();

    (0..n_ctg)
        .flat_map(|ctg_idx0| (0..n_ctg).map(move |ctg_idx1| (ctg_idx0, ctg_idx1)))
        .for_each(|(ctg_idx0, ctg_idx1)| {
            if ctg_idx0 > ctg_idx1 {
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
                dist_map.insert((ctg_idx0, ctg_idx1), dist);
            }
        });

    let mut dist_mat = vec![];
    (0..n_ctg - 1).for_each(|i| {
        (i + 1..n_ctg).for_each(|j| {
            dist_mat.push(*dist_map.get(&(i, j)).unwrap());
        })
    });
    let dend = linkage(&mut dist_mat, n_ctg, Method::Average);

    let steps = dend.steps().to_vec();
    let mut node_data = FxHashMap::<usize, (String, Vec<usize>, f32)>::default();
    (0..n_ctg).for_each(|ctg_idx| {
        node_data.insert(ctg_idx, (format!("{}", ctg_idx), vec![ctg_idx], 0.0_f32));
    });

    let mut last_node_id = 0_usize;
    steps.iter().enumerate().for_each(|(c, s)| {
        let (node_string1, nodes1, height1) = node_data.remove(&s.cluster1).unwrap();
        let (node_string2, nodes2, height2) = node_data.remove(&s.cluster2).unwrap();
        let new_node_id = c + n_ctg;
        let mut nodes = Vec::<usize>::new();
        let new_node_string = if nodes1.len() > nodes2.len() {
            nodes.extend(nodes1);
            nodes.extend(nodes2);
            format!(
                "({}:{}, {}:{})",
                node_string1,
                s.dissimilarity - height1,
                node_string2,
                s.dissimilarity - height2
            )
        } else {
            nodes.extend(nodes2);
            nodes.extend(nodes1);
            format!(
                "({}:{}, {}:{})",
                node_string2,
                s.dissimilarity - height2,
                node_string1,
                s.dissimilarity - height1
            )
        };
        node_data.insert(new_node_id, (new_node_string, nodes, s.dissimilarity));
        last_node_id = new_node_id;
    });

    let mut tree_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("nwk"))
            .expect("can't create the nwk file"),
    );

    let emptyp_string = ("".to_string(), vec![], 0.0);
    let (tree_string, nodes, _) = node_data.get(&last_node_id).unwrap_or(&emptyp_string);
    writeln!(tree_file, "{};", tree_string).expect("can't write the nwk file");

    let mut dendrogram_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("ddg"))
            .expect("can't create the dendrogram file"),
    );
    let mut node_position_size = FxHashMap::<usize, ((f32, f32), usize)>::default();
    let mut position = 0.0_f32;
    nodes.iter().for_each(|&ctg_idx| {
        node_position_size.insert(ctg_idx, ((position, 0.0), 1));
        writeln!(dendrogram_file, "L\t{}\t{}", ctg_idx, ctg_data[ctg_idx].0)
            .expect("can't write the dendrogram file");
        position += 1.0;
    });
    steps.into_iter().enumerate().for_each(|(c, s)| {
        let ((pos0, _), size0) = *node_position_size.get(&s.cluster1).unwrap();
        let ((pos1, _), size1) = *node_position_size.get(&s.cluster2).unwrap();

        let pos = ((size0 as f32) * pos0 + (size1 as f32) * pos1) / ((size0 + size1) as f32);
        writeln!(
            dendrogram_file,
            "I\t{}\t{}\t{}\t{}\t{}",
            c + n_ctg,
            s.cluster1,
            s.cluster2,
            s.size,
            s.dissimilarity,
        )
        .expect("can't write the dendrogram file");
        node_position_size.insert(c + n_ctg, ((pos, s.dissimilarity), s.size));
    });
    let mut node_positions = node_position_size
        .into_iter()
        .collect::<Vec<(usize, ((f32, f32), usize))>>();
    node_positions.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    node_positions
        .into_iter()
        .for_each(|(vid, ((pos, h), size))| {
            writeln!(dendrogram_file, "P\t{}\t{}\t{}\t{}", vid, pos, h, size)
                .expect("can't write the dendrogram file");
        });
    Ok(())
}
