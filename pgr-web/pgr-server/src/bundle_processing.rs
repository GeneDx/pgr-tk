use std::io::{BufWriter, Write};
use std::sync::Arc;

use pgr_db::ext::{get_principal_bundle_decomposition, SeqIndexDB};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use svg::node::{self, element, Node};
use svg::Document;

static CMAP: [&str; 97] = [
    "#870098", "#00aaa5", "#3bff00", "#ec0000", "#00a2c3", "#00f400", "#ff1500", "#0092dd",
    "#00dc00", "#ff8100", "#007ddd", "#00c700", "#ffb100", "#0038dd", "#00af00", "#fcd200",
    "#0000d5", "#009a00", "#f1e700", "#0000b1", "#00a55d", "#d4f700", "#4300a2", "#00aa93",
    "#a1ff00", "#dc0000", "#00aaab", "#1dff00", "#f40000", "#009fcb", "#00ef00", "#ff2d00",
    "#008ddd", "#00d700", "#ff9900", "#0078dd", "#00c200", "#ffb900", "#0025dd", "#00aa00",
    "#f9d700", "#0000c9", "#009b13", "#efed00", "#0300aa", "#00a773", "#ccf900", "#63009e",
    "#00aa98", "#84ff00", "#e10000", "#00a7b3", "#00ff00", "#f90000", "#009bd7", "#00ea00",
    "#ff4500", "#0088dd", "#00d200", "#ffa100", "#005ddd", "#00bc00", "#ffc100", "#0013dd",
    "#00a400", "#f7dd00", "#0000c1", "#009f33", "#e8f000", "#1800a7", "#00aa88", "#c4fc00",
    "#78009b", "#00aaa0", "#67ff00", "#e60000", "#00a4bb", "#00fa00", "#fe0000", "#0098dd",
    "#00e200", "#ff5d00", "#0082dd", "#00cc00", "#ffa900", "#004bdd", "#00b400", "#ffc900",
    "#0000dd", "#009f00", "#f4e200", "#0000b9", "#00a248", "#dcf400", "#2d00a4", "#00aa8d",
    "#bcff00",
];
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MatchSummary {
    pub q_bgn: u32,
    pub q_end: u32,
    pub t_bgn: u32,
    pub t_end: u32,
    pub num_hits: usize,
    pub reversed: bool,
}

#[derive(Serialize, Deserialize)]
pub struct TargetMatchPrincipalBundles {
    pub query: SequenceQuerySpec,
    pub match_summary: Vec<(u32, Vec<MatchSummary>)>, // (q_id, vec[(q_bgn, q_end, t_bgn, t_end, num_hits, reversed)])
    pub sid_ctg_src: Vec<(u32, String, String)>,
    pub bundle_bed_records: Vec<Vec<PrincipalBundleBedRecord>>,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct PrincipalBundleBedRecord {
    pub ctg: String,
    pub bgn: u32,
    pub end: u32,
    pub b_id: u32,
    pub b_size: usize,
    pub b_direction: u32,
    pub b_bgn: usize,
    pub b_end: usize,
    pub r_type: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ShmmrSpec {
    pub w: u32,
    pub k: u32,
    pub r: u32,
    pub min_span: u32,
    pub sketch: bool,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct SequenceQuerySpec {
    pub source: String,
    pub ctg: String,
    pub bgn: usize,
    pub end: usize,
    pub padding: usize,
    pub merge_range_tol: usize,
    //pub pb_shmmr_spec: ShmmrSpec,
    pub w: u32,
    pub k: u32,
    pub r: u32,
    pub min_span: u32,
    pub sketch: bool,
    pub min_cov: usize,
    pub min_branch_size: usize,
    pub bundle_length_cutoff: usize,
    pub bundle_merge_distance: usize,
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

pub fn get_target_and_principal_bundle_decomposition(
    seq_query_spec: &SequenceQuerySpec,
    seq_db: Arc<SeqIndexDB>,
) -> Option<TargetMatchPrincipalBundles> {
    let sample_name = seq_query_spec.source.clone();
    let ctg_name = seq_query_spec.ctg.clone();
    let padding = seq_query_spec.padding;
    let merge_range_tol = seq_query_spec.merge_range_tol;
    let seq_len = match seq_db
        .seq_index
        .as_ref()
        .unwrap()
        .get(&(ctg_name.clone(), Some(sample_name.clone())))
    {
        None => 0,
        Some(value) => value.1,
    };

    let _q_seq_len = seq_query_spec.end - seq_query_spec.bgn;
    let q_seq_bgn = if padding > seq_query_spec.bgn {
        0
    } else {
        seq_query_spec.bgn - padding
    };
    let q_seq_end = if seq_query_spec.end + padding > seq_len as usize {
        seq_len as usize
    } else {
        seq_query_spec.end + padding
    };

    let sub_seq = seq_db
        .get_sub_seq(sample_name.clone(), ctg_name.clone(), q_seq_bgn, q_seq_end)
        .unwrap();

    /*
    println!(
        "DBG: sub_seq_len {:?} {} {}",
        sub_seq.len(),
        q_seq_bgn,
        q_seq_end
    );
     */

    let query_results = seq_db.query_fragment_to_hps_from_mmap_file(
        sub_seq.clone(),
        0.25,
        Some(128),
        Some(128),
        Some(128),
        Some(0),
    );

    let aln_range = if let Some(qr) = query_results {
        let mut sid_to_alns = FxHashMap::default();
        qr.into_iter().for_each(|(sid, alns)| {
            let mut aln_lens = vec![];
            let mut f_count = 0_usize;
            let mut r_count = 0_usize;
            alns.into_iter().for_each(|(_score, aln)| {
                if aln.len() > 2 {
                    aln_lens.push(aln.len());
                    for hp in &aln {
                        if hp.0 .2 == hp.1 .2 {
                            f_count += 1;
                        } else {
                            r_count += 1;
                        }
                    }
                    let orientation = if f_count > r_count { 0_u32 } else { 1_u32 };
                    let e = sid_to_alns.entry(sid).or_insert_with(Vec::new);
                    e.push((aln, orientation))
                }
            })
        });

        let mut aln_range = FxHashMap::default();
        sid_to_alns.into_iter().for_each(|(sid, alns)| {
            alns.into_iter().for_each(|(aln, orientation)| {
                let mut target_coordinates = aln
                    .iter()
                    .map(|v| (v.1 .0, v.1 .1))
                    .collect::<Vec<(u32, u32)>>();
                target_coordinates.sort();
                let bgn = target_coordinates[0].0;
                let end = target_coordinates[target_coordinates.len() - 1].1;
                let e = aln_range.entry(sid).or_insert_with(Vec::new);
                e.push((bgn, end, end - bgn, orientation, aln));
            })
        });

        // merge aln_range
        let aln_range = aln_range
            .into_iter()
            .map(|(sid, rgns)| {
                let mut f_rgns = rgns
                    .iter()
                    .filter(|&v| v.3 == 0)
                    .cloned()
                    .collect::<Vec<_>>();

                let mut r_rgns = rgns
                    .iter()
                    .filter(|&v| v.3 == 1)
                    .cloned()
                    .collect::<Vec<_>>();

                f_rgns.sort();
                r_rgns.sort();

                let mut out_rgns = vec![];
                let mut last_rgn: (u32, u32, u32, u32, Vec<_>) = (0, 0, 0, 0, vec![]);
                f_rgns.into_iter().for_each(|r| {
                    if last_rgn.4.is_empty() {
                        last_rgn = r;
                    } else {
                        let l_bgn = last_rgn.0;
                        let l_end = last_rgn.1;
                        assert!(l_end > l_bgn);
                        let r_bgn = r.0;
                        let r_end = r.1;
                        if (r_bgn as i64) - (l_end as i64) < merge_range_tol as i64 {
                            let bgn = l_bgn;
                            let end = if r_end > l_end { r_end } else { l_end };
                            let len = end - bgn;
                            let orientation = last_rgn.3;
                            let mut aln = last_rgn.4.clone();
                            aln.extend(r.4);
                            last_rgn = (bgn, end, len, orientation, aln);
                        } else {
                            out_rgns.push(last_rgn.clone());
                            last_rgn = r;
                        }
                    }
                });
                if last_rgn.2 > 0 {
                    //not empty
                    out_rgns.push(last_rgn);
                };

                let mut last_rgn: (u32, u32, u32, u32, Vec<_>) = (0, 0, 0, 0, vec![]);
                r_rgns.into_iter().for_each(|r| {
                    if last_rgn.4.is_empty() {
                        last_rgn = r;
                    } else {
                        let l_bgn = last_rgn.0;
                        let l_end = last_rgn.1;
                        assert!(l_end > l_bgn);
                        let r_bgn = r.0;
                        let r_end = r.1;
                        if (r_bgn as i64) - (l_end as i64) < merge_range_tol as i64 {
                            let bgn = l_bgn;
                            let end = if r_end > l_end { r_end } else { l_end };
                            let len = end - bgn;
                            let orientation = last_rgn.3;
                            let mut aln = last_rgn.4.clone();
                            aln.extend(r.4);
                            last_rgn = (bgn, end, len, orientation, aln);
                        } else {
                            out_rgns.push(last_rgn.clone());
                            last_rgn = r;
                        }
                    }
                });
                if last_rgn.2 > 0 {
                    //not empty
                    out_rgns.push(last_rgn);
                };

                (sid, out_rgns)
            })
            .collect::<FxHashMap<_, _>>();
        Some(aln_range)
    } else {
        None
    };

    let sid_ctg_src = if let Some(&ref aln_range) = aln_range.as_ref() {
        aln_range
            .keys()
            .into_iter()
            .map(|sid| {
                let (ctg, src, _ctg_len) = seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
                let src = (*src).as_ref().unwrap_or(&"N/A".to_string()).clone();
                (*sid, ctg.clone(), src)
            })
            .collect::<Vec<(u32, String, String)>>()
    } else {
        vec![]
    };

    let (sub_seq_range_for_fasta, match_summary) = if let Some(aln_range) = aln_range {
        let mut sub_seq_range_for_fasta = Vec::<(u32, u32, u32, u32, String)>::new();
        let match_summary = aln_range
            .into_iter()
            .map(|(sid, rgns)| {
                let (ctg, _src, _ctg_len) = seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
                let hits = rgns
                    .into_iter()
                    .map(|(b, e, _, orientation, mut aln)| {
                        aln.sort();
                        let q_bgn = aln[0].0 .0;
                        let q_end = aln[aln.len() - 1].0 .1;

                        MatchSummary {
                            q_bgn,
                            q_end,
                            t_bgn: b,
                            t_end: e,
                            num_hits: aln.len(),
                            reversed: orientation == 1,
                        }
                    })
                    .filter(|v| {
                        (v.num_hits > 100)
                            & ((v.t_end - v.t_bgn) as f32 / (v.q_end - v.q_bgn) as f32 > 0.6)
                    })
                    .collect::<Vec<MatchSummary>>();
                hits.iter().for_each(|v| {
                    sub_seq_range_for_fasta.push((
                        sid,
                        v.t_bgn,
                        v.t_end,
                        if v.reversed { 1 } else { 0 },
                        ctg.clone(),
                    ))
                });
                (sid, hits)
            })
            .collect::<Vec<(u32, Vec<MatchSummary>)>>();
        (sub_seq_range_for_fasta, match_summary)
    } else {
        (vec![], vec![])
    };

    let seq_list = sub_seq_range_for_fasta
        .par_iter()
        .map(|(sid, b, e, orientation, target_seq_name)| {
            let target_seq = seq_db
                .get_sub_seq_by_id(*sid, *b as usize, *e as usize)
                .unwrap();
            let target_seq = if *orientation == 1 {
                pgr_db::fasta_io::reverse_complement(&target_seq)
            } else {
                target_seq
            };
            (target_seq_name.into(), target_seq)
        })
        .collect::<Vec<(String, Vec<u8>)>>();

    let mut new_seq_db = SeqIndexDB::new();
    let shmmr_spec = ShmmrSpec {
        w: seq_query_spec.w,
        k: seq_query_spec.k,
        r: seq_query_spec.r,
        min_span: seq_query_spec.r,
        sketch: seq_query_spec.sketch,
    };

    new_seq_db
        .load_from_seq_list(
            seq_list,
            "Memory".into(),
            shmmr_spec.w,
            shmmr_spec.k,
            shmmr_spec.r,
            shmmr_spec.min_span,
        )
        .expect("can't load seq_db");

    let (principal_bundles_with_id, vertex_to_bundle_id_direction_pos) = new_seq_db
        .get_principal_bundles_with_id(
            seq_query_spec.min_cov,
            seq_query_spec.min_branch_size,
            None,
        );

    let bid_to_size = principal_bundles_with_id
        .iter()
        .map(|v| (v.0, v.2.len()))
        .collect::<FxHashMap<usize, usize>>();

    let sid_smps =
        get_principal_bundle_decomposition(&vertex_to_bundle_id_direction_pos, &new_seq_db);
    let sid_smps: FxHashMap<u32, Vec<_>> = sid_smps.into_iter().collect();

    let mut seq_info = new_seq_db
        .seq_info
        .unwrap()
        .into_iter()
        .map(|(k, v)| (k, v))
        .collect::<Vec<_>>();

    let mut repeat_count = FxHashMap::<u32, Vec<u32>>::default();
    let mut non_repeat_count = FxHashMap::<u32, Vec<u32>>::default();
    seq_info.sort_by_key(|k| k.1 .0.clone());

    let bundle_bed_records = seq_info
        .iter()
        .map(|(sid, sdata)| {
            let (ctg, _src, _len) = sdata;
            let smps = sid_smps.get(&sid).unwrap();
            let smp_partitions = group_smps_by_principle_bundle_id(
                smps,
                seq_query_spec.bundle_length_cutoff,
                seq_query_spec.bundle_merge_distance,
            );
            let mut ctg_bundle_count = FxHashMap::<usize, usize>::default();
            smp_partitions.iter().for_each(|p| {
                let bid = p[0].1;
                *ctg_bundle_count.entry(bid).or_insert_with(|| 0) += 1;
            });
            let pbundle_bed_records = smp_partitions
                .into_iter()
                .map(|p| {
                    let b = p[0].0 .2;
                    let e = p[p.len() - 1].0 .3 + shmmr_spec.k;
                    let bid = p[0].1;
                    let direction = p[0].2;
                    let is_repeat;
                    if *ctg_bundle_count.get(&bid).unwrap_or(&0) > 1 {
                        repeat_count
                            .entry(*sid)
                            .or_insert_with(|| vec![])
                            .push(e - b - shmmr_spec.k);
                        is_repeat = "R";
                    } else {
                        non_repeat_count
                            .entry(*sid)
                            .or_insert_with(|| vec![])
                            .push(e - b - shmmr_spec.k);
                        is_repeat = "U";
                    };
                    PrincipalBundleBedRecord {
                        ctg: ctg.clone(),
                        bgn: b,
                        end: e,
                        b_id: bid as u32,
                        b_size: bid_to_size[&bid],
                        b_direction: direction,
                        b_bgn: p[0].3,
                        b_end: p[p.len() - 1].3,
                        r_type: is_repeat.to_string(),
                    }
                })
                .collect::<Vec<PrincipalBundleBedRecord>>();
            pbundle_bed_records
        })
        .collect::<Vec<Vec<PrincipalBundleBedRecord>>>();

    Some(TargetMatchPrincipalBundles {
        query: (*seq_query_spec).clone(),
        match_summary,
        sid_ctg_src,
        bundle_bed_records,
    })
}

pub fn pb_data_to_html_string(targets: &TargetMatchPrincipalBundles) -> String {
    let mut target_lenths = targets
        .match_summary
        .iter()
        .flat_map(|v| v.1.iter().map(|v| v.t_end - v.t_bgn).collect::<Vec<u32>>())
        .collect::<Vec<u32>>();

    target_lenths.sort();
    let max_length = target_lenths
        .get(target_lenths.len() - 1 )
        .unwrap_or(&200000);


    let ctg_data_vec = targets.bundle_bed_records.iter().map(|v| {
        let b_segements = v
            .iter()
            .map(|r| (r.bgn, r.end, r.b_id, r.b_direction))
            .collect::<Vec<_>>();
        let ctg = if let Some(r) = v.get(0) {
            r.ctg.clone()
        } else {
            "NA".to_string()
        };
        (ctg.clone(), ctg.clone(), b_segements)
    });

    let track_scaling = 1.0;
    let stroke_width = 1.0;
    let left_padding = 50.0;
    let no_tooltips = false;
    let highlight_repeats = 1.2;
    let mut y_offset = 0.0_f32;
    let track_range = *max_length as f32 * 1.05;
    //let track_range = 1200000.0;
    let track_panel_width = 1200.0;
    let annotation_panel_width = 800.0;
    let tree_width = 0.0;
    let h_factor = 1.5;
    let scaling_factor = track_panel_width as f32 / (track_range + 2.0 * left_padding) as f32;
    let delta_y = 16.0_f32 * track_scaling;

    let mut bundle_class_styles = FxHashMap::<String, String>::default();

    let ctg_with_svg_paths: Vec<(String, (Vec<element::Group>, element::Text))>  = ctg_data_vec
        .map(|(ctg, annotation, bundle_segment)| {
            let mut bundle_segment_count = FxHashMap::<u32, usize>::default();
            bundle_segment.iter().for_each(|&(_bgn, _end, bundle_id, _direction)| {
                let e = bundle_segment_count.entry(bundle_id).or_insert_with(|| 0);
                *e += 1;
            });

            // let offset = *ctg_to_offset.get(&ctg).unwrap_or(&0);
            let offset = 0;
            let paths: Vec<element::Group> = bundle_segment
                .into_iter()
                .map(|(bgn0, end0, bundle_id, direction)| {
                    let mut bgn = (bgn0 as i64 + offset) as f32 * scaling_factor;
                    let mut end = (end0 as i64 + offset) as f32 * scaling_factor;
                    if direction == 1 {
                        (bgn, end) = (end, bgn);
                    }

                    let arror_end = end as f32;
                    let halfwidth = 5.0 * track_scaling;
                    let end =
                        if direction == 0 {
                            if end as f32 - halfwidth < bgn {
                                bgn
                            } else {
                                end as f32 - halfwidth
                            }
                        } else if end as f32 + halfwidth > bgn {
                            bgn
                        } else {
                            end as f32 + halfwidth
                        };

                    let bottom0 = -halfwidth * 0.6;
                    let top0 = halfwidth * 0.6;
                    let bottom1 = -halfwidth * 0.8;
                    let top1 = halfwidth * 0.8;
                    let center = 0 as f32;

                    let bundle_class = format!("bundle_{bundle_id:05}");
                    let bundle_rep_class = format!("bundle_{bundle_id:05} repeat");

                    let bundle_color = CMAP[((bundle_id * 57) % 59) as usize];
                    let stroke_color = CMAP[93 - ((bundle_id * 31) % 47) as usize];
                    let css_string = format!(
r#".{bundle_class} {{fill:{bundle_color}; stroke:{stroke_color}; stroke-width:{stroke_width}; fill-opacity:0.5}}"#);
                    bundle_class_styles.entry(bundle_class.clone()).or_insert(css_string);

                    let bundle_class = if *bundle_segment_count.get(&bundle_id).unwrap_or(&0) > 1 && highlight_repeats > 1.0001 {
                        bundle_rep_class
                    } else {
                        bundle_class
                    };

                    let path_str = format!(
					"M {bgn} {bottom0} L {bgn} {top0} L {end} {top0} L {end} {top1} L {arror_end} {center} L {end} {bottom1} L {end} {bottom0} Z");
                    let mut p = element::Path::new()
                        .set("d", path_str)
                        .set("class", "bundle ".to_string() + bundle_class.as_str());
                    let mut g = element::Group::new().set("transform", format!("translate({left_padding} {y_offset})"));
                    if !no_tooltips { // it may be good idea to disable it for every large region visualization
                        p.append(element::Title::new().add(node::Text::new(format!("{}:{}-{}:{}", ctg, bgn0, end0, bundle_id ))));
                    };
                    g.append(p);
                    g
                })
                .collect();

            let text = element::Text::new()
                .set("x", 20.0 + left_padding + track_range as f32 * scaling_factor)
                .set("y", y_offset + 2.0)
                .set("font-size", "10px")
                .set("font-family", "monospace")
                .add(node::Text::new(annotation));
            y_offset += delta_y;
            (ctg, (paths, text))
        })
        .collect();

    // start to construct the SVG element

    let mut document = Document::new()
        .set(
            "viewBox",
            (
                -tree_width,
                -32,
                tree_width + track_panel_width as f32 + annotation_panel_width as f32,
                24.0 + y_offset,
            ),
        )
        .set(
            "width",
            tree_width + track_panel_width as f32 + annotation_panel_width as f32,
        )
        .set("height", 56.0 + y_offset)
        .set("preserveAspectRatio", "none")
        .set("id", "bundleViwer");

    // insert CSS
    let stroke_width_rep = stroke_width * highlight_repeats;
    let stroke_width_hover = stroke_width * 2.0;
    let stroke_width_hover_rep = stroke_width_rep * 2.0;
    let mut css_strings = vec![
        format!(".repeat {{stroke-width:{stroke_width_rep};}}"),
        format!(".bundle:hover {{ stroke-width:{stroke_width_hover};}}"),
        format!(".repeat:hover {{ stroke-width:{stroke_width_hover_rep};}}"),
        format!(".region {{ stroke-opacity: 0.5 }};"),
    ];
    css_strings.extend(bundle_class_styles.values().cloned());
    let h_factor = h_factor;
    css_strings.push(format!(
        r#"path.highlighted {{transform: scaleY({h_factor}); fill-opacity:1}}"#
    ));
    let style = element::Style::new(css_strings.join("\n")).set("type", "text/css");
    document.append(style);

    let track_tick_interval: Option<usize> = None;
    let track_tick_interval = if let Some(tick_interval) = track_tick_interval {
        tick_interval
    } else {
        let mut tick_interval = 1_usize;
        let mut tmp = track_range as f32;
        tmp = tmp * 0.1;
        while tmp > 1.01 {
            tick_interval *= 10;
            tmp = tmp * 0.1;
        }
        tick_interval
    };

    let mut tickx = track_tick_interval;
    loop {
        if tickx > track_range as usize {
            break;
        }
        let x = tickx as f32 * scaling_factor + left_padding;
        let tick_path_str = format!("M {x} -16 L {x} -20");
        let tick_path = element::Path::new()
            .set("stroke", "#000")
            .set("fill", "none")
            .set("stroke-width", 1)
            .set("d", tick_path_str);
        document.append(tick_path);
        tickx += track_tick_interval;
    }

    let text = element::Text::new()
        .set("x", 20.0 + left_padding + scaling_factor)
        .set("y", -14)
        .set("font-size", "10px")
        .set("font-family", "sans-serif")
        .add(node::Text::new(format!("{} bps", track_range)));
    document.append(text);

    // insert the bundle paths
    ctg_with_svg_paths
        .into_iter()
        .for_each(|(_ctg, (paths, text))| {
            // println!("{}", ctg);
            document.append(text);
            paths.into_iter().for_each(|path| document.append(path));
        });

    let out_str = Vec::new();
    let mut out_file = BufWriter::new(out_str);
    let msg = "can't write the HTML doc";
    writeln!(out_file, "<html><body>").expect(msg);
    let jscript = r#"
<script>
document.addEventListener('readystatechange', event => {
    if (event.target.readyState === "complete") {
        var bundles = document.getElementsByClassName("bundle");
        for (let i = 0; i < bundles.length; i++) {
            bundles[i].onclick = function (e) {
                // alert(e.target.classList);
                let is_highlighted = false;
                let bundle_id = "";
                for (let cidx = 0; cidx < e.target.classList.length; cidx++) {
                    if (e.target.classList[cidx] == "highlighted") {
                        is_highlighted = true;
                    }
                    if (e.target.classList[cidx].match("bundle_")) {
                        bundle_id = e.target.classList[cidx]
                    }
                };
                var bundles2 = document.getElementsByClassName(bundle_id);
                for (let j = 0; j < bundles2.length; j++) {
                    if (is_highlighted) {
                        bundles2[j].classList.remove("highlighted");
                    } else {
                        bundles2[j].classList.add("highlighted");
                    }
                }
            };
        }
    }
});
</script>
"#;
    writeln!(out_file, "{}", jscript).expect(msg);
    let mut svg_elment = BufWriter::new(Vec::new());
    svg::write(&mut svg_elment, &document).unwrap();
    writeln!(
        out_file,
        "{}",
        String::from_utf8_lossy(&svg_elment.into_inner().unwrap())
    )
    .expect(msg);
    writeln!(out_file, "</body></html>").expect(msg);

    let _ = out_file.flush();
    let out_str = out_file.into_inner().unwrap();
    String::from_utf8_lossy(&out_str[..]).to_string()
}
