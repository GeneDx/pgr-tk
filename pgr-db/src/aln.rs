// use rayon::prelude::*;
use crate::seq_db::{self, FragmentHit};
use crate::shmmrutils::{self, ShmmrSpec};
use log::debug;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::HashSet;
use wavefront_aln::*;

pub type HitPair = ((u32, u32, u8), (u32, u32, u8)); //(bgn1, end1, orientation1),  (bgn2, end2, orientation2)

pub fn sparse_aln(
    sp_hits: &mut Vec<HitPair>,
    max_span: u32,
    penalty: f32,
) -> Vec<(f32, Vec<HitPair>)> {
    // given a set of hits in the form of (bgn1, end1, orientation1),  (bgn2, end2, orientation2)
    // perform (banded) dynamic programming to group them into list of hit chains
    sp_hits.sort_by(|a, b| a.0 .0.partial_cmp(&b.0 .0).unwrap());
    let mut v_s = FxHashMap::<HitPair, f32>::default(); // score for each vertex
    let mut best_pre_v = FxHashMap::<HitPair, Option<HitPair>>::default(); // look up for the best pre-vertex
    assert!(sp_hits.len() > 1);
    let first_hp = sp_hits[0];
    v_s.insert(first_hp, first_hp.0 .1 as f32 - first_hp.0 .0 as f32); // the score of the first node is just its length
    best_pre_v.insert(first_hp, None);

    (1..sp_hits.len()).for_each(|i| {
        let hp = sp_hits[i];
        let mut best_v = Option::<HitPair>::None;
        let mut best_s = 0_f32;
        let mut j = i;
        let mut span_set = HashSet::<(u32, u32, u8)>::new();
        loop {
            if j == 0 {
                break;
            };
            j -= 1;

            let pre_hp = sp_hits[j];
            if pre_hp.0 == hp.0 {
                continue;
            }; // don't connect node with the same left coordinate
            span_set.insert(pre_hp.0);
            let p_s = v_s.get(&pre_hp).unwrap_or(&0_f32);
            let mut s: f32 = *p_s + (hp.0 .1 as f32 - hp.0 .0 as f32);

            if hp.0 .2 == hp.1 .2 {
                // same orientation
                s -= penalty
                    * ((hp.0 .0 as f32 - pre_hp.0 .1 as f32).abs()
                        + (hp.1 .0 as f32 - pre_hp.1 .1 as f32).abs());
            } else {
                // opposite orientation
                s -= penalty
                    * ((hp.0 .0 as f32 - pre_hp.0 .1 as f32).abs()
                        + (hp.1 .1 as f32 - pre_hp.1 .0 as f32).abs());
            }

            if s > best_s {
                best_s = s;
                best_v = Some(pre_hp);
            }

            if span_set.len() >= max_span as usize {
                break;
            };
        }

        if best_s > 0_f32 {
            v_s.insert(hp, best_s);
            best_pre_v.insert(hp, best_v);
        } else {
            v_s.insert(hp, hp.0 .1 as f32 - hp.0 .0 as f32);
            best_pre_v.insert(hp, None);
        }
    });

    let mut unvisited_v = FxHashSet::<HitPair>::default();
    unvisited_v.extend(sp_hits.iter());
    let mut out = Vec::<(f32, Vec<HitPair>)>::new();
    while !unvisited_v.is_empty() {
        let mut best_s = 0_f32; // global best score
        let mut best_v: Option<HitPair> = None; // global best vertex
                                                // println!("DBG un-visit len; {}", unvisited_v.len());
        unvisited_v.iter().for_each(|hp| {
            let s = v_s.get(hp).unwrap_or(&0_f32);
            if *s > best_s {
                best_s = *s;
                best_v = Some(*hp);
            }
        });
        let mut track = Vec::<HitPair>::new();
        let mut v = best_v;
        while v.is_some() {
            let hp = v.unwrap();
            if !unvisited_v.contains(&hp) {
                break;
            };
            track.push(hp);
            v = *best_pre_v.get(&hp).unwrap_or(&None);
        }
        if track.is_empty() {
            continue;
        };
        track.reverse();
        track.iter().for_each(|hp| {
            // let s = v_s.get(hp).unwrap_or(&0_f32);
            // println!("H {} {} {} {} {} {} {}", hp.0.0, hp.0.1, hp.0.2, hp.1.0, hp.1.1, hp.1.2, s );
            unvisited_v.remove(hp);
        });
        let bgn_s = v_s.get(&track[0]).unwrap_or(&0_f32);
        out.push((best_s - bgn_s, track));
    }
    out
}

pub type TargetHitPairLists = Vec<(u32, Vec<(f32, Vec<HitPair>)>)>; // target_id, Vec<(score, HitPairs)>

#[allow(clippy::too_many_arguments)]
pub fn query_fragment_to_hps(
    raw_query_hits: Vec<FragmentHit>,
    frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
    penalty: f32,
    max_count: Option<u32>,
    query_max_count: Option<u32>,
    target_max_count: Option<u32>,
    max_aln_span: Option<u32>,
) -> TargetHitPairLists {
    let mut shmmr_pair_hash_count = FxHashMap::<(u64, u64), u32>::default();
    let mut query_shmmr_pair_hash_count = FxHashMap::<(u64, u64), u32>::default();
    let mut target_shmer_pair_count = FxHashMap::<(u64, u64, u32), u32>::default();

    seq_db::pair_shmmrs(&shmmrutils::sequence_to_shmmrs(0, frag, shmmr_spec, false))
        .iter()
        .for_each(|shmmr_pair| {
            let entry = query_shmmr_pair_hash_count
                .entry((shmmr_pair.0.hash(), shmmr_pair.1.hash()))
                .or_insert(0);
            *entry += 1;
        });

    raw_query_hits.iter().for_each(
        |(shmmr_pair_hash, _query_position, frag_signature): &(
            (u64, u64),
            _,
            Vec<seq_db::FragmentSignature>,
        )| {
            //let sp = d.0;
            // count shimmer pair hits
            let entry = shmmr_pair_hash_count.entry(*shmmr_pair_hash).or_insert(0);
            *entry += 1;

            frag_signature
                .iter()
                .for_each(|(_frg_id, seq_id, _bgn, _end, _orientation)| {
                    // count shimmer pair on target hits
                    // v = frg_id, seq_id, bgn, end, orientation(to shimmer pair)
                    let key = (shmmr_pair_hash.0, shmmr_pair_hash.1, *seq_id);
                    let entry = target_shmer_pair_count.entry(key).or_insert(0);
                    *entry += 1;
                })
        },
    );

    let mut target_squence_id_to_hits =
        FxHashMap::<u32, Vec<((u32, u32, u8), (u32, u32, u8))>>::default();
    raw_query_hits.into_iter().for_each(
        |(shmmr_pair, query_position, frag_signature): (
            (u64, u64),
            _,
            Vec<seq_db::FragmentSignature>,
        )| {
            let count = *shmmr_pair_hash_count.get(&shmmr_pair).unwrap_or(&0);
            let max_count = max_count.unwrap_or(128);
            if count > max_count {
                return;
            };
            let max_count_query = query_max_count.unwrap_or(128);
            if count > max_count_query {
                return;
            };
            let left_frag_coordinate = query_position;
            frag_signature
                .iter()
                .for_each(|&(_frg_id, sid, pos0, pos1, orientation)| {
                    let count = *target_shmer_pair_count
                        .get(&(shmmr_pair.0, shmmr_pair.1, sid))
                        .unwrap_or(&0);
                    let max_count_target = target_max_count.unwrap_or(128);
                    if count > max_count_target {
                        return;
                    };
                    let e = target_squence_id_to_hits.entry(sid).or_default();
                    let right_frag_coordinate = (pos0, pos1, orientation);
                    e.push((left_frag_coordinate, right_frag_coordinate));
                });
        },
    );

    let max_aln_span = max_aln_span.unwrap_or(8);

    target_squence_id_to_hits
        .into_iter()
        .filter(|(_sid, hps)| hps.len() > 1)
        .map(|(sid, mut hps)| (sid, sparse_aln(&mut hps, max_aln_span, penalty)))
        .collect::<Vec<_>>()
}

pub fn wfa_align_bases(
    target_str: &str,
    query_str: &str,
    max_wf_length: u32,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
) -> Option<(String, String)> {
    let capacity = std::cmp::max(1024, std::cmp::max(target_str.len(), query_str.len()) >> 5);
    let mut wfs = WaveFronts::new_with_capacity(
        target_str,
        query_str,
        max_wf_length,
        mismatch_penalty,
        open_penalty,
        extension_penalty,
        capacity,
    );
    if wfs.step_all(Some(1024)) == WaveFrontStepResult::ReachEnd {
        Some(wfs.backtrace())
    } else {
        None
    }
}

pub fn wfa_aln_pair_map(aln_target_str: &str, aln_query_str: &str) -> Vec<(u32, u32, char)> {
    let paired = std::iter::zip(aln_target_str.as_bytes(), aln_query_str.as_bytes());
    let mut t_pos = 0_u32;
    let mut q_pos = 0_u32;
    paired
        .into_iter()
        .map(|(&tb, &qb)| {
            let mut t = '-';
            let new_t_pos = if tb == b'-' {
                t = 'I';
                t_pos
            } else {
                t_pos + 1
            };
            let new_q_pos = if qb == b'-' {
                t = 'D';
                q_pos
            } else {
                q_pos + 1
            };
            if tb == qb {
                t = 'M';
            };
            if tb != qb && tb != b'-' && qb != b'-' {
                t = 'X';
            };
            let out = (t_pos, q_pos, t);
            t_pos = new_t_pos;
            q_pos = new_q_pos;
            out
        })
        .collect::<Vec<_>>()
}

pub fn get_variants_from_aln_pair_map(
    aln_pairs: &[(u32, u32, char)],
    target_str: &str,
    query_str: &str,
) -> Vec<(u32, u32, char, String, String)> {
    let mut current_variant = Vec::<(char, char, char)>::new();
    let aggregate_variants = |previous_match: &(u32, u32, char, char, char),
                              current_variant: &Vec<(char, char, char)>|
     -> Option<(u32, u32, char, String, String)> {
        let t_variant_segment = String::from_iter(current_variant.iter().map(|v| v.0));
        let q_variant_segment = String::from_iter(current_variant.iter().map(|v| v.1));
        let v_type = current_variant[0].2;
        match v_type {
            'X' => Some((
                previous_match.0 + 1,
                previous_match.1 + 1,
                'X',
                t_variant_segment,
                q_variant_segment,
            )),
            'I' => Some((
                previous_match.0,
                previous_match.1,
                'I',
                [previous_match.3.to_string(), t_variant_segment].join(""),
                [previous_match.4.to_string(), q_variant_segment].join(""),
            )),
            'D' => Some((
                previous_match.0,
                previous_match.1,
                'D',
                [previous_match.3.to_string(), t_variant_segment].join(""),
                [previous_match.4.to_string(), q_variant_segment].join(""),
            )),
            _ => None,
        }
    };
    let mut previous_match = (0_u32, 0_u32, 'U', '-', '-');
    let mut variants = Vec::<Option<(u32, u32, char, String, String)>>::new();

    aln_pairs.iter().for_each(|&(t_pos, q_pos, t)| match t {
        'M' => {
            let t_char = target_str.as_bytes()[t_pos as usize] as char;
            let q_char = query_str.as_bytes()[q_pos as usize] as char;
            if !current_variant.is_empty() {
                variants.push(aggregate_variants(&previous_match, &current_variant));
                current_variant.clear();
            };
            previous_match = (t_pos, q_pos, 'M', t_char, q_char);
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, t_char, q_char, t);
        }
        'X' => {
            let t_char = target_str.as_bytes()[t_pos as usize] as char;
            let q_char = query_str.as_bytes()[q_pos as usize] as char;
            if !current_variant.is_empty() {
                if current_variant.last().unwrap().2 == t {
                    current_variant.push((t_char, q_char, t));
                } else {
                    variants.push(aggregate_variants(&previous_match, &current_variant));
                    current_variant.clear();
                    current_variant.push((t_char, q_char, t));
                }
            } else {
                current_variant.push((t_char, q_char, t));
            }
        }
        'I' => {
            let q_char = query_str.as_bytes()[q_pos as usize] as char;
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, '-', q_char, t);
            if !current_variant.is_empty() {
                if current_variant.last().unwrap().2 == t {
                    current_variant.push(('-', q_char, t));
                } else {
                    variants.push(aggregate_variants(&previous_match, &current_variant));
                    current_variant.clear();
                    current_variant.push(('-', q_char, t));
                }
            } else {
                current_variant.push(('-', q_char, t));
            }
        }
        'D' => {
            let t_char = target_str.as_bytes()[t_pos as usize] as char;
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, t_char, '-', t);
            if !current_variant.is_empty() {
                if current_variant.last().unwrap().2 == t {
                    current_variant.push((t_char, '-', t));
                } else {
                    variants.push(aggregate_variants(&previous_match, &current_variant));
                    current_variant.clear();
                    current_variant.push((t_char, '-', t));
                }
            } else {
                current_variant.push((t_char, '-', t));
            }
        }
        _ => {}
    });
    if !current_variant.is_empty() {
        variants.push(aggregate_variants(&previous_match, &current_variant));
    };
    variants.into_iter().flatten().collect::<Vec<_>>()
}

#[cfg(test)]
mod test {

    #[test]

    fn sparse_aln_test() {
        use crate::aln::{sparse_aln, HitPair};
        use std::fs::File;
        use std::io::{BufRead, BufReader};
        let f = BufReader::new(File::open("./test/test_data/test_hits").unwrap());
        let mut hp = Vec::<HitPair>::new();
        f.lines().for_each(|s| {
            if let Ok(s) = s {
                let s = s.split_ascii_whitespace();
                let out = s
                    .into_iter()
                    .map(|s| s.parse::<u32>().unwrap())
                    .collect::<Vec<u32>>();
                assert_eq!(out.len(), 6);
                hp.push((
                    (out[0], out[1], out[2] as u8),
                    (out[3], out[4], out[5] as u8),
                ));
            }
        });

        let out = sparse_aln(&mut hp, 8, 0.5_f32);
        out.iter().for_each(|(s, v)| println!("{} {}", s, v.len()));
    }

    #[test]
    fn test_wfa_align_bases() {
        use crate::aln::{get_variants_from_aln_pair_map, wfa_align_bases, wfa_aln_pair_map};
        use log::debug;
        use simple_logger::SimpleLogger;
        SimpleLogger::new().init().unwrap();
        let t_str = "ACATACATGTGTGTGAAAAATATATAAGTAAAAAAAATGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCATGAAACCCCATACATGAAAGTTGCATGAA";
        let q_str = "ACATACATGTGAAATATAATAAAAGTTGCATGAAAAAACATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCAAAAGTTGCATGAAACATACATGAAAATGAAAAAACATACATGAAAGTTGCATGAA";
        if let Some((t_aln_str, q_aln_str)) = wfa_align_bases(t_str, q_str, 20, 2, 2, 1) {
            println!("{}", t_aln_str);
            println!("{}", q_aln_str);
            let aln_pairs = wfa_aln_pair_map(&t_aln_str, &q_aln_str);
            let variants = get_variants_from_aln_pair_map(&aln_pairs, t_str, q_str);
            variants.into_iter().for_each(|(t_pos, q_pos, t, s1, s2)| {
                println!("{} {} {} {} {}", t_pos, q_pos, t, s1, s2);
            });
        }
    }
}
