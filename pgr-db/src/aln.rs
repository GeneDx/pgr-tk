// use rayon::prelude::*;
use crate::seq_db::{self, FragmentHit};
use crate::shmmrutils::{self, ShmmrSpec};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::HashSet;

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

    let mut target_squence_id_to_hits = FxHashMap::<u32, Vec<((u32, u32, u8), (u32, u32, u8))>>::default();
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
            frag_signature.iter().for_each(|&(_frg_id, sid, pos0, pos1, orientation)| {
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

#[test]

fn sparse_aln_test() {
    use crate::aln::sparse_aln;
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
