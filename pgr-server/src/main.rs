use axum::{
    extract::{Extension, Path, Query},
    http::{
        header::{self, HeaderMap, HeaderName},
        HeaderValue, Method,
    },
    response::{Html, IntoResponse},
    routing::{get, post},
    Json, Router,
};
use pgr_db::ext::*;
use pgr_db::{
    aln::{self, HitPair},
    fasta_io::reverse_complement,
};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::net::SocketAddr;
use std::sync::Arc;
use tower::ServiceBuilder;
use tower_http::cors::Any;
use tower_http::cors::CorsLayer;
use tower_http::trace::TraceLayer;
use tracing::Span;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

#[derive(Clone, Debug, Deserialize)]
pub struct ShmmrSpec {
    pub w: u32,
    pub k: u32,
    pub r: u32,
    pub min_span: u32,
    pub sketch: bool,
}

#[derive(Deserialize)]

struct SequenceQuerySpec {
    source: String,
    ctg: String,
    bgn: usize,
    end: usize,
    padding: usize,
    merge_range_tol: usize,
    pb_shmmr_spec: ShmmrSpec,
    min_cov: usize,
    min_branch_size: usize,
    bundle_length_cutoff: usize,
    bundle_merge_distance: usize,
}

#[derive(Serialize)]
struct TargetRanges {
    query_src_ctg: (String, String),
    matches: Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>,
    sid_ctg_src: Vec<(u32, String, String)>,
    //principal_bundle_decomposition: Vec<(u32, Vec<SmpsWithBundleLabel>)>,
}

type MatchSummary = (u32, u32, u32, u32, usize, bool); //(q_bgn, q_end, t_bgn, t_end, num_hits, reversed)
type PrincipalBundleBedRecord = (String, u32, u32, usize, usize, u32, usize, usize, String);

#[derive(Serialize)]
struct TargetRangesSimplified {
    query_src_ctg: (String, String),
    match_summary: Vec<(u32, Vec<MatchSummary>)>, // (q_id, vec[(q_bgn, q_end, t_bgn, t_end, num_hits, reversed)])
    sid_ctg_src: Vec<(u32, String, String)>,
    bundle_bed_records: Vec<PrincipalBundleBedRecord>,
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

#[tokio::main]
async fn main() {
    tracing_subscriber::registry()
        .with(tracing_subscriber::EnvFilter::new(
            std::env::var("RUST_LOG")
                .unwrap_or_else(|_| "example_tracing_aka_logging=debug,tower_http=debug".into()),
        ))
        .with(tracing_subscriber::fmt::layer())
        .init();

    let mut seq_db = SeqIndexDB::new();
    let _ = seq_db.load_from_agc_index(
        //"/wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-small_panel".to_string(),
        "/wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-v0".to_string(),
    );
    let seq_db = Arc::new(seq_db);
    // build our application with a route
    let app = Router::new()
        .route(
            "/",
            get({
                let seq_db = seq_db.clone();
                move || handler(seq_db)
            }),
        )
        .route(
            "/query_sdb",
            post({
                let seq_db = seq_db.clone();
                move |params| query_sdb_with(params, seq_db)
            }),
        )
        .layer(
            CorsLayer::new()
                .allow_origin(Any)
                //.allow_origin("http://127.0.0.1:8080".parse::<HeaderValue>().unwrap())
                .allow_methods(Any)
                .allow_headers(Any),
        )
        .layer(ServiceBuilder::new().layer(TraceLayer::new_for_http()));

    // run it
    let addr = SocketAddr::from(([127, 0, 0, 1], 3000));
    println!("listening on {}", addr);
    axum::Server::bind(&addr)
        .serve(app.into_make_service())
        .await
        .unwrap();
}

/*
async fn handler(seq_db: Arc<SeqIndexDB>) -> impl IntoResponse {
    let n_ctg = 0;
    let mut headers = HeaderMap::new();
    headers.insert(header::CONTENT_TYPE, "text/plain".parse().unwrap());
    headers.insert(header::URI, "http://127.0.0.1:3000".parse().unwrap());
    let rtn = format!("Hello, World! {}", n_ctg);
    (headers, rtn)
}
*/

async fn handler(seq_db: Arc<SeqIndexDB>) -> Json<usize> {
    let n_ctg = seq_db.seq_index.as_ref().unwrap().len();
    Json(n_ctg)
}

async fn query_sdb_with(
    Json(payload): Json<Option<SequenceQuerySpec>>,
    seq_db: Arc<SeqIndexDB>,
) -> Json<Option<TargetRangesSimplified>> {
    if payload.is_none() {
        return Json(None);
    };

    let payload = payload.unwrap();

    let sample_name = payload.source;
    let ctg_name = payload.ctg;
    let padding = payload.padding;
    let merge_range_tol = payload.merge_range_tol;
    let seq_len = match seq_db
        .seq_index
        .as_ref()
        .unwrap()
        .get(&(ctg_name.clone(), Some(sample_name.clone())))
    {
        None => 0,
        Some(value) => value.1,
    };

    let q_seq_len = payload.end - payload.bgn;
    let q_seq_bgn = if padding > payload.bgn {
        0
    } else {
        payload.bgn - padding
    };
    let q_seq_end = if payload.end + padding > seq_len as usize {
        seq_len as usize
    } else {
        payload.end + padding
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
                let (ctg, src, _ctg_len) = seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
                rgns.iter().for_each(|(b, e, _, orientation, _)| {
                    sub_seq_range_for_fasta.push((sid, *b, *e, *orientation, ctg.clone()));
                });
                let hits = rgns
                    .into_iter()
                    .map(|(b, e, _, orientation, mut aln)| {
                        aln.sort();
                        let q_bgn = aln[0].0 .0;
                        let q_end = aln[aln.len() - 1].0 .1;
                        (q_bgn, q_end, b, e, aln.len(), orientation == 1)
                    })
                    .collect::<Vec<MatchSummary>>();
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
    let shmmr_spec = payload.pb_shmmr_spec;
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

    let (principal_bundles_with_id, vertex_to_bundle_id_direction_pos) =
        new_seq_db.get_principal_bundles_with_id(payload.min_cov, payload.min_branch_size, None);
    
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

    let bundle_bed_records = seq_info.iter().flat_map(|(sid, sdata)| {
        let (ctg, _src, _len) = sdata;
        let smps = sid_smps.get(&sid).unwrap();
        let smp_partitions = group_smps_by_principle_bundle_id(
            smps,
            payload.bundle_length_cutoff,
            payload.bundle_merge_distance,
        );
        let mut ctg_bundle_count = FxHashMap::<usize, usize>::default();
        smp_partitions.iter().for_each(|p| {
            let bid = p[0].1;
            *ctg_bundle_count.entry(bid).or_insert_with(|| 0) += 1;
        });
        let pbundle_bed_records = smp_partitions.into_iter().map(|p| {
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

            (
                ctg.clone(),
                b,
                e,
                bid,
                bid_to_size[&bid],
                direction,
                p[0].3,
                p[p.len() - 1].3,
                is_repeat.to_string(),
            )
        }).collect::<Vec<PrincipalBundleBedRecord>>();
        pbundle_bed_records
    }).collect::<Vec<PrincipalBundleBedRecord>>();

    Json(Some(TargetRangesSimplified {
        query_src_ctg: (sample_name, ctg_name),
        match_summary,
        sid_ctg_src,
        bundle_bed_records
    }))
}
