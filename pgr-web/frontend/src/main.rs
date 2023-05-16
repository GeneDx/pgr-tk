// main.rs
use log;
use wasm_logger;

use dioxus::prelude::*;
use reqwest;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;

//use pgr_db::aln::{self, HitPair};
//type HitPair = ((u32, u32, u8), (u32, u32, u8)); //(bgn1, end1, orientation1),  (bgn2, end2, orientation2)

//type SmpBundleTuple = ((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>);
//type SmpsWithBundleLabel = Vec<SmpBundleTuple>;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MatchSummary {
    pub q_bgn: u32,
    pub q_end: u32,
    pub t_bgn: u32,
    pub t_end: u32,
    pub num_hits: usize,
    pub reversed: bool,
}

#[derive(Deserialize, Clone, Debug)]

pub struct TargetMatchPrincipalBundles {
    pub query: SequenceQuerySpec,
    pub match_summary: Vec<(u32, Vec<MatchSummary>)>, // (q_id, vec[(q_bgn, q_end, t_bgn, t_end, num_hits, reversed)])
    pub sid_ctg_src: Vec<(u32, String, String)>,
    pub bundle_bed_records: Vec<Vec<PrincipalBundleBedRecord>>,
}

#[derive(Deserialize, Clone, Debug)]
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

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ShmmrSpec {
    pub w: u32,
    pub k: u32,
    pub r: u32,
    pub min_span: u32,
    pub sketch: bool,
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct SequenceQuerySpec {
    pub source: String,
    pub ctg: String,
    pub bgn: usize,
    pub end: usize,
    pub padding: usize,
    pub merge_range_tol: usize,
    //pub pb_shmmr_spec: ShmmrSpec,
    // flatten this out, make it easier for URL query string
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

#[derive(Clone)]
struct QueryState(String);

pub fn base_url() -> String {
    web_sys::window().unwrap().location().origin().unwrap()
}

fn main() {
    dioxus_web::launch(app);
    wasm_logger::init(wasm_logger::Config::default());
}

fn app(cx: Scope) -> Element {

    let roi_json = include_str!("data/ROIs.json");
    let rois: FxHashMap<String, SequenceQuerySpec> = serde_json::from_str(roi_json).unwrap();

    let query = use_state(cx, || SequenceQuerySpec {
        source: "hg19_tagged.fa".to_string(),
		ctg: "chr1_hg19".to_string(),
		bgn: 104198140,
		end: 104207173,
		padding: 150000,
		merge_range_tol: 120000,
		w: 48,
		k: 56,
		r: 4,
		min_span: 12,
		sketch: false,
		min_cov: 2,
		min_branch_size: 8,
		bundle_length_cutoff: 500,
		bundle_merge_distance: 10000
    });
    let targets = use_state(cx, || <Option<TargetMatchPrincipalBundles>>::None);
    let query_state = use_state(cx, || "Please send a query".to_string());  

    let mut kvs = rois
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect::<Vec<_>>();

    kvs.sort_by_key(|v| v.0.clone());
    let labels = kvs.iter().map(|(k, _v)| k.clone()).collect::<Vec<_>>(); 

    let selected_label = use_state(cx, || labels.get(0).unwrap_or(&"".to_string()).clone());

   
    
    cx.render(
        rsx! {
            div { class: "flex justify-center p-5",
                p { class: "text-2xl", "PanGenome Research Tool Kit: Principal Bundle Decomposition Demo" }
            }
            div { class: "container justify-center mx-auto w-full p-5",

                div { class: "flex flex-row",

                    div { class: "basis-4/6", id: "query_results", query_results { targets: targets } }

                    div { class: "basis-2/6",
                        div { class: "p-1", id: "query_status", "Status: {query_state}" }
                        div { class: "p-1", query_preset { labels: labels, selected_label: selected_label } }
                        button {
                            class: "p-1",
                            id: "query_button",
                            //disabled: "false",
                            class: "middle none center w-full rounded-lg px-2 py-1.5 bg-blue-600 text-white",
                            onclick: move |_evt| {
                                let query_name = selected_label.current().as_ref().clone();
                                let query0 = rois.get(&query_name).unwrap();
                                query.set(query0.clone());
                            },
                            "Set Query Parameters"
                        }
                        div { class: "p-1", id: "set_parameters", set_parameters { query: query } }
                        div { class: "flex flex-row p-1",
                            div { class: "basis-1/2 p-1",
                                update_query { query: query, targets: targets, query_state: query_state }
                            }
                            br {}
                            div { class: "basis-1/2 p-1", id: "get_html", get_html { query: query } }
                        }
                    }
                }
            }
        }
    )
}

#[inline_props]
fn query_preset<'a>(cx: Scope<'a>, labels: Vec<String>, selected_label: &'a UseState<String>) -> Element<'a> {
    cx.render(
        rsx! {
            div { class: "flex flex-row p-0",
                div { class: "basis-2/4", "Query Preset:" }
                select {
                    class: "basis-2/4",
                    name: "ROI_selector",
                    id: "ROI_selector",
                    class: "form-select appearance-none  w-full px-3 py-1.5 focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none",
                    oninput: |evt| {
                        selected_label.set(evt.value.clone());
                    },
                    labels.iter().map(|k| {
                        rsx! { 
                            option {
                                value: "{k}",
                                "{k}"
                            }
                        }
                    })
                }
            }
        }) 
}

fn get_targets<'a, T>(cx: Scope<'a, T>, 
    query: &'a SequenceQuerySpec,
    targets: &'a UseState<Option<TargetMatchPrincipalBundles>>, 
    query_state: &'a UseState<String>) -> (){
    let query = query.clone();
    let targets = targets.to_owned(); 
    let query_state = query_state.to_owned(); 
    
    cx.spawn({
            
        async move {
            let client = reqwest::Client::new();
            let url = base_url() + "/api/post_query_for_json_data";
            let response = client
                .post(url)
                .json(&query)
                .send()
                .await
                .unwrap()
                .json::<Option<TargetMatchPrincipalBundles>>()
                .await;
            match response {
                Ok(val) => {
                    targets.set(val);
                    query_state.set("Query results fetched".into());
                },
                Err(e) => {
                    log::debug!("{:?}",e);
                }
            };
        }
    })
}

#[inline_props]
pub fn query_results<'a>(
    cx: Scope<'a>,
    targets: &'a UseState<Option<TargetMatchPrincipalBundles>>
) -> Element<'a> {

    let targets = targets.current().as_ref().clone();
    if targets.is_none() {
        log::debug!("target none");
        let r = rsx! { div { class: "p-4", id: "query_results_title" } };
        return cx.render(r);
    }

    let targets = targets.unwrap();

    let sid_to_ctg_src = targets
        .sid_ctg_src
        .iter()
        .map(|v| {
            let (sid, ctg_name, src) = v;
            (*sid, (ctg_name, src))
        })
        .collect::<HashMap<u32, (&String, &String)>>();

    let query = targets.query.clone();
    let ctg = query.ctg;
    let bgn = query.bgn;
    let end = query.end;

    cx.render ( {
        
    rsx!{
        div { class: "grid p-2  grid-cols-1 justify-center space-y-2",
            div { class: "flex flex-col min-w-[1280px]  max-h-screen",
                h2 { class: "px-8 py-2", p { "Returned Hits for Query: {ctg}:{bgn}-{end}" } }
                div { class: "px-8 content-center overflow-auto min-w-[1280px] max-h-[60px]" }
            }
            //hr { class: "my-2 h-px bg-gray-700 border-0 dark:bg-gray-700" }
            div { class: "flex flex-col px-8 py-1",
                div { class: "flex-grow overflow-auto max-h-[650px]",
                    table { class: "relative w-full",
                        thead {
                            tr {
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "sid"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "contig"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "source"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "hit count"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "query span"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "query len"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "target span"
                                }
                                th { class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300",
                                    "target len"
                                }
                            }
                        }
                        tbody { class: "divide-y",
                            rsx!(targets.match_summary.iter().map(|v| {
                                        let sid = v.0;
                                        let (ctg, src) = *sid_to_ctg_src.get(&sid).unwrap();
                                        let style_classes = "px-1 py-2 text-center";
                                        let hit_summary = v.1.iter().map(move |match_summary| {
                                            let ms = match_summary;
                                            let q_span = format!("{}-{}", ms.q_bgn, ms.q_end);
                                            let t_span = format!("{}-{}", ms.t_bgn, ms.t_end);
                                            let q_len = if ms.q_end > ms.q_bgn { ms.q_end - ms.q_bgn } else { ms.q_bgn - ms.q_end };
                                            let t_len = if ms.t_end > ms.t_bgn { ms.t_end - ms.t_bgn } else { ms.t_bgn - ms.t_end };
                                            let n_hits = ms.num_hits;
                                            rsx!( tr {
                                                td { class: "{style_classes}", "{sid}"}  
                                                td { class: "{style_classes}", "{ctg}"} 
                                                td { class: "{style_classes}", "{src}"}
                                                td { class: "{style_classes}", "{n_hits}"} 
                                                td { class: "{style_classes}", "{q_span}"} 
                                                td { class: "{style_classes}", "{q_len}"} 
                                                td { class: "{style_classes}", "{t_span}"}
                                                td { class: "{style_classes}", "{t_len}"}
                                                } )
                                        });
                                            
                                    rsx!( hit_summary )
                                    }))
                        }
                    }
                }
            }
        }
    }}
    )
}

macro_rules! set_parameter {
    ($fn_name:ident, $field: ident, $type: ty) => {
        #[inline_props]
        fn $fn_name<'a>(cx: Scope<'a>, query: &'a UseState<SequenceQuerySpec>)  -> Element<'a> {
            let val = query.$field.clone(); 
            cx.render ( 
                rsx! { td {
                    input {
                        value: "{val}",
                        oninput: move |evt| {
                            let val = evt.value.clone().parse::<$type>();
                            if let Ok(val) = val {
                                let mut new_query = (*query.get()).clone();
                                new_query.$field = val;
                                query.set(new_query);
                            }
                        }
                    }
                }}
            )
        } 
    };
}

set_parameter!(set_parameter_source, source, String);
set_parameter!(set_parameter_ctg, ctg, String);
set_parameter!(set_parameter_bgn, bgn, usize);
set_parameter!(set_parameter_end, end, usize);
set_parameter!(set_parameter_padding, padding, usize);
//set_parameter!(set_parameter_merge_range_tol, merge_range_tol, usize);

set_parameter!(set_parameter_w, w, u32);
set_parameter!(set_parameter_k, k, u32);
set_parameter!(set_parameter_r, r, u32);
set_parameter!(set_parameter_min_span, min_span, u32);

set_parameter!(set_parameter_min_cov, min_cov, usize);
set_parameter!(set_parameter_min_branch_size, min_branch_size, usize);
set_parameter!(set_parameter_bundle_length_cutoff, bundle_length_cutoff, usize);
set_parameter!(set_parameter_bundle_merge_distance, bundle_merge_distance, usize);

#[inline_props]
fn set_parameters<'a>(cx: Scope<'a>, query:&'a UseState<SequenceQuerySpec> ) -> Element<'a> {
 
    cx.render ( 
        rsx!{
            div {
                table {
                    thead {
                        tr {
                            th { class: "px-5 py-2", "parameter" }
                            th { class: "px-5 py-2", "value" }
                        }
                    }
                    tbody {

                        tr {
                            td { class: "px-5 py-1", "source" }
                            set_parameter_source { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "query_ctg" }
                            set_parameter_ctg { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "begin coordinate" }
                            set_parameter_bgn { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "end coordinate" }
                            set_parameter_end { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "flanking size" }
                            set_parameter_padding { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "w" }
                            set_parameter_w { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "k" }
                            set_parameter_k { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "r" }
                            set_parameter_r { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "min span" }
                            set_parameter_min_span { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "min cov" }
                            set_parameter_min_cov { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "min branch size" }
                            set_parameter_min_branch_size { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "bundle length cutoff" }
                            set_parameter_bundle_length_cutoff { query: query }
                        }

                        tr {
                            td { class: "px-5 py-1", "bundle merge distance" }
                            set_parameter_bundle_merge_distance { query: query }
                        }
                    }
                }
            }
        })
} 


#[inline_props]
pub fn get_html<'a>( cx: Scope<'a>, query: &'a UseState<SequenceQuerySpec> ) -> Element<'a> {

    let query = query.current().as_ref().clone();
    let query_url = {
        let qstr = serde_qs::to_string(&query).unwrap();
        base_url() + "/api/get_html_by_query/?" + &qstr[..]
    };

    cx.render ( {
        rsx!{
            button {
                id: "get_html_button",
                class: "middle none center w-full rounded-lg px-2 py-1.5 bg-blue-600 text-white",
                
                a { class: "w-full", href: "{query_url}", target: "_blank", p { "Get HTML" } }
            }
        }}
    )
}

#[inline_props]
pub fn update_query<'a>( cx: Scope<'a>, 
    query: &'a UseState<SequenceQuerySpec>,
    targets: &'a UseState<Option<TargetMatchPrincipalBundles>> ,
    query_state: &'a UseState<String>) -> Element<'a> {

    let query = query.to_owned();
    let targets = targets.to_owned();
    let query_state = query_state.to_owned(); 

    cx.render ( {
        rsx!{
            button {
                id: "query_button",
                class: "middle none center w-full rounded-lg px-2 py-1.5 bg-blue-600 text-white",
                //disabled: "false",
                onclick: move |_evt| {
                    let query0 = query.get();
                    get_targets(cx, query0, &targets, &query_state);
                    query_state.set("getting query results".to_string());
                },
                "Update"
            }
        }}
    )
}



