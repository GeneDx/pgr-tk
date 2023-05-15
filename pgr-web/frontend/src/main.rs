// main.rs
use log;
use wasm_logger;

use dioxus::prelude::*;
use reqwest;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use wasm_bindgen::JsCast;
use web_sys::console;

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

    let get_seleted_query_name = || {
        let window = web_sys::window().expect("global window does not exists");
        let document = window.document().expect("expecting a document on window");
        let roi_selector: web_sys::HtmlSelectElement = document
            .get_element_by_id(&"ROI_selector")
            .unwrap()
            .dyn_into()
            .unwrap();
        let options = roi_selector.options();
        let selected_value = options
            .get_with_index(options.selected_index().unwrap() as u32)
            .unwrap()
            .get_attribute("value")
            .unwrap();
        let query_name = selected_value.clone();
        query_name
    };

    cx.render(
        rsx! {
            div { class: "flex flex-row p-4",
                div { class: "basis-2/6", h2 { "PanGenome Research Tool Kit: Principal Bundle Decomposition Demo" } }
                div { class: "basis-3/6",
                    div { class: "flex flex-row p-4",
                        div { "Query Preset:" }
                        select {
                            name: "ROI_selector",
                            id: "ROI_selector",
                            class: "form-select appearance-none  w-full px-3 py-1.5 focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none",
                            kvs.iter().map(|(k, _v)| {
                            rsx! { 
                                option {
                                    value: "{k}",
                                    "{k}"
                                }
                            }
                        })
                        }
                        button {
                            class: "basis-1/5 mb-3 xl:w-32",
                            id: "query_button",
                            //disabled: "false",
                            class: "inline-block px-6 py-1.5 bg-blue-600 text-white rounded",
                            onclick: move |_evt| {
                                console::log_1(&"clicked".into());
                                let query_name = get_seleted_query_name();
                                let query0 = rois.get(&query_name).unwrap();
                                get_targets(cx, query0, targets, query_state);
                                query_state.set("getting query results".to_string());
                                query.set(query0.clone());
                            },
                            "Query"
                        }
                    }
                }
            }

            div { id: "query_status",
                cx.render(
                    rsx! {
                        div { 
                        class: "p-4",
                        "status: {query_state}"
                        }
                    }
                )
            }
            div { class: "flex flex-row p-4", id: "query_results", query_results(cx, targets.clone()) }
            div { class: "flex flex-row p-4",
                div { class: "basis-3/6" }
                div { class: "basis-2/6", id: "set_parameters", set_parameters(cx, query.clone()) }
                div { class: "basis-1/6",
                    div { update_query(cx, query.clone(),  targets.clone(),  query_state.clone()) }

                    br {}

                    div { id: "get_html", get_html(cx, query.clone()) }
                }
            }
        }
    )
}

 

fn get_targets(cx: Scope, 
    query: &SequenceQuerySpec,
    targets: &UseState<Option<TargetMatchPrincipalBundles>>, 
    query_state: &UseState<String>) -> (){
    let query = query.clone();
    let targets = targets.to_owned(); 
    let query_state = query_state.to_owned(); 
    
    cx.spawn({
            
        log::debug!("query");
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

pub fn query_results(
    cx: Scope,
    targets: UseState<Option<TargetMatchPrincipalBundles>>
) -> Element {

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
            div { class: "overflow-x-auto sm:-mx-6 lg:-mx-8",
                div { class: "flex flex-col min-w-[1280px]  max-h-screen",
                    rsx!(
                        h2 {class: "px-8 py-2", "Principal Bundle Decomposition, Query: {ctg}:{bgn}-{end}"}
                        div {
                            class: "px-8 content-center overflow-auto min-w-[1280px] max-h-[450px]",
                            //val.principal_bundle_decomposition.iter().flat_map(|(sid, ctg_name, r)| {
                            //    track(cx, ctg_name.clone(), track_size, (*sid, r.clone()))
                            //}
                        }
                       )
                }
                hr { class: "my-2 h-px bg-gray-700 border-0 dark:bg-gray-700" }
                div { class: "px-8 py-1",
                    div { class: "flex-grow overflow-auto max-h-[250px]",
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
                                            
                                    rsx!( hit_summary)
                                    }))
                            }
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
        fn $fn_name(cx: Scope, query: UseState<SequenceQuerySpec>)  -> Element {
            let query = query.to_owned();
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

fn set_parameters(cx: Scope, query: UseState<SequenceQuerySpec> ) -> Element {
 
    cx.render ( 
        rsx!{
            div {
                table {
                    thead {
                        tr {
                            th { "parameter" }
                            th { "value" }
                        }
                    }
                    tbody {

                        tr {
                            td { "source" }
                            set_parameter_source(cx, query.to_owned())
                        }

                        tr {
                            td { "query_ctg" }
                            set_parameter_ctg(cx, query.to_owned())
                        }

                        tr {
                            td { "begin coordinate" }
                            set_parameter_bgn(cx, query.to_owned())
                        }

                        tr {
                            td { "end coordinate" }
                            set_parameter_end(cx, query.to_owned())
                        }

                        tr {
                            td { "flanking size" }
                            set_parameter_padding(cx, query.to_owned())
                        }

                        tr {
                            td { "w" }
                            set_parameter_w(cx, query.to_owned())
                        }

                        tr {
                            td { "k" }
                            set_parameter_k(cx, query.to_owned())
                        }

                        tr {
                            td { "r" }
                            set_parameter_r(cx, query.to_owned())
                        }

                        tr {
                            td { "min span" }
                            set_parameter_min_span(cx, query.to_owned())
                        }

                        tr {
                            td { "min cov" }
                            set_parameter_min_cov(cx, query.to_owned())
                        }

                        tr {
                            td { "min branch size" }
                            set_parameter_min_branch_size(cx, query.to_owned())
                        }

                        tr {
                            td { "bundle length cutoff" }
                            set_parameter_bundle_length_cutoff(cx, query.to_owned())
                        }

                        tr {
                            td { "bundle merge distance" }
                            set_parameter_bundle_merge_distance(cx, query.to_owned())
                        }
                    }
                }
            }
        })
} 


pub fn get_html( cx: Scope, query: UseState<SequenceQuerySpec> ) -> Element {

    let query = query.current().as_ref().clone();
    let query_url = {
        let qstr = serde_qs::to_string(&query).unwrap();
        base_url() + "/api/get_html_by_query/?" + &qstr[..]
    };

    cx.render ( {
        rsx!{
            div { class: "basis-1/4",
                button {
                    id: "get_html_button",
                    class: "inline-block px-6 py-1.5 bg-blue-600 text-white rounded",
                    a { href: "{query_url}", target: "_blank", "Get HTML" }
                }
            }
        }}
    )
}


pub fn update_query( cx: Scope, 
    query: UseState<SequenceQuerySpec>,
    targets: UseState<Option<TargetMatchPrincipalBundles>> ,
    query_state: UseState<String>) -> Element {

    let query = query.to_owned();
    let targets = targets.to_owned();
    let query_state = query_state.to_owned(); 

    cx.render ( {
        rsx!{
            button {
                id: "query_button",
                //disabled: "false",
                class: "inline-block px-6 py-1.5 bg-blue-600 text-white rounded",
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



