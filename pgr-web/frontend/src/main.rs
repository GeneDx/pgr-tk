// main.rs
use log;
use wasm_logger;

use dioxus::prelude::*;
use reqwest;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_json;
use std::{collections::HashMap, rc::Rc};
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
    pub query_src_ctg: (String, String),
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
}

fn app(cx: Scope) -> Element {
    wasm_logger::init(wasm_logger::Config::default());

    
    let roi_json = include_str!("data/ROIs.json");
    let rois: FxHashMap<String, SequenceQuerySpec> = serde_json::from_str(roi_json).unwrap();

    let query = use_state(cx, || <Option<SequenceQuerySpec>>::None);
    let query_name = use_state(cx, || <Option<String>>::None);
    let targets = use_state(cx, || <Option<TargetMatchPrincipalBundles>>::None);
    let query_state = use_state(cx, || "Please send a query".to_string());


    let get_targets = move |_| {
        cx.spawn({
    
            log::debug!("query");
            
            let query_name = query_name.to_owned();
            let query = query.to_owned();
            let targets= targets.to_owned();
            let rois2: FxHashMap<String, SequenceQuerySpec> =
            serde_json::from_str(roi_json).unwrap();
           
            let query_state = query_state.to_owned();
            query_state.set("Fetching query results".to_string());
            async move {
                console::log_1(&"clicked".into()); 
                let window = web_sys::window().expect("global window does not exists");    
                let document = window.document().expect("expecting a document on window");
                
                let roi_selector: web_sys::HtmlSelectElement = document.get_element_by_id(&"ROI_selector").unwrap().dyn_into().unwrap();
                let options =  roi_selector.options();
                let selected_value = options.get_with_index(options.selected_index().unwrap() as u32).unwrap().get_attribute("value").unwrap();
                let new_query =rois2.get(&selected_value).unwrap().clone(); 
                query.modify(move |_| Some(new_query));
                query_name.modify(move |_| Some(selected_value.clone()));

                let client = reqwest::Client::new();
                let qn: Option<String> = query_name.current().as_ref().clone();
                let q = if qn.is_none() {
                    None
                } else {
                    Some(rois2.get(&qn.unwrap()).unwrap())
                };
                let url = base_url() + "/api/post_query_for_json_data";
                let response = client
                    .post(url)
                    .json(&q)
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
    };

    let mut kvs = rois
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect::<Vec<_>>();
    kvs.sort_by_key(|v| v.0.clone());
    cx.render(
        rsx! {
            div { class: "flex flex-row p-4", 
                div { class: "basis-2/4",
                    h2 {"PanGenome Research Tool Kit: Principal Bundle Decomposition Demo"}
                }
                div { class: "basis-1/4 mb-3 xl:w-96",
                 
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
                }

                div { class: "basis-1/4",
                    button { 
                        id: "query_button",
                        //disabled: "false",
                        class: "inline-block px-6 py-1.5 bg-blue-600 text-white rounded",
                        onclick: get_targets ,
                        
                        "Query" 
                    }
                }
            }

            div { id: "query_status",     
                cx.render(
                    rsx! {
                        div { 
                        class: "p-4",
                        "{query_state}"
                        }
                    }
                ) 
            }

            div { id: "query_results",
                query_results(cx, query.clone(), targets.clone())
    
            }
        }
    )
}


pub fn query_results(
    cx: Scope,
    query: UseState<Option<SequenceQuerySpec>>,
    targets: UseState<Option<TargetMatchPrincipalBundles>>
) -> Element {

    let query = query.to_owned();
    if query.is_none() {
        log::debug!("query none");
        let r = rsx! { div { class: "p-4", id: "query_results_title" } };
        return cx.render(r);
    }

    let targets = targets.current().as_ref().clone();
    if targets.is_none() {
        log::debug!("target none");
        let r = rsx! { div { class: "p-4", id: "query_results_title" } };
        return cx.render(r);
    }

    let val = targets.unwrap();

    let sid_to_ctg_src = val
        .sid_ctg_src
        .iter()
        .map(|v| {
            let (sid, ctg_name, src) = v;
            (*sid, (ctg_name, src))
        })
        .collect::<HashMap<u32, (&String, &String)>>();

    let query = query.as_ref().unwrap().clone();
    let qstr = serde_qs::to_string(&query.clone()).unwrap();
    let ctg = query.ctg.clone();
    let bgn = query.bgn;
    let end = query.end;
    let query_url = base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..];
    let query_url = use_state(cx, move || {query_url});
    log::debug!("query: {:?}", query);
    let query_spec = use_state(cx, || query.clone()); 

    log::debug!("query_spec: {:?}", query_spec);
    log::debug!("rendering query_results");
    cx.render (
    rsx!{
        div { class: "grid p-2  grid-cols-1 justify-center space-y-2",
            div { class: "overflow-x-auto sm:-mx-6 lg:-mx-8",
                div {class: "flex flex-col min-w-[1280px]  max-h-screen",
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
                    hr {class: "my-2 h-px bg-gray-700 border-0 dark:bg-gray-700"}
                    div {class: "px-8 py-1",
                        div {class: "flex-grow overflow-auto max-h-[250px]",
                            table { class: "relative w-full",
                                thead {
                                    tr{
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "sid"} 
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "contig"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "source"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "hit count"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "query span"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "query len"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "target span"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "target len"}
                                    }
                                }
                                tbody {
                                    class: "divide-y",
                                    rsx!(val.match_summary.iter().map(|v| {
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

            div {
                span { "w:" }
                input {
                    value: "{query_spec.w}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.w = val.parse::<u32>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "k:" }
                input {
                    value: "{query_spec.k}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.k = val.parse::<u32>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "r:" }
                input {
                    value: "{query_spec.r}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.r = val.parse::<u32>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "min_span:" }
                input {
                    value: "{query_spec.min_span}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.min_span = val.parse::<u32>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "min_cov:" }
                input {
                    value: "{query_spec.min_cov}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.min_cov = val.parse::<usize>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "min_branch_size:" }
                input {
                    value: "{query_spec.min_branch_size}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.min_branch_size = val.parse::<usize>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "bundle_length_cutoff:" }
                input {
                    value: "{query_spec.bundle_length_cutoff}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.bundle_length_cutoff = val.parse::<usize>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }

            div {
                span { "bundle_merge_distance:" }
                input {
                    value: "{query_spec.bundle_merge_distance}",
                    oninput: move |evt| query_spec.modify(|v| {let val = evt.value.clone();
                        let mut v = v.clone();
                        v.bundle_merge_distance = val.parse::<usize>().expect("not an integer");
                        let qstr = serde_qs::to_string(&v.clone()).unwrap();
                        query_url.set(base_url() + "/api/get_html_by_query/?" + &qstr.clone()[..]);
                        v
                    }),
                }
            }


            div { class: "basis-1/4",
                button { 
                    id: "get_html_button",
                    class: "inline-block px-6 py-1.5 bg-blue-600 text-white rounded",
                    a { href: "{query_url}", target:"_blank", "Get Principal Bundle Decomposition HTML" }
                
                }
            } 
        }
    )
}
