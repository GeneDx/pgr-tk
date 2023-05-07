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
use pgr_server::bundle_processing::*;
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
        .route(
            "/get_html_by_query_with",
            get({
                let seq_db = seq_db.clone();
                move |params| get_html_by_query_with(params, seq_db)
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
    let addr = SocketAddr::from(([0, 0, 0, 0], 3000));
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
    Json(seq_query_spec): Json<Option<SequenceQuerySpec>>,
    seq_db: Arc<SeqIndexDB>,
) -> Json<Option<TargetMatchPrincipalBundles>> {
    if seq_query_spec.is_none() {
        return Json(None);
    };

    let seq_query_spec = seq_query_spec.unwrap();
    println!("{:?}", seq_query_spec);
    Json(get_target_and_principal_bundle_decomposition(
        &seq_query_spec, seq_db,
    ))
}


async fn get_html_by_query_with(
    Json(seq_query_spec): Json<Option<SequenceQuerySpec>>,
    seq_db: Arc<SeqIndexDB>,
) -> Html<String> {
    if seq_query_spec.is_none() {
        return Html("".into());
    };

    let seq_query_spec = seq_query_spec.unwrap();
    println!("{:?}", seq_query_spec);

    let data  = get_target_and_principal_bundle_decomposition(
        &seq_query_spec, seq_db);
    let output = pb_data_to_html_string(&data.unwrap()) ;
    
    Html(output)

}
