use axum::{
    body::{boxed, Body},
    extract::Query,
    http::{Response, StatusCode},
    response::Html,
    routing::{get, post},
    Json, Router,
};
use clap::Parser;
use pgr_db::ext::*;
use pgr_server::bundle_processing::*;
use std::net::SocketAddr;
use std::{
    net::{IpAddr, Ipv6Addr},
    path::PathBuf,
    str::FromStr,
    sync::Arc,
};
use tokio::fs;
use tower::{ServiceBuilder, ServiceExt};
use tower_http::cors::Any;
use tower_http::cors::CorsLayer;
use tower_http::services::ServeDir;
use tower_http::trace::TraceLayer;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

#[derive(Parser, Debug)]
#[clap(name = "pgr-server", about = "Experimental Server")]
struct Opt {
    /// set the listen addr
    #[clap(short = 'a', long = "addr", default_value = "::1")]
    addr: String,

    /// set the listen port
    #[clap(short = 'p', long = "port", default_value = "5000")]
    port: u16,

    /// set the directory where static files are to be found
    #[clap(long = "static-dir", default_value = "./dist")]
    static_dir: String,

    /// set the listen port
    #[clap(
        short = 'd',
        long = "data-path-prefix",
        default_value = "./pgr-tk-HGRP-y1-evaluation-set-v0"
    )]
    data_path_prefix: String,

    /// set the listen port
    #[clap(short = 'f', long = "frg-file")]
    frg_file: bool,
}

#[tokio::main]
async fn main() {
    let opt = Opt::parse();

    tracing_subscriber::registry()
        .with(tracing_subscriber::EnvFilter::new(
            std::env::var("RUST_LOG")
                .unwrap_or_else(|_| "example_tracing_aka_logging=debug,tower_http=debug".into()),
        ))
        .with(tracing_subscriber::fmt::layer())
        .init();

    let mut seq_db = SeqIndexDB::new();
   
    if opt.frg_file {
        let _ = seq_db.load_from_frg_index(opt.data_path_prefix); 
    } else {
        #[cfg(feature = "with_agc")]
        let _ = seq_db.load_from_agc_index(opt.data_path_prefix);

        #[cfg(not(feature = "with_agc"))]
        panic!("This command is compiled with only frg file support, please specify `--frg-file");
    }


    let seq_db = Arc::new(seq_db);
    // build our application with a route
    let app = Router::new()
        .route(
            "/api/get_number_of_ctgs",
            get({
                let seq_db = seq_db.clone();
                move || get_number_of_ctgs(seq_db)
            }),
        )
        .route(
            "/api/post_query_for_json_data",
            post({
                let seq_db = seq_db.clone();
                move |params| post_query_for_json_data(params, seq_db)
            }),
        )
        .route(
            "/api/get_html_by_query",
            get({
                let seq_db = seq_db.clone();
                move |params| get_html_by_query(params, seq_db)
            }),
        )
        .fallback(get(|req| async move {
            match ServeDir::new(&opt.static_dir).oneshot(req).await {
                Ok(res) => {
                    let status = res.status();
                    match status {
                        StatusCode::NOT_FOUND => {
                            let index_path = PathBuf::from(&opt.static_dir).join("index.html");
                            let index_content = match fs::read_to_string(index_path).await {
                                Err(_) => {
                                    return Response::builder()
                                        .status(StatusCode::NOT_FOUND)
                                        .body(boxed(Body::from("index file not found")))
                                        .unwrap()
                                }
                                Ok(index_content) => index_content,
                            };

                            Response::builder()
                                .status(StatusCode::OK)
                                .body(boxed(Body::from(index_content)))
                                .unwrap()
                        }
                        _ => res.map(boxed),
                    }
                }
                Err(err) => Response::builder()
                    .status(StatusCode::INTERNAL_SERVER_ERROR)
                    .body(boxed(Body::from(format!("error: {err}"))))
                    .expect("error response"),
            }
        }))
        .layer(
            CorsLayer::new()
                .allow_origin(Any)
                //.allow_origin("http://127.0.0.1:8080".parse::<HeaderValue>().unwrap())
                .allow_methods(Any)
                .allow_headers(Any),
        )
        .layer(ServiceBuilder::new().layer(TraceLayer::new_for_http()));

    // run it
    let addr = SocketAddr::from((
        IpAddr::from_str(opt.addr.as_str()).unwrap_or(IpAddr::V6(Ipv6Addr::LOCALHOST)),
        opt.port,
    ));
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

async fn get_number_of_ctgs(seq_db: Arc<SeqIndexDB>) -> Json<usize> {
    let n_ctg = seq_db.seq_index.as_ref().unwrap().len();
    Json(n_ctg)
}

async fn post_query_for_json_data(
    Json(seq_query_spec): Json<Option<SequenceQuerySpec>>,
    seq_db: Arc<SeqIndexDB>,
) -> Json<Option<TargetMatchPrincipalBundles>> {
    if seq_query_spec.is_none() {
        return Json(None);
    };

    let seq_query_spec = seq_query_spec.unwrap();
    println!("{:?}", seq_query_spec);
    Json(get_target_and_principal_bundle_decomposition(
        &seq_query_spec,
        seq_db,
    ))
}

async fn get_html_by_query(
    Query(seq_query_spec): Query<SequenceQuerySpec>,
    seq_db: Arc<SeqIndexDB>,
) -> Html<String> {
    //if seq_query_spec.is_none() {
    //    return Html("<html><body>No Query Yet</body></html>".into());
    //};

    //let seq_query_spec = seq_query_spec.unwrap();
    println!("{:?}", seq_query_spec);

    let data = get_target_and_principal_bundle_decomposition(&seq_query_spec, seq_db);
    let output = pb_data_to_html_string(&data.unwrap());

    Html(output)
}
