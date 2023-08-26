const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Deserialize;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, Read};
use std::path::{self, Path};
use svg::node::{self, element, Node};
use svg::Document;

#[allow(dead_code)] // need the standard names for deserization if they are not use 
#[derive(Deserialize, Clone, Debug)]
struct CtgMapRec {
    t_name: String,
    ts: u32,
    te: u32,
    q_name: String,
    qs: u32,
    qe: u32,
    ctg_len: u32,
    orientation: u32,
    ctg_orientation: u32,
    t_dup: bool,
    t_ovlp: bool,
    q_dup: bool,
    q_ovlp: bool,
}

#[derive(Deserialize)]
struct CtgMapSet {
    records: Vec<CtgMapRec>,
    target_length: Vec<(u32, String, u32)>,
    query_length: Vec<(u32, String, u32)>,
}

/// generate align block plot from ctgmap.json file
#[derive(Parser, Debug)]
#[clap(name = "pgr-generate-chr-aln-plot")]
#[clap(author, version)]
#[clap(about, long_about = None)]

struct CmdOptions {
    /// path to a ctgmap.json file
    ctgmap_json_path: String,
    /// the prefix of the output files
    output_path: String,
}

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

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    let mut ctgmap_json_file = BufReader::new(
        File::open(Path::new(&args.ctgmap_json_path)).expect("can't open the input file"),
    );
    let mut buffer = Vec::new();
    ctgmap_json_file.read_to_end(&mut buffer)?;
    let mut ctgmap_set: CtgMapSet = serde_json::from_str(&String::from_utf8_lossy(&buffer[..]))
        .expect("can't parse the ctgmap.json file");

    ctgmap_set.query_length.sort();
    ctgmap_set.target_length.sort();
    let mut ctg_target_hit_len = FxHashMap::<String, FxHashMap<String, u32>>::default();

    let query_length = ctgmap_set
        .query_length
        .iter()
        .map(|v| (v.1.clone(), v.2))
        .collect::<FxHashMap<_, _>>();

    ctgmap_set.records.iter().for_each(|r| {
        if r.q_dup {
            return;
        };
        let e = ctg_target_hit_len.entry(r.q_name.clone()).or_default();
        let e2 = e.entry(r.t_name.clone()).or_default();
        *e2 += (r.qe as i32 - r.qs as i32).unsigned_abs();
    });

    let mut ctg2tgt = FxHashMap::<String, String>::default();
    let mut tgt2ctg = FxHashMap::<String, String>::default();

    ctg_target_hit_len.into_iter().for_each(|(ctg, tgt_len)| {
        let mut tgt_len = tgt_len.into_iter().collect::<Vec<_>>();
        if !tgt_len.is_empty() {
            tgt_len.sort_by(|a, b| b.1.cmp(&a.1));
            let tgt = tgt_len[0].0.clone();
            //println!("DBG: {} {:?}", ctg, tgt_len);
            ctg2tgt.insert(ctg.clone(), tgt.clone());
            tgt2ctg.insert(tgt.clone(), ctg.clone());
        };
    });

    let mut tgt_to_records = FxHashMap::<String, Vec<CtgMapRec>>::default();
    ctgmap_set.records.iter().for_each(|r| {
        if r.q_dup {
            return;
        };
        if *ctg2tgt.get(&r.q_name).unwrap() != r.t_name {
            return;
        }
        let e = tgt_to_records.entry(r.t_name.clone()).or_default();
        e.push((*r).clone());
    });

    let target_padding = 1.5e6;
    let mut offset = 0_f64;
    let target_aln_blocks = ctgmap_set
        .target_length
        .iter()
        .flat_map(|(id, t_name, t_len)| {
            let mut q_len_sum = 0.0;
            let mut q_set = FxHashSet::<String>::default();
            tgt_to_records
                .get(t_name)
                .unwrap_or(&vec![])
                .iter()
                .for_each(|record| {
                    let q_len = query_length.get(&record.q_name).unwrap();
                    if !q_set.contains(&record.q_name) {
                        q_set.insert(record.q_name.clone());
                        q_len_sum += *q_len as f64;
                    };
                });

            if let Some(records) = tgt_to_records.get(t_name) {
                let out = (*id, t_name.clone(), *t_len, offset, records);
                offset += (*t_len as f64).max(q_len_sum) + target_padding;
                Some(out)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    // start to construct the SVG element
    let mut document = Document::new()
        .set("viewBox", (-100, -50, 1300, 250 * 14))
        .set("width", 1400)
        .set("height", 250 * 14)
        .set("preserveAspectRatio", "none")
        .set("id", "WholeGenomeViwer");

    let scaling_factor = 1100.0 / offset;

    target_aln_blocks
        .iter()
        .for_each(|target_aln_block_records| {
            let t_offset = target_aln_block_records.3;

            let b = t_offset * scaling_factor;
            let e = (t_offset + target_aln_block_records.2 as f64) * scaling_factor;
            let w = 4.0 + ((target_aln_block_records.0 + 1) % 2) as f64 * 1.5;
            let path_str = format!("M {b} 6 L {e} 6");
            let path = element::Path::new()
                .set("stroke", "#000")
                .set("stroke-width", format!("{w}"))
                .set("opacity", "0.7")
                .set("d", path_str);
            document.append(path);

            let text = element::Text::new()
                .set("x", b)
                .set("y", 0)
                .set("font-size", "6px")
                .set("font-family", "monospace")
                .add(node::Text::new(target_aln_block_records.1.clone()));
            document.append(text);

            let mut best_query_block = FxHashMap::<String, CtgMapRec>::default();
            target_aln_block_records.4.iter().for_each(|record| {
                let e = best_query_block
                    .entry(record.q_name.clone())
                    .or_insert(record.clone());
                if (e.qs as i32 - e.qe as i32).abs() < (record.qs as i32 - record.qe as i32).abs() {
                    *e = record.clone();
                }
            });

            let mut best_query_block = best_query_block.values().collect::<Vec<_>>();
            best_query_block.sort_by_key(|&v| v.ts);
            let mut q_offset = 0.0;
            let mut q_offset_map = FxHashMap::<String, f64>::default();
            best_query_block.into_iter().for_each(|record| {
                let q_len = query_length.get(&record.q_name).unwrap();
                if !q_offset_map.contains_key(&record.q_name) {
                    q_offset_map.insert(record.q_name.clone(), q_offset);

                    let b = (t_offset + q_offset) * scaling_factor;
                    let e = (t_offset + q_offset + *q_len as f64) * scaling_factor;
                    let y = 95.0;
                    let path_str = format!("M {b} {y} L {e} {y}");
                    let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
                    let path = element::Path::new()
                        .set("stroke", color)
                        .set("stroke-width", "5")
                        .set("opacity", "0.7")
                        .set("d", path_str);
                    document.append(path);

                    q_offset += *q_len as f64;
                };
            });

            target_aln_block_records.4.iter().for_each(|record| {
                if record.t_dup {
                    return;
                };

                let q_len = query_length.get(&record.q_name).unwrap();

                let ts = record.ts as f64 + t_offset;
                let te = record.te as f64 + t_offset;

                let (qs, qe) = if record.ctg_orientation == 1 {
                    (q_len - record.qe, q_len - record.qs)
                } else {
                    (record.qs, record.qe)
                };

                // let qs = record.qs;
                // let qe = record.qe;
                let (qs, qe) = if record.orientation != record.ctg_orientation {
                    (qe, qs)
                } else {
                    (qs, qe)
                };
                let offset = q_offset_map.get(&record.q_name).unwrap();
                let qs = qs as f64 + t_offset + offset;
                let qe = qe as f64 + t_offset + offset;
                let ts = ts * scaling_factor;
                let te = te * scaling_factor;
                let qs = qs * scaling_factor;
                let qe = qe * scaling_factor;
                // println!("{:?}", record);
                // println!("{} {} {} {}", ts, te, qs, qe);

                let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];

                let path_str = format!("M {ts} 10 L {te} 10 L {qe} 90 L {qs} 90");
                let path = element::Path::new()
                    .set("fill", color)
                    .set("stroke", "#000")
                    .set("stroke-width", "0.5")
                    .set("opacity", "0.7")
                    .set("d", path_str);
                document.append(path);
            });
        });

    // per chromo plot
    let mut y_offset = 200.0;
    let scaling_factor = 1100.0 / offset * 12.0;

    target_aln_blocks
        .iter()
        .for_each(|target_aln_block_records| {
            let t_offset = 0.0;

            let b = t_offset * scaling_factor;
            let e = (t_offset + target_aln_block_records.2 as f64) * scaling_factor;
            let w = 4.0 + ((target_aln_block_records.0 + 1) % 2) as f64 * 1.5;
            let y = y_offset + 6.0;
            let path_str = format!("M {b} {y} L {e} {y}");
            let path = element::Path::new()
                .set("stroke", "#000")
                .set("stroke-width", format!("{w}"))
                .set("opacity", "0.7")
                .set("d", path_str);
            document.append(path);

            let text = element::Text::new()
                .set("x", 0.0)
                .set("y", y - 10.0)
                .set("font-size", "20px")
                .set("font-family", "monospace")
                .add(node::Text::new(target_aln_block_records.1.clone()));
            document.append(text);

            let mut best_query_block = FxHashMap::<String, CtgMapRec>::default();
            target_aln_block_records.4.iter().for_each(|record| {
                let e = best_query_block
                    .entry(record.q_name.clone())
                    .or_insert(record.clone());
                if (e.qs as i32 - e.qe as i32).abs() < (record.qs as i32 - record.qe as i32).abs() {
                    *e = record.clone();
                }
            });

            let mut best_query_block = best_query_block.values().collect::<Vec<_>>();
            best_query_block.sort_by_key(|&v| v.ts);
            let mut q_offset = 0.0;
            let mut q_offset_map = FxHashMap::<String, f64>::default();
            best_query_block.into_iter().for_each(|record| {
                let q_len = query_length.get(&record.q_name).unwrap();
                if !q_offset_map.contains_key(&record.q_name) {
                    q_offset_map.insert(record.q_name.clone(), q_offset);

                    let b = (t_offset + q_offset) * scaling_factor;
                    let e = (t_offset + q_offset + *q_len as f64) * scaling_factor;
                    let y = 95.0 + y_offset;
                    let path_str = format!("M {b} {y} L {e} {y}");
                    let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
                    let mut path = element::Path::new()
                        .set("stroke", color)
                        .set("stroke-width", "5")
                        .set("opacity", "0.7")
                        .set("d", path_str);
                    path.append(element::Title::new().add(node::Text::new(record.q_name.clone())));
                    document.append(path);

                    q_offset += *q_len as f64;
                };
            });

            target_aln_block_records.4.iter().for_each(|record| {
                if record.t_dup {
                    return;
                };

                let q_len = query_length.get(&record.q_name).unwrap();

                let ts = record.ts as f64 + t_offset;
                let te = record.te as f64 + t_offset;

                let (qs, qe) = if record.ctg_orientation == 1 {
                    (q_len - record.qe, q_len - record.qs)
                } else {
                    (record.qs, record.qe)
                };

                // let qs = record.qs;
                // let qe = record.qe;
                let (qs, qe) = if record.orientation != record.ctg_orientation {
                    (qe, qs)
                } else {
                    (qs, qe)
                };
                let offset = q_offset_map.get(&record.q_name).unwrap();
                let qs = qs as f64 + t_offset + offset;
                let qe = qe as f64 + t_offset + offset;
                let ts = ts * scaling_factor;
                let te = te * scaling_factor;
                let qs = qs * scaling_factor;
                let qe = qe * scaling_factor;
                // println!("{:?}", record);
                // println!("{} {} {} {}", ts, te, qs, qe);

                let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
                let y = 10.0 + y_offset;
                let y2 = 90.0 + y_offset;
                let path_str = format!("M {ts} {y} L {te} {y} L {qe} {y2} L {qs} {y2}");
                let mut path = element::Path::new()
                    .set("fill", color)
                    .set("stroke", "#000")
                    .set("stroke-width", "0.5")
                    .set("opacity", "0.7")
                    .set("d", path_str);
                let orientation = if record.orientation == 0 { '+' } else { '-' };
                path.append(element::Title::new().add(node::Text::new(format!(
                    "{}:{}-{} @ {}:{}-{} {}",
                    record.t_name, record.ts, record.te,
                    record.q_name, record.qs, record.qe, 
                    orientation
                ))));

                document.append(path);
            });
            y_offset += 125.0;
        });

    let out_path = path::Path::new(&args.output_path);
    svg::save(out_path, &document).unwrap();

    Ok(())
}
