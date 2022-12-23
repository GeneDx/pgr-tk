const VERSION_STRING: &'static str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader};
use std::{fs::File, path};
use svg::node::{element, Node};
use svg::Document;

#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-svg")]
#[clap(author, version)]
#[clap(about = "generate SVG from a principal bundle bed file", long_about = None)]
struct CmdOptions {
    bed_file_path: String,
    output_prefix: String,
    #[clap(long, default_value_t = 100000)]
    track_range: usize,
    #[clap(long, default_value_t = 1600)]
    track_length: usize,
    #[clap(long)]
    left_padding: Option<usize>,
    #[clap(long, default_value_t = 0.5)]
    stroke_width: f32,
}

static cmap: [&str; 97] = [
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

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let bed_file_path = path::Path::new(&args.bed_file_path);
    let bed_file = BufReader::new(File::open(bed_file_path)?);
    let mut ctg_data = FxHashMap::<String, Vec<_>>::default();
    bed_file.lines().into_iter().for_each(|line| {
        let line = line.unwrap();
        let bed_fields = line.split("\t").collect::<Vec<&str>>();
        let ctg: String = bed_fields[0].to_string();
        let bgn: u32 = bed_fields[1].parse().unwrap();
        let end: u32 = bed_fields[2].parse().unwrap();
        let pbundle_fields = bed_fields[3].split(":").collect::<Vec<&str>>();
        let bundle_id: u32 = pbundle_fields[0].parse().unwrap();
        let bundle_dir: u32 = pbundle_fields[2].parse().unwrap();
        let e = ctg_data.entry(ctg).or_insert(vec![]);
        e.push((bgn, end, bundle_id, bundle_dir));
    });

    let left_padding = if args.left_padding.is_some() {
        args.left_padding.unwrap()
    } else {
        10000
    };

    let scaling_factor = args.track_length as f32 / (args.track_range + 2 * left_padding) as f32;
    let left_padding = left_padding as f32 * scaling_factor as f32;
    let stroke_width = 0.5;

    ctg_data
        .into_iter()
        .enumerate()
        .for_each(|(id, (ctg, bundle_segment))| {
            let mut document = Document::new()
                .set("viewBox", (0, -16, args.track_length, 24))
                .set("width", args.track_length)
                .set("height", 40)
                .set("preserveAspectRatio", "none");

            bundle_segment
                .into_iter()
                .for_each(|(bgn, end, bundle_id, direction)| {
                    let mut bgn = bgn as f32 * scaling_factor + left_padding;
                    let mut end = end as f32 * scaling_factor + left_padding;
                    if direction == 1 {
                        (bgn, end) = (end, bgn);
                    }

                    let bundle_color = cmap[(bundle_id % 97) as usize];
                    let stroke_color = cmap[((bundle_id * 47) % 43) as usize];
                    let arror_end = end as f32;
                    let end = if direction == 0 {
                        if end as f32 - 5.0 < bgn {
                            bgn
                        } else {
                            end as f32 - 5.0
                        }
                    } else {
                        if end as f32 + 5.0 > bgn {
                            bgn
                        } else {
                            end as f32 + 5.0
                        }
                    };

                    let path_str = format!(
					"M {bgn} -3 L {bgn} 3 L {end} 3 L {end} 4 L {arror_end} 0 L {end} -4 L {end} -3 Z"
				);

                    let path = element::Path::new()
                        .set("fill", bundle_color)
                        .set("stroke", stroke_color)
                        .set("stroke-width", stroke_width)
                        .set("d", path_str);
                    document.append(path);
                });

            let ext = format!("{:03}.svg", id);
            let out_path = path::Path::new(&args.output_prefix).with_extension(ext);
            svg::save(out_path, &document).unwrap();
        });
    Ok(())
}
