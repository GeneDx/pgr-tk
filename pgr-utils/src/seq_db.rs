#![allow(dead_code)]

use crate::fasta_io::FastaReader;
use crate::seqmap::Shmmrs;
use crate::shmmrutils::{sequence_to_shmmrs, ShmmrSpec, MM128};
use flate2::bufread::MultiGzDecoder;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};

pub struct SeqDB {
    pub seqs: Vec<Vec<u8>>,
    pub shmmrs: Shmmrs,
    pub filepath: String,
    pub lengths: FxHashMap<String, usize>,
    pub id2name: FxHashMap<u32, String>,
    pub name2id: FxHashMap<String, u32>,
    pub sequences_loaded: bool,
    pub shmmrs_built: bool,
}

impl SeqDB {
    pub fn new(filepath: String) -> Self {
        let seqs = Vec::<Vec<u8>>::default();
        let shimmers_db = Vec::<Vec<MM128>>::new();
        let seqlen = FxHashMap::<String, usize>::default();
        let name2id = FxHashMap::<String, u32>::default();
        let id2name = FxHashMap::<u32, String>::default();

        SeqDB {
            filepath: filepath,
            seqs: seqs,
            shmmrs: shimmers_db,
            lengths: seqlen,
            id2name,
            name2id,
            sequences_loaded: false,
            shmmrs_built: false,
        }
    }

    pub fn load_sequences(&mut self) -> Result<(), std::io::Error> {
        let file = File::open(&self.filepath)?;

        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!(
                    "input file: {} detected as gz-compressed file",
                    self.filepath
                );
                is_gzfile = true;
            }
        }

        reader.seek(SeekFrom::Start(0))?;
        if is_gzfile {
            let fastx_buf = BufReader::new(MultiGzDecoder::new(&mut reader));
            let mut fastx_reader = FastaReader::new(fastx_buf, &self.filepath, 1 << 14, true)?;
            let mut sid = 0;
            while let Some(rec) = fastx_reader.next_rec() {
                let rec = rec.unwrap();
                let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                self.id2name.insert(sid, seqname.clone());
                self.name2id.insert(seqname.clone(), sid);
                self.lengths.insert(seqname.clone(), rec.seq.len());
                self.seqs.push(rec.seq);
                sid += 1;
            }
        } else {
            let mut fastx_reader = FastaReader::new(reader, &self.filepath, 1 << 14, true).unwrap();
            let mut sid = 0;
            // unfortunatly we, need to repeat the code here as the type fastx_reader is different from the above
            while let Some(rec) = fastx_reader.next_rec() {
                let rec = rec.unwrap();
                let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                self.id2name.insert(sid, seqname.clone());
                self.name2id.insert(seqname.clone(), sid);
                self.lengths.insert(seqname.clone(), rec.seq.len());
                self.seqs.push(rec.seq);
                sid += 1;
            }
        }
        self.sequences_loaded = true;
        Ok(())
    }

    pub fn build_shmmrs(&mut self, w: u32, k: u32, r: u32) -> () {
        let e_seqs = self
            .seqs
            .iter()
            .enumerate()
            .collect::<Vec<(usize, &Vec<u8>)>>();

        let mut out = e_seqs
            .par_chunks(1)
            .into_par_iter()
            .flat_map(|xv| {
                let mut out = Vec::<(usize, Vec<MM128>)>::new();
                for x in xv.iter() {
                    out.push((
                        x.0,
                        sequence_to_shmmrs(
                            x.0 as u32,
                            x.1,
                            ShmmrSpec {
                                w,
                                k,
                                r,
                                min_span: 32,
                                sketch: false,
                            },
                        ),
                    ));
                }
                out
            })
            .collect::<Vec<(usize, Vec<MM128>)>>();

        out.par_sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        self.shmmrs.clear();
        for (_sid, shmmrs) in out {
            let mut shmmrs_f = Vec::<MM128>::new();

            if shmmrs.len() > 0 {
                let mut p_shmmr = shmmrs[0];
                shmmrs_f.push(p_shmmr);
                shmmrs[1..shmmrs.len()].iter().for_each(|&shmmr| {
                    if shmmr.x >> 8 != p_shmmr.x >> 8 {
                        // filter out neighboring identical mmers, note this operation is not strand skew symmtrical
                        p_shmmr = shmmr.clone();
                        shmmrs_f.push(p_shmmr);
                    }
                });
            }
            self.shmmrs.push(shmmrs_f);
        }

        self.shmmrs_built = true;
        ()
    }

    pub fn get_subseq_by_name(&self, seqname: String, s: usize, e: usize) -> Option<String> {
        if let Some(sid) = self.name2id.get(&seqname) {
            if e > *self.lengths.get(&seqname).unwrap() {
                return None;
            } else {
                let str = String::from_utf8_lossy(&self.seqs.get(*sid as usize).unwrap()[s..e])
                    .to_string();
                Some(str)
            }
        } else {
            None
        }
    }

    pub fn get_subseq_by_id(&self, sid: u32, b: usize, e: usize) -> Option<String> {
        if sid < self.seqs.len() as u32 {
            let seqname = self.id2name.get(&sid).unwrap();
            if e > *self.lengths.get(seqname).unwrap() {
                None
            } else {
                let str = String::from_utf8_lossy(&self.seqs[sid as usize][b..e]).to_string();
                Some(str)
            }
        } else {
            None
        }
    }
}
