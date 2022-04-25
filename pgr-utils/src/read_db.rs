#![allow(dead_code)]

// src/lib.rs
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

// use super::seqmap::{self, MapIntervalRecord};
use super::shmmrutils::{sequence_to_shmmrs, MM128};
use super::seqmap::Shmmrs;
use memmap::{Mmap, MmapOptions};
use rayon::prelude::*;

trait HasShmmer {
    fn get_shmmr(&self) -> &Shmmrs;
}
pub struct ReadDB {
    pub seqs: Mmap,
    pub shmmrs: Shmmrs,
    pub seqdb_filepath: String,
    pub seqidx_filepath: String,
    pub seqstart: FxHashMap<String, usize>,
    pub seqlen: FxHashMap<String, usize>,
    pub id2seqname: FxHashMap<u32, String>,
    pub seqname2id: FxHashMap<String, u32>,
    pub index_loaded: bool,
    pub shmmrs_built: bool,
}

impl HasShmmer for ReadDB {
    fn get_shmmr(&self) -> &Shmmrs {
        self.shmmrs.as_ref()
    }
}

impl ReadDB {
    pub fn new(seqdb_filepath: String, seqidx_filepath: String) -> Self {
        let file = File::open(seqdb_filepath.clone()).unwrap();
        let seqdb_mmap = unsafe { MmapOptions::new().map(&file).unwrap() };
        let shmmrs_db = Vec::<Vec<MM128>>::new();
        let mut seqlen = FxHashMap::<String, usize>::default();
        let mut seqstart = FxHashMap::<String, usize>::default();
        let mut seqname2id = FxHashMap::<String, u32>::default();
        let mut id2seqname = FxHashMap::<u32, String>::default();

        let indexfile = File::open(seqidx_filepath.clone()).unwrap();
        let reader = BufReader::new(indexfile);
        for line in reader.lines() {
            if let Ok(rec) = line {
                let v: Vec<&str> = rec.split_whitespace().collect();
                let rid: u32 = v[0].parse().unwrap();
                let readname = String::from(v[1]);
                let len: usize = v[2].parse().unwrap();
                let start: usize = v[3].parse().unwrap();
                seqlen.insert(readname.clone(), len);
                seqstart.insert(readname.clone(), start);
                id2seqname.insert(rid, readname.clone());
                seqname2id.insert(readname, rid);
            }
        }

        ReadDB {
            seqs: seqdb_mmap,
            seqdb_filepath: seqdb_filepath,
            seqidx_filepath: seqidx_filepath,
            shmmrs: shmmrs_db,
            seqstart: seqstart,
            seqlen: seqlen,
            id2seqname: id2seqname,
            seqname2id: seqname2id,
            index_loaded: true,
            shmmrs_built: false,
        }
    }

    pub fn build_shmmrs_parallel(&mut self, w: u32, k: u32, r: u32) -> () {
        let mut seqidx = Vec::<(u32, usize, usize)>::new();
        let rids = self.id2seqname.iter().collect::<Vec<(&u32, &String)>>();
        for (rid, rname) in rids {
            let s = *self.seqstart.get(rname).unwrap();
            let l = *self.seqlen.get(rname).unwrap();
            seqidx.push((*rid, s, s + l));
        }
        let seqdb_mmap = self.seqs.as_ref();
        let seqs = std::sync::Arc::new(seqdb_mmap);
        //let base_map = &[b'\0', b'A', b'C', b'\0', b'G', b'\0', b'\0', b'\0', b'T'];
        let base_map = &[b'A', b'C', b'G', b'T'];

        let mut out = seqidx
            .par_iter()
            .map(move |&x| {
                let seq = seqs[x.1..x.2]
                    .iter()
                    .map(|c| base_map[(*c & 0b0011) as usize])
                    .collect::<Vec<u8>>();
                (
                    x.0 as usize,
                    sequence_to_shmmrs(x.0 as u32, &seq, w, k, r),
                )
            })
            .collect::<Vec<(usize, Vec<MM128>)>>();

        out.par_sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        self.shmmrs.clear();
        for (_, shmmrs) in out {
            self.shmmrs.push(shmmrs);
        }
        self.shmmrs_built = true;
        ()
    }

    pub fn get_seqname_by_id(&self, sid: u32) -> String {
        self.id2seqname.get(&sid).unwrap().clone()
    }

    pub fn get_id_by_seqname(&self, seqname: String) -> u32 {
        self.seqname2id.get(&seqname).unwrap().clone()
        
    }

    pub fn get_seqlen_by_id(&self, sid: u32) -> usize {
        let seqname = self.id2seqname.get(&sid).unwrap();
        self.seqlen.get(seqname).unwrap().clone()
    }

    pub fn get_seqlen_by_seqname(&self, seqname: String) -> usize {
        self.seqlen.get(&seqname).unwrap().clone()
    }

    pub fn get_seq_by_name(&self, seqname: String) -> String {
        //let base_map = &[b'\0', b'A', b'C', b'\0', b'G', b'\0', b'\0', b'\0', b'T'];
        let base_map = &[b'A', b'C', b'G', b'T'];
        let b = *self.seqstart.get(&seqname).unwrap();
        let l = *self.seqlen.get(&seqname).unwrap();
        let seq = &self.seqs[b..(b + l)]
            .iter()
            .map(|c| base_map[(c & 0b0011) as usize])
            .collect::<Vec<u8>>();
        String::from_utf8_lossy(&seq).to_string()
    }

    pub fn get_base_by_id(&self, sid: u32, pos: usize, strand: u32) -> char {
        let sname = self.id2seqname.get(&sid).unwrap().clone();
        let base_map = &[b'A', b'C', b'G', b'T'];
        let b = *self.seqstart.get(&sname).unwrap();
        let l = *self.seqlen.get(&sname).unwrap();
        let base: u8;
        if strand == 0 {
            base = base_map[ (&self.seqs[b..(b + l)][pos] & 0b0011) as usize ];
        } else {
            base = base_map[ ( ((&self.seqs[b..(b + l)][pos]) >> 4) & 0b0011) as usize ];
        }
        base as char
    }

    pub fn get_seq_by_id(&self, sid: u32) -> String {
        let sname = self.id2seqname.get(&sid).unwrap();
        self.get_seq_by_name(sname.clone())
    }

    pub fn get_all_seqnames(&self) -> Vec<String> {
        let mut names = Vec::<String>::new();
        for n in self.seqname2id.iter() {
            names.push(n.0.clone());
        }
        names
    }

    pub fn get_all_ids(&self) -> Vec<u32> {
        let mut ids = Vec::<u32>::new();
        for n in self.seqname2id.iter() {
            ids.push(*n.1);
        }
        ids
    }

    pub fn get_shmmer_locs_by_id(&self, sid: u32) -> Vec<u32> {
        let mut locs = Vec::<u32>::new();
        for shmer in &self.shmmrs[sid as usize] {
            locs.push(((shmer.y & 0xFFFFFFFF) >> 1) as u32);
        }
        locs
    }
}