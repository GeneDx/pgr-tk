// src/lib.rs
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::PyString;
use pyo3::wrap_pyfunction;

use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
use peregrine_utils::fasta_io::FastaReader;
use peregrine_utils::seqmap::{self, MapIntervalRecord};
use peregrine_utils::multi_seqmap;

use peregrine_utils::shmmrutils::{sequence_to_shmmrs, MM128};
use pyo3::Python;
use rustc_hash::FxHashMap;
use std::fs::File;
use flate2::bufread::MultiGzDecoder;
use std::io::{BufReader, BufRead, Read, Seek, SeekFrom};
type Shmmrs = peregrine_utils::seqmap::Shmmrs;
use memmap::{Mmap, MmapOptions};
use rayon::prelude::*;

use peregrine_utils::seqs2variants;

trait HasShmmer {
    fn get_shmmr(&self) -> &Shmmrs;
}

#[pyclass]
struct SeqDB {
    seqs: Vec<Vec<u8>>,
    shmmrs: Shmmrs,
    #[pyo3(get)]
    filepath: String,
    lengths: FxHashMap<String, usize>,
    id2name: FxHashMap<u32, String>,
    name2id: FxHashMap<String, u32>,
    #[pyo3(get)]
    sequences_loaded: bool,
    #[pyo3(get)]
    shmmrs_built: bool,
}

impl HasShmmer for SeqDB {
    fn get_shmmr(&self) -> &Shmmrs {
        self.shmmrs.as_ref()
    }
}
#[pyclass]
struct MapIntervals {
    intervals: seqmap::MapIntervals,
}

#[pymethods]
impl MapIntervals {
    #[new]
    fn new() -> Self {
        MapIntervals {
            intervals: seqmap::MapIntervals::default(),
        }
    }
}

#[pyclass]
#[derive(Clone)]
struct AlnSegment {
    #[pyo3(get, set)]
    t: u8,
    #[pyo3(get, set)]
    ref_loc: (u32, u32, u32),
    #[pyo3(get, set)]
    tgt_loc: (u32, u32, u32),
}

#[pyclass]
#[derive(Clone)]
pub struct AlnMap {
    #[pyo3(get, set)]
    pmap: Vec<(u32, u32)>,
    #[pyo3(get, set)]
    ref_a_seq: Vec<u8>,
    #[pyo3(get, set)]
    tgt_a_seq: Vec<u8>,
    #[pyo3(get, set)]
    aln_seq: Vec<u8>,
}

#[pymethods]
impl SeqDB {
    #[new]
    fn new(filepath: String) -> Self {
        let seqs = Vec::<Vec<u8>>::default();
        let shmmrs_db = Vec::<Vec<MM128>>::new();
        let seqlen = FxHashMap::<String, usize>::default();
        let name2id = FxHashMap::<String, u32>::default();
        let id2name = FxHashMap::<u32, String>::default();

        SeqDB {
            filepath: filepath,
            seqs: seqs,
            shmmrs: shmmrs_db,
            lengths: seqlen,
            id2name,
            name2id,
            sequences_loaded: false,
            shmmrs_built: false,
        }
    }

    fn load_sequences(&mut self) -> PyResult<()> {
        let file = File::open(&self.filepath)?;
        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                is_gzfile = true;
            }
        }

        reader.seek(SeekFrom::Start(0))?;
        if is_gzfile {
            let fastx_buf = BufReader::new(MultiGzDecoder::new(&mut reader));
            let mut fastx_reader = FastaReader::new(fastx_buf, &self.filepath)?;
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
            let mut fastx_reader = FastaReader::new(reader, &self.filepath).unwrap();
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

    fn build_shmmrs(&mut self, w: u32, k: u32, r: u32) -> PyResult<()> {
        let mut len = 0_usize;
        self.shmmrs.clear();
        for (sid, seq) in self.seqs.iter().enumerate() {
            //let sid = self.name2id.get(sname).unwrap();
            let shmmers = sequence_to_shmmrs(sid as u32, seq, w, k, r);
            self.shmmrs.push(shmmers);

            // allow to catch Python interuption whne processed enough sequences
            len += seq.len();
            if len > 10_000_000 {
                Python::with_gil(|py| -> PyResult<()> {
                    py.check_signals()?;
                    Ok(())
                })?;
                len = 0;
            }
        }
        self.shmmrs_built = true;
        Ok(())
    }

    fn build_shmmrs_parallel(&mut self, w: u32, k: u32, r: u32) -> () {
        let e_seqs = self
            .seqs
            .iter()
            .enumerate()
            .collect::<Vec<(usize, &Vec<u8>)>>();

        let mut out = e_seqs
            .par_iter()
            .map(|&x| (x.0, sequence_to_shmmrs(x.0 as u32, x.1, w, k, r)))
            .collect::<Vec<(usize, Vec<MM128>)>>();

        out.par_sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        self.shmmrs.clear();
        for (_, shmmrs) in out {
            self.shmmrs.push(shmmrs);
        }
        self.shmmrs_built = true;
        ()
    }

    fn get_name_by_id(&self, sid: u32) -> PyResult<String> {
        let name = self.id2name.get(&sid).unwrap();
        Ok(name.clone())
    }

    fn get_id_by_name(&self, name: String) -> PyResult<u32> {
        let sid = self.name2id.get(&name).unwrap();
        Ok(*sid)
    }

    fn get_subseq_by_name(&self, name: String, s: usize, e: usize) -> PyResult<String> {
        let sid = self.name2id.get(&name).unwrap();
        let str = String::from_utf8_lossy(&self.seqs.get(*sid as usize).unwrap()[s..e]).to_string();
        Ok(str)
    }

    fn get_subseq_by_id(&self, sid: u32, s: usize, e: usize) -> PyResult<String> {
        let str = String::from_utf8_lossy(&self.seqs[sid as usize][s..e]).to_string();
        Ok(str)
    }

    fn get_len_by_id(&self, sid: u32) -> PyResult<usize> {
        let name = self.id2name.get(&sid).unwrap();
        let len = self.lengths.get(name).unwrap();
        Ok(*len)
    }

    fn get_len_by_name(&self, name: String) -> PyResult<usize> {
        let len = self.lengths.get(&name).unwrap();
        Ok(*len)
    }

    fn get_seq_by_name(&self, name: String) -> PyResult<String> {
        let sid = self.name2id.get(&name).unwrap();
        let str = String::from_utf8_lossy(&self.seqs.get(*sid as usize).unwrap()).to_string();
        Ok(str)
    }

    fn get_seq_by_id(&self, sid: u32) -> PyResult<String> {
        let str = String::from_utf8_lossy(&self.seqs[sid as usize]).to_string();
        Ok(str)
    }

    fn get_all_names(&self) -> PyResult<Vec<String>> {
        let mut names = Vec::<String>::new();
        for n in self.name2id.iter() {
            names.push(n.0.clone());
        }
        Ok(names)
    }

    fn get_all_ids(&self) -> PyResult<Vec<u32>> {
        let mut ids = Vec::<u32>::new();
        for n in self.name2id.iter() {
            ids.push(*n.1);
        }
        Ok(ids)
    }

    fn get_shmmer_locs_by_id(&self, sid: u32) -> PyResult<Vec<u32>> {
        let mut locs = Vec::<u32>::new();
        for shmer in &self.shmmrs[sid as usize] {
            locs.push(((shmer.y & 0xFFFFFFFF) >> 1) as u32);
        }
        Ok(locs)
    }
}

#[pyclass]
struct ReadDB {
    seqs: Mmap,
    shmmrs: Shmmrs,
    #[pyo3(get)]
    seqdb_filepath: String,
    #[pyo3(get)]
    seqidx_filepath: String,
    starts: FxHashMap<String, usize>,
    lengths: FxHashMap<String, usize>,
    id2name: FxHashMap<u32, String>,
    name2id: FxHashMap<String, u32>,
    #[pyo3(get)]
    index_loaded: bool,
    #[pyo3(get)]
    shmmrs_built: bool,
}

impl HasShmmer for ReadDB {
    fn get_shmmr(&self) -> &Shmmrs {
        self.shmmrs.as_ref()
    }
}

#[pymethods]
impl ReadDB {
    #[new]
    fn new(seqdb_filepath: String, seqidx_filepath: String) -> Self {
        let file = File::open(seqdb_filepath.clone()).unwrap();
        let seqdb_mmap = unsafe { MmapOptions::new().map(&file).unwrap() };
        let shmmrs_db = Vec::<Vec<MM128>>::new();
        let mut seqlen = FxHashMap::<String, usize>::default();
        let mut seqstart = FxHashMap::<String, usize>::default();
        let mut name2id = FxHashMap::<String, u32>::default();
        let mut id2name = FxHashMap::<u32, String>::default();

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
                id2name.insert(rid, readname.clone());
                name2id.insert(readname, rid);
            }
        }

        ReadDB {
            seqs: seqdb_mmap,
            seqdb_filepath: seqdb_filepath,
            seqidx_filepath: seqidx_filepath,
            //seqs: seqdb_mmap,
            shmmrs: shmmrs_db,
            starts: seqstart,
            lengths: seqlen,
            id2name,
            name2id,
            index_loaded: true,
            shmmrs_built: false,
        }
    }

    fn build_shmmrs_parallel(&mut self, w: u32, k: u32, r: u32) -> () {
        let mut seqidx = Vec::<(u32, usize, usize)>::new();
        let rids = self.id2name.iter().collect::<Vec<(&u32, &String)>>();
        for (rid, rname) in rids {
            let s = *self.starts.get(rname).unwrap();
            let l = *self.lengths.get(rname).unwrap();
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
                (x.0 as usize, sequence_to_shmmrs(x.0 as u32, &seq, w, k, r))
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

    fn get_name_by_id(&self, sid: u32) -> PyResult<String> {
        let name = self.id2name.get(&sid).unwrap();
        Ok(name.clone())
    }

    fn get_id_by_name(&self, name: String) -> PyResult<u32> {
        let sid = self.name2id.get(&name).unwrap();
        Ok(*sid)
    }

    fn get_seqlen_by_id(&self, sid: u32) -> PyResult<usize> {
        let name = self.id2name.get(&sid).unwrap();
        let len = self.lengths.get(name).unwrap();
        Ok(*len)
    }

    fn get_seqlen_by_name(&self, name: String) -> PyResult<usize> {
        let len = self.lengths.get(&name).unwrap();
        Ok(*len)
    }

    fn get_seq_by_name(&self, name: String) -> PyResult<String> {
        //let base_map = &[b'\0', b'A', b'C', b'\0', b'G', b'\0', b'\0', b'\0', b'T'];
        let base_map = &[b'A', b'C', b'G', b'T'];
        let b = *self.starts.get(&name).unwrap();
        let l = *self.lengths.get(&name).unwrap();
        let seq = &self.seqs[b..(b + l)]
            .iter()
            .map(|c| base_map[(c & 0b0011) as usize])
            .collect::<Vec<u8>>();
        let str = String::from_utf8_lossy(&seq).to_string();
        Ok(str)
    }

    fn get_seq_by_id(&self, sid: u32) -> PyResult<String> {
        let sname = self.id2name.get(&sid).unwrap();
        let str = self.get_seq_by_name(sname.clone())?;
        Ok(str)
    }

    fn get_all_names(&self) -> PyResult<Vec<String>> {
        let mut names = Vec::<String>::new();
        for n in self.name2id.iter() {
            names.push(n.0.clone());
        }
        Ok(names)
    }

    fn get_all_ids(&self) -> PyResult<Vec<u32>> {
        let mut ids = Vec::<u32>::new();
        for n in self.name2id.iter() {
            ids.push(*n.1);
        }
        Ok(ids)
    }

    fn get_shmmer_locs_by_id(&self, sid: u32) -> PyResult<Vec<u32>> {
        let mut locs = Vec::<u32>::new();
        for shmer in &self.shmmrs[sid as usize] {
            locs.push(((shmer.y & 0xFFFFFFFF) >> 1) as u32);
        }
        Ok(locs)
    }
}

#[pyfunction]
fn get_shmmr_dots(
    seq0: &PyString,
    seq1: &PyString,
    w: u32,
    k: u32,
    r: u32,
) -> (Vec<u32>, Vec<u32>) {
    let mut x = Vec::<u32>::new();
    let mut y = Vec::<u32>::new();
    let seq0v = seq0.to_string_lossy().as_bytes().to_vec();
    let seq1v = seq1.to_string_lossy().as_bytes().to_vec();
    let shmmr0 = sequence_to_shmmrs(0, &seq0v, w, k, r);
    let shmmr1 = sequence_to_shmmrs(1, &seq1v, w, k, r);
    let mut basemmer_x = FxHashMap::<u64, Vec<u32>>::default();

    for m in shmmr0 {
        let hash = m.x >> 8;
        let pos = ((m.y & 0xFFFFFFFF) >> 1) as u32;
        basemmer_x.entry(hash).or_insert_with(|| vec![]).push(pos);
    }

    for m in shmmr1 {
        let hash = m.x >> 8;
        let py = ((m.y & 0xFFFFFFFF) >> 1) as u32;
        if basemmer_x.contains_key(&hash) {
            for px in basemmer_x.get(&hash).unwrap() {
                x.push(*px);
                y.push(py);
            }
        }
    }

    (x, y)
}

#[pyfunction]
fn generate_shmmr_map(seqdb0: &SeqDB, seqdb1: &SeqDB, max_hits: usize) -> PyResult<MapIntervals> {
    let all_itvl = seqmap::generate_shmmr_map(&seqdb0.shmmrs, &seqdb1.shmmrs, max_hits);
    Ok(MapIntervals {
        intervals: all_itvl,
    })
}

#[pyfunction]
fn build_shmmer_map_from_query_results(mqr: Vec<MapIntervalRecord>) -> PyResult<MapIntervals> {
    let all_itvl = seqmap::build_shmmer_map_from_query_results(&mqr);

    Ok(MapIntervals {
        intervals: all_itvl,
    })
}

#[pyfunction]
fn generate_shmmr_map_from_read_db(
    seqdb0: &SeqDB,
    readdb: &ReadDB,
    max_hits: usize,
) -> PyResult<MapIntervals> {
    let all_itvl = seqmap::generate_shmmr_map(&seqdb0.shmmrs, &readdb.shmmrs, max_hits);
    Ok(MapIntervals {
        intervals: all_itvl,
    })
}

#[pyfunction]
fn write_shmmr_map(mi: &MapIntervals, filename: &PyString) -> PyResult<()> {
    seqmap::write_shmmr_map(&mi.intervals, &String::from(filename.to_string_lossy()))?;
    Ok(())
}

#[pyfunction]
fn load_shmmr_map(mi: &mut MapIntervals, filename: &PyString) -> PyResult<()> {
    seqmap::load_shmmr_map(&mut mi.intervals, &String::from(filename.to_string_lossy()))?;
    Ok(())
}

#[pyfunction]
fn map_interval_query(
    ivtl: &MapIntervals,
    sid: u32,
    bgn: u32,
    end: u32,
) -> PyResult<Vec<MapIntervalRecord>> {
    let mut rtn = Vec::<MapIntervalRecord>::new();
    if ivtl.intervals.contains_key(&sid) {
        let mut q_res: Vec<_> = ivtl
            .intervals
            .get(&sid)
            .unwrap()
            .query(bgn..end)
            .map(|x| x.value)
            .collect();
        q_res.sort();
        for e in q_res {
            rtn.push(e);
            //println!("q: {} {} {} r:{} {} {} {} {} {} {}", sid, bgn, end, e[0], e[1], e[2], e[3], e[4], e[5], e[6]);
        }
    }
    Ok(rtn)
}

#[pyfunction]
fn generate_deltas(seqdb0: &SeqDB, seqdb1: &SeqDB, m: MapIntervalRecord, k: u32) -> PyResult<()> {
    let mut rev = false;
    if m[6] == 1 {
        //strand
        rev = true;
    }
    let sid0 = m[0] as usize;
    let b0 = m[1] as usize;
    let e0 = m[2] as usize;
    let sid1 = m[3] as usize;
    let b1 = m[4] as usize;
    let e1 = m[5] as usize;
    let subseq0;
    let subseq1;
    let strand: u8;
    let bb0;
    let ee0;
    let bb1;
    let ee1;
    // we need to pad the sequence with the minimize on both end. the match_reads() need at least 8 base match at the beginning
    match rev {
        true => {
            bb0 = b0 - (k - 1) as usize;
            ee0 = e0;
            subseq0 = seqdb0.seqs[sid0 as usize][bb0..ee0].to_vec();

            bb1 = b1 - (k - 2) as usize;
            ee1 = e1 + 1;

            let subseq_tmp = seqdb1.seqs[sid1 as usize][bb1..ee1].to_vec();
            subseq1 = seqmap::rc_dna_seq(subseq_tmp);
            strand = 1;
            /*
            println!("DEBUB s0 {}", String::from_utf8(subseq0.clone()).unwrap());
            println!("DEBUB s1 {}", String::from_utf8(subseq1.clone()).unwrap());
            */
        }
        false => {
            bb0 = b0 - (k - 1) as usize;
            ee0 = e0;
            subseq0 = seqdb0.seqs[sid0 as usize][bb0..ee0].to_vec();
            bb1 = b1 - (k - 1) as usize;
            ee1 = e1;
            subseq1 = seqdb1.seqs[sid1 as usize][bb1..ee1].to_vec();
            strand = 0;
        }
    }

    let v = seqmap::get_deltas(
        &subseq0,
        &subseq1,
        sid0 as u32,
        sid1 as u32,
        bb0 as u32,
        bb1 as u32,
        strand,
    );

    let sname0 = &seqdb0.id2name[&(sid0 as u32)];
    let sname1 = &seqdb1.id2name[&(sid1 as u32)];
    println!(
        "M {} {} {} {} {} {} {} {} {} {}",
        sname0,
        sid0,
        b0,
        e0,
        sname1,
        sid1,
        b1,
        e1,
        strand,
        v.len()
    );
    /*
    v.iter().for_each(|e| {
        println!(
            "D {} {} {} {} {} {} {} {} {}",
            sname0, e.sid0, e.pos0, sname1, e.sid1, e.pos1, e.strand1, e.dk, e.base
        );
    });
    */

    Ok(())
}

#[pyfunction]
fn generate_deltas_from_read_db(
    seqdb0: &SeqDB,
    seqdb1: &ReadDB,
    m: MapIntervalRecord,
) -> PyResult<()> {
    let mut rev = false;
    if m[6] == 1 {
        rev = true;
    }
    let sid0 = m[0] as usize;
    let b0 = m[1] as usize;
    let e0 = m[2] as usize;
    let subseq0 = seqdb0.seqs[sid0 as usize][b0..e0].to_vec();

    let sid1 = m[3];
    let b1 = m[4] as usize;
    let e1 = m[5] as usize;

    let sname = seqdb1.id2name.get(&sid1).unwrap();
    //let base_map = &[b'\0', b'A', b'C', b'\0', b'G', b'\0', b'\0', b'\0', b'T'];
    let base_map = &[b'A', b'C', b'G', b'T'];
    let s = *seqdb1.starts.get(sname).unwrap();

    let subseq1;
    let strand: u8;
    match rev {
        true => {
            let bb1 = b1 - (56 - 1); //kmer size - 1
            let ee1 = e1 - (56 - 1);
            let subseq_tmp = seqdb1.seqs[s + bb1..s + ee1]
                .iter()
                .map(|c| base_map[(c & 0x0011) as usize])
                .collect::<Vec<u8>>();
            subseq1 = seqmap::rc_dna_seq(subseq_tmp);
            strand = 1;
        }
        false => {
            subseq1 = seqdb1.seqs[s + b1..s + e1]
                .iter()
                .map(|c| base_map[(c & 0b0011) as usize])
                .collect::<Vec<u8>>();
            strand = 0;
        }
    }
    let v = seqmap::get_deltas(
        &subseq0,
        &subseq1,
        sid0 as u32,
        sid1 as u32,
        b0 as u32,
        b1 as u32,
        strand,
    );

    let sname0 = &seqdb0.id2name[&(sid0 as u32)];
    let sname1 = &seqdb1.id2name[&(sid1 as u32)];
    println!(
        "M {} {} {} {} {} {} {} {} {} {}",
        sname0,
        sid0,
        b0,
        e0,
        sname1,
        sid1,
        b1,
        e1,
        strand,
        v.len()
    );
    /*
    v.iter().for_each(|e| {
        println!(
            "D {} {} {} {} {} {} {} {} {}",
            sname0, e.sid0, e.pos0, sname1, e.sid1, e.pos1, e.strand1, e.dk, e.base
        );
    });
    */
    Ok(())
}

#[pyfunction]
fn find_match_chain(mqr: Vec<MapIntervalRecord>, max_count:Option<u32>) -> PyResult<Vec<MapIntervalRecord>> {
    let matches = &mqr;
    let mqr = seqmap::find_match_chain(matches, max_count);
    Ok(mqr)
}


#[pyfunction]
fn find_match_chain_m(mqr: Vec<MapIntervalRecord>) -> PyResult<Vec<MapIntervalRecord>> {
    let matches = &mqr;
    let mqr = multi_seqmap::find_match_chain(matches);
    Ok(mqr)
}

#[pyfunction]
fn get_cigar(seq0: &PyString, seq1: &PyString) -> PyResult<(isize, String, Vec<u8>)> {
    let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);
    let pattern = seq0.to_string();
    let text = seq1.to_string();
    let mut penalties = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };
    let pat_len = pattern.as_bytes().len();
    let text_len = text.as_bytes().len();

    let mut wavefronts =
        AffineWavefronts::new_reduced(pat_len, text_len, &mut penalties, 100, 100, &alloc);
    wavefronts
        .align(pattern.as_bytes(), text.as_bytes())
        .unwrap();

    let score = wavefronts.edit_cigar_score(&mut penalties);
    //let cigar = wavefronts.cigar_bytes_raw();
    let cigar = wavefronts.cigar_bytes();
    let cg_str = std::str::from_utf8(&cigar).unwrap();

    Ok((score, cg_str.to_string(), wavefronts.cigar_bytes_raw()))
}

#[pyfunction]
fn map_seqs_with_db(
    seq0_id_length: Vec<[u32;2]>,
    shmmrmap: &MapIntervals,
) -> PyResult<Vec<MapIntervalRecord>> {
    Ok(seqmap::map_seqs_with_db(seq0_id_length, &shmmrmap.intervals))
}

#[pyfunction]
fn get_aln_segements(
    ref_id: u32,
    ref_seq: &PyString,
    tgt_id: u32,
    tgt_seq: &PyString,
) -> PyResult<Vec<AlnSegment>> {
    let ref_seq = ref_seq.to_string();
    let tgt_seq = tgt_seq.to_string();
    let aln_segs = seqs2variants::get_aln_segements(ref_id, &ref_seq, tgt_id, &tgt_seq);

    match aln_segs {
        Ok(segs) => Ok(segs
            .iter()
            .map(|seg| {
                let t = match seg.t {
                    seqs2variants::AlnSegType::Match => b'M',
                    seqs2variants::AlnSegType::Mismatch => b'X',
                    seqs2variants::AlnSegType::Insertion => b'I',
                    seqs2variants::AlnSegType::Deletion => b'D',
                    seqs2variants::AlnSegType::Unspecified => b'?',
                };
                AlnSegment {
                    t: t,
                    ref_loc: (seg.ref_loc.id, seg.ref_loc.bgn, seg.ref_loc.len),
                    tgt_loc: (seg.tgt_loc.id, seg.tgt_loc.bgn, seg.tgt_loc.len),
                }
            })
            .collect()),
        Err(_) => Err(exceptions::PyException::new_err("alignment failed")),
    }
}

#[pyfunction]
fn get_aln_map(aln_segs: Vec<AlnSegment>, s0: &PyString, s1: &PyString) -> PyResult<AlnMap> {
    let s0 = s0.to_string();
    let s1 = s1.to_string();
    let aln_segs = aln_segs
        .iter()
        .map(|s| seqs2variants::AlnSegment {
            ref_loc: seqs2variants::SeqLocus {
                id: s.ref_loc.0,
                bgn: s.ref_loc.1,
                len: s.ref_loc.2,
            },
            tgt_loc: seqs2variants::SeqLocus {
                id: s.tgt_loc.0,
                bgn: s.tgt_loc.1,
                len: s.tgt_loc.2,
            },
            t: match s.t {
                b'M' => seqs2variants::AlnSegType::Match,
                b'X' => seqs2variants::AlnSegType::Mismatch,
                b'I' => seqs2variants::AlnSegType::Insertion,
                b'D' => seqs2variants::AlnSegType::Deletion,
                _ => seqs2variants::AlnSegType::Unspecified,
            },
        })
        .collect::<seqs2variants::AlnSegments>();

    let aln_map = seqs2variants::get_aln_map(&aln_segs, &s0, &s1).unwrap();

    Ok(AlnMap {
        pmap: aln_map.pmap,
        ref_a_seq: aln_map.ref_a_seq,
        tgt_a_seq: aln_map.tgt_a_seq,
        aln_seq: aln_map.aln_seq,
    })
}

#[pymodule]
fn pgr(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<SeqDB>()?;
    m.add_class::<ReadDB>()?;
    m.add_class::<MapIntervals>()?;
    m.add_function(wrap_pyfunction!(get_shmmr_dots, m)?)?;
    m.add_function(wrap_pyfunction!(generate_shmmr_map, m)?)?;
    m.add_function(wrap_pyfunction!(generate_shmmr_map_from_read_db, m)?)?;
    m.add_function(wrap_pyfunction!(map_interval_query, m)?)?;
    m.add_function(wrap_pyfunction!(write_shmmr_map, m)?)?;
    m.add_function(wrap_pyfunction!(load_shmmr_map, m)?)?;
    m.add_function(wrap_pyfunction!(generate_deltas, m)?)?;
    m.add_function(wrap_pyfunction!(generate_deltas_from_read_db, m)?)?;
    m.add_function(wrap_pyfunction!(find_match_chain, m)?)?;
    m.add_function(wrap_pyfunction!(find_match_chain_m, m)?)?;
    m.add_function(wrap_pyfunction!(build_shmmer_map_from_query_results, m)?)?;
    m.add_function(wrap_pyfunction!(get_cigar, m)?)?;
    m.add_function(wrap_pyfunction!(get_aln_segements, m)?)?;
    m.add_function(wrap_pyfunction!(get_aln_map, m)?)?;
    m.add_function(wrap_pyfunction!(map_seqs_with_db, m)?)?;
    Ok(())
}
