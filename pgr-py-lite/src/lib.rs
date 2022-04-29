// src/lib.rs
use pgr_db::agc_io;
use pgr_db::aln::{self, HitPair};
use pgr_db::seq_db;
use pyo3::exceptions;
// use pgr_utils::fasta_io;
use pgr_db::shmmrutils::{sequence_to_shmmrs, ShmmrSpec};
// use pyo3::exceptions;
use pyo3::prelude::*;
// use pyo3::types::PyString;
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
use pgr_utils::seqs2variants;
use pyo3::types::PyString;
use pyo3::wrap_pyfunction;
use pyo3::Python;
use rustc_hash::FxHashMap;

#[pyclass]
#[derive(Clone)]
struct ShmmrFragMap {
    pub shmmr_spec: Option<ShmmrSpec>,
    #[pyo3(get)]
    pub shmmr_to_frags: seq_db::ShmmrToFrags,
}

#[pymethods]
impl ShmmrFragMap {
    #[new]
    pub fn new() -> Self {
        let shmmr_to_frags = seq_db::ShmmrToFrags::default();
        let shmmr_spec = None;
        ShmmrFragMap {
            shmmr_spec,
            shmmr_to_frags,
        }
    }

    pub fn load_from_mdb(&mut self, filename: String) -> () {
        let (shmmr_spec, new_map) = seq_db::read_mdb_file(filename).unwrap();
        self.shmmr_to_frags = new_map;
        self.shmmr_spec = Some(shmmr_spec);
    }

    pub fn load_from_fastx(
        &mut self,
        filepath: String,
        w: u32,
        k: u32,
        r: u32,
        min_span: u32,
    ) -> PyResult<()> {
        let spec = ShmmrSpec {
            w,
            k,
            r,
            min_span,
            sketch: false,
        };
        let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
        sdb.load_index_from_fastx(filepath)?;
        self.shmmr_to_frags = sdb.frag_map;
        self.shmmr_spec = Some(spec);
        Ok(())
    }

    pub fn load_from_seq_list(
        &mut self,
        source: Option<String>,
        seq_list: Vec<(u32, String, Vec<u8>)>,
        w: u32,
        k: u32,
        r: u32,
        min_span: u32,
    ) -> PyResult<()> {
        let spec = ShmmrSpec {
            w,
            k,
            r,
            min_span,
            sketch: false,
        };
        let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
        let seq_vec = seq_list
            .into_iter()
            .map(|v| (v.0, source.clone(), v.1, v.2))
            .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
        sdb.load_index_from_seq_vec(&seq_vec);
        self.shmmr_to_frags = sdb.frag_map;
        self.shmmr_spec = Some(spec);
        Ok(())
    }

    pub fn query_fragment(
        &self,
        seq: Vec<u8>,
    ) -> PyResult<Vec<((u64, u64), (u32, u32, u8), Vec<seq_db::FragmentSignature>)>> {
        let shmmr_spec = &self.shmmr_spec.as_ref().unwrap();
        let res: Vec<((u64, u64), (u32, u32, u8), Vec<seq_db::FragmentSignature>)> =
            seq_db::query_fragment(&self.shmmr_to_frags, &seq, shmmr_spec);
        Ok(res)
    }

    pub fn query_fragment_to_hps(
        &self,
        seq: Vec<u8>,
        penality: f32,
    ) -> PyResult<Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>> {
        let shmmr_spec = &self.shmmr_spec.as_ref().unwrap();
        let res = aln::query_fragment_to_hps(&self.shmmr_to_frags, &seq, shmmr_spec, penality);
        Ok(res)
    }

    pub fn get_shmmr_spec(&self) -> PyResult<Option<(u32, u32, u32, u32, bool)>> {
        if let Some(spec) = self.shmmr_spec.as_ref() {
            Ok(Some((spec.w, spec.k, spec.r, spec.min_span, spec.sketch)))
        } else {
            Ok(None)
        }
    }

    pub fn get_shmmr_pair_list(&mut self) -> PyResult<Vec<(u64, u64, u32, u32, u32, u8)>> {
        let py_out = self
            .shmmr_to_frags
            .iter()
            .flat_map(|v| {
                v.1.iter()
                    .map(|vv| (v.0 .0, v.0 .1, vv.1, vv.2, vv.3, vv.4))
                    .collect::<Vec<(u64, u64, u32, u32, u32, u8)>>()
            })
            .collect::<Vec<(u64, u64, u32, u32, u32, u8)>>();
        Ok(py_out)
    }
}

#[pyclass(unsendable)] // lock in one thread (see https://github.com/PyO3/pyo3/blob/main/guide/src/class.md)
struct AGCFile {
    agc_file: agc_io::AGCFile,
    #[pyo3(get)]
    pub ctg_lens: FxHashMap<(String, String), usize>,
}

#[pymethods]
impl AGCFile {
    #[new]
    pub fn new(filepath: String) -> Self {
        let agc_file = agc_io::AGCFile::new(filepath);
        let mut ctg_lens = FxHashMap::<(String, String), usize>::default();
        agc_file.ctg_lens.iter().for_each(|(k, v)| {
            ctg_lens.insert((k.0.clone(), k.1.clone()), *v);
        });
        AGCFile { agc_file, ctg_lens }
    }

    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> PyResult<Vec<u8>> {
        Ok(self.agc_file.get_sub_seq(sample_name, ctg_name, bgn, end))
    }

    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> PyResult<Vec<u8>> {
        Ok(self.agc_file.get_seq(sample_name, ctg_name))
    }
}

#[pyfunction]

pub fn sparse_aln(
    sp_hits: Vec<HitPair>,
    max_span: u32,
    penality: f32,
) -> PyResult<Vec<(f32, Vec<HitPair>)>> {
    let mut hp = sp_hits.clone();
    Ok(aln::sparse_aln(&mut hp, max_span, penality))
}

#[pyfunction]
fn get_shmmr_dots(
    seq0: Vec<u8>,
    seq1: Vec<u8>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
) -> (Vec<u32>, Vec<u32>) {
    let mut x = Vec::<u32>::new();
    let mut y = Vec::<u32>::new();
    //let seq0v = seq0.to_string_lossy().as_bytes().to_vec();
    //let seq1v = seq1.to_string_lossy().as_bytes().to_vec();
    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };

    let shmmr0 = sequence_to_shmmrs(0, &seq0, &shmmr_spec);
    let shmmr1 = sequence_to_shmmrs(1, &seq1, &shmmr_spec);
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

fn pgrlite(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<ShmmrFragMap>()?;
    m.add_class::<AGCFile>()?;
    m.add_function(wrap_pyfunction!(sparse_aln, m)?)?;
    m.add_function(wrap_pyfunction!(get_shmmr_dots, m)?)?;
    m.add_function(wrap_pyfunction!(get_cigar, m)?)?;
    m.add_function(wrap_pyfunction!(get_aln_segements, m)?)?;
    m.add_function(wrap_pyfunction!(get_aln_map, m)?)?;
    Ok(())
}
