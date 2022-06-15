// src/lib.rs
pub const VERSION_STRING: &'static str = env!("VERSION_STRING");

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
use pgr_db::seqs2variants;
use pyo3::types::PyString;
use pyo3::wrap_pyfunction;
use pyo3::Python;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Get the revision (git-hashtag) of the build
#[pyfunction]
pub fn pgr_lib_version() -> PyResult<String> {
    Ok(VERSION_STRING.to_string())
}

/// A class that stores pangenomics indices and seqeunces with multiple backend storage options (AGC, fasta file, memory)
/// Large set of genomic sequenes, a user should use AGC backend. A binary file provides the command ``pgr-mdb``
/// which can read an AGC to create the index file. For example, we can create the index files from an AGC file::
///
///     # create a file that contains a list of file that contains a set of files from which we want to build the indices
///  
///     $ echo HPRC-y1-rebuild-04252022.agc > filelist
///  
///     # using pgr-mdb to create the index files, for 97 haplotyed genome assembly from HPRC year one release,
///     # it takes about 30 to 40 min to create the index files
///
///     $ pgr-mdb filelist HPRC-y1-rebuild-04252022
///
///     # two index files will be created by the pgr-mdb command
///     # one with a suffix .mdb and another one with a suffix .midx
///     # when we use the load_from_agc_index() method, all three files, e.g., genomes.agc, genomes.mdb and
///     # genomes.midx should have the same prefix as the parameter used to call  load_from_agc_index() method
///
/// One can also create index and load the seqeunces from a fasta file using ```load_from_fastx()``` methods.
/// Curretnly, this might be a good option for mid-size dataset (up to a couple of hundred magebases).
///
/// Or, a user can load the sequnece from memory using a Python list. This is convinient when one needs to
/// rebuild the SHIMMER index with different parameters for a different resolution.
///  
/// Once the index is built, the database can be queried quickly by using the ``query_fragment()`` or
/// the ``query_fragment_to_hps()`` method.
///  
#[pyclass]
#[derive(Clone)]
struct SeqIndexDB {
    /// Rust internal: store the specification of the shmmr specifcation
    pub shmmr_spec: Option<ShmmrSpec>,
    /// Rust internal: store the sequences
    pub seq_db: Option<seq_db::CompactSeqDB>,
    /// Rust internal: store the agc file and the index
    pub agc_db: Option<(agc_io::AGCFile, seq_db::ShmmrToFrags)>,
    /// a dictionary maps (ctg_name, source) -> (id, len)
    #[pyo3(get)]
    pub seq_index: Option<HashMap<(String, Option<String>), (u32, u32)>>,
    /// a dictionary maps id -> (ctg_name, source, len)
    #[pyo3(get)]
    pub seq_info: Option<HashMap<u32, (String, Option<String>, u32)>>,
}

#[pymethods]
impl SeqIndexDB {
    /// constructor, take no argument
    #[new]
    pub fn new() -> Self {
        SeqIndexDB {
            seq_db: None,
            agc_db: None,
            shmmr_spec: None,
            seq_index: None,
            seq_info: None,
        }
    }

    /// use AGC file for sequences and the index created from an AGC file
    ///
    /// Parameters
    /// ----------
    ///
    /// prefix: string
    ///     the prefix to the `.agc`, `.mdb` and `.midx` files
    ///
    /// Returns
    /// -------
    ///
    /// None or I/O Error
    ///
    #[pyo3(text_signature = "($self, prefix)")]
    pub fn load_from_agc_index(&mut self, prefix: String) -> PyResult<()> {
        // let (shmmr_spec, new_map) = seq_db::read_mdb_file(prefix.to_string() + ".mdb").unwrap();
        let (shmmr_spec, new_map) =
            seq_db::read_mdb_file_parallel(prefix.to_string() + ".mdb").unwrap();
        let agc_file = agc_io::AGCFile::new(prefix.to_string() + ".agc")?;
        self.agc_db = Some((agc_file, new_map));
        self.seq_db = None;
        self.shmmr_spec = Some(shmmr_spec);

        let mut seq_index = HashMap::<(String, Option<String>), (u32, u32)>::new();
        let mut seq_info = HashMap::<u32, (String, Option<String>, u32)>::new();
        let midx_file = BufReader::new(File::open(prefix.to_string() + ".midx")?);
        midx_file
            .lines()
            .into_iter()
            .try_for_each(|line| -> Result<(), std::io::Error> {
                let line = line.unwrap();
                let mut line = line.as_str().split("\t");
                let sid = line.next().unwrap().parse::<u32>().unwrap();
                let len = line.next().unwrap().parse::<u32>().unwrap();
                let ctg_name = line.next().unwrap().to_string();
                let source = line.next().unwrap().to_string();
                seq_index.insert((ctg_name.clone(), Some(source.clone())), (sid, len));
                seq_info.insert(sid, (ctg_name, Some(source), len));
                Ok(())
            })?;
        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        Ok(())
    }
    /// load and create the index created from a fasta / fastq file
    ///
    /// Parameters
    /// ----------
    ///
    ///filepath : string
    ///     the path the fasta or fastq file
    ///
    /// w : int
    ///     the window size of the shimmer index, default to 80
    ///
    /// k : int
    ///     the k-mer size of the shimmer index, default to 56
    ///
    /// r : int
    ///     the reduction factor of the shimmer index, default to 4
    ///
    /// min_span : int
    ///     the min_span ofr the shimmer index, default to 8
    ///
    /// Returns
    /// -------
    ///
    /// None or I/O Error
    ///     None
    ///
    #[pyo3(text_signature = "($self, w, k, r, min_span)")]
    #[args(w = "80", k = "56", r = "4", min_span = "8")]
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
        sdb.load_seqs_from_fastx(filepath)?;
        self.shmmr_spec = Some(spec);
        let mut seq_index = HashMap::<(String, Option<String>), (u32, u32)>::new();
        let mut seq_info = HashMap::<u32, (String, Option<String>, u32)>::new();
        sdb.seqs.iter().for_each(|v| {
            seq_index.insert((v.name.clone(), v.source.clone()), (v.id, v.len as u32));
            seq_info.insert(v.id, (v.name.clone(), v.source.clone(), v.len as u32));
        });
        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        self.seq_db = Some(sdb);
        self.agc_db = None;
        Ok(())
    }

    /// load and create the index created from a python list
    ///
    /// Parameters
    /// ----------
    ///
    /// seq_list : list
    ///     a list of tuple of the form (squence_id : int, sequence_name : string, sequence: list of bytes)
    ///
    /// source : string
    ///     a string indicating the source of the sequence, default to "Memory"
    ///
    /// w : int
    ///     the window size of the shimmer index, default to 80
    ///
    /// k : int
    ///     the k-mer size of the shimmer index, default to 56
    ///
    /// r : int
    ///     the reduction factor of the shimmer index, default to 4
    ///
    /// min_span : int
    ///     the min_span ofr the shimmer index, default to 8
    ///
    /// Returns
    /// -------
    ///
    /// None or I/O Error
    ///     None
    ///
    #[pyo3(text_signature = "($self, seq_list, source, w, k, r, min_span)")]
    #[args(source = "\"Memory\"", w = "80", k = "56", r = "4", min_span = "8")]
    pub fn load_from_seq_list(
        &mut self,
        seq_list: Vec<(String, Vec<u8>)>,
        source: Option<&str>,
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
        let source = Some(source.unwrap().to_string());
        let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
        let seq_vec = seq_list
            .into_iter()
            .enumerate()
            .map(|(sid, v)| (sid as u32, source.clone(), v.0, v.1))
            .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
        sdb.load_seqs_from_seq_vec(&seq_vec);

        self.shmmr_spec = Some(spec);
        let mut seq_index = HashMap::<(String, Option<String>), (u32, u32)>::new();
        let mut seq_info = HashMap::<u32, (String, Option<String>, u32)>::new();
        sdb.seqs.iter().for_each(|v| {
            seq_index.insert((v.name.clone(), v.source.clone()), (v.id, v.len as u32));
            seq_info.insert(v.id, (v.name.clone(), v.source.clone(), v.len as u32));
        });
        self.seq_index = Some(seq_index);
        self.seq_info = Some(seq_info);
        self.seq_db = Some(sdb);
        self.agc_db = None;
        Ok(())
    }

    /// use a fragement of sequence to query the database to get all hits
    ///
    /// Parameters
    /// ----------
    ///
    /// seq : list
    ///     the sequnece in bytes used for query
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///   a list of hits in the format (shimmer_pair, query_fragement, target_fragments), where
    ///     - shimmer_pair: (int, int), tuple of the shimmer_pair
    ///     - query_fragment: (int, int, int) = (start_coordinate, end_coordinate, orientation)
    ///     - target_fragments: a list of ``FragmentSignature``: (frg_id, seq_id, bgn, end,
    ///       orientation(to the shimmer pair)) defined as::
    ///
    ///           pub type FragmentSignature = (u32, u32, u32, u32, u8);
    ///       
    ///
    #[pyo3(text_signature = "($self, seq)")]
    pub fn query_fragment(
        &self,
        seq: Vec<u8>,
    ) -> PyResult<Vec<((u64, u64), (u32, u32, u8), Vec<seq_db::FragmentSignature>)>> {
        let shmmr_spec = &self.shmmr_spec.as_ref().unwrap();
        let shmmr_to_frags = self.get_shmmr_map_internal();
        let res: Vec<((u64, u64), (u32, u32, u8), Vec<seq_db::FragmentSignature>)> =
            seq_db::query_fragment(shmmr_to_frags, &seq, shmmr_spec);
        Ok(res)
    }

    /// use a fragement of sequence to query the database to get all hits
    ///
    /// sparese dynamic programming is performed to long chain of alignment
    ///  
    /// Parameters
    /// ----------
    /// seq : list of bytes
    ///    a list of bytes representing the DNA sequence
    ///
    /// penality : float
    ///    the gap penalty factor used in sparse dyanmic programming for finding the hits
    ///
    /// merge_range_tol : int
    ///    a parameter used to merge the alignment ranges
    ///
    /// max_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` for sparse dynamic programming
    ///
    /// max_query_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` in the query sequence for sparse dynamic programming
    ///
    /// max_query_count : int
    ///    only use the shimmer pairs that less than the ``max_count`` in the target sequence for sparse dynamic programming
    ///
    /// max_aln_span : int
    ///    the size of span used in the sparse dynamic alignment for finding the hits
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     a list of tuples of
    ///     (``target_sequence_id``, (``score``, ``list_of_the_hit_pairs``)), where
    ///     the ``list_of_the_hit_pairs`` is a list of tuples of
    ///     ((``query_start``, ``query_end``, ``query_orientation``),
    ///     (``target_start``, ``target_end``, ``target_orientation``))
    #[pyo3(
        text_signature = "($self, seq, penality, max_count, max_query_count, max_target_count, max_aln_span)"
    )]
    pub fn query_fragment_to_hps(
        &self,
        seq: Vec<u8>,
        penality: f32,
        max_count: Option<u32>,
        max_count_query: Option<u32>,
        max_count_target: Option<u32>,
        max_aln_span: Option<u32>,
    ) -> PyResult<Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>> {
        let shmmr_spec = &self.shmmr_spec.as_ref().unwrap();
        let shmmr_to_frags = self.get_shmmr_map_internal();
        let res = aln::query_fragment_to_hps(
            shmmr_to_frags,
            &seq,
            shmmr_spec,
            penality,
            max_count,
            max_count_query,
            max_count_target,
            max_aln_span,
        );
        Ok(res)
    }

    /// count the number of shimmer hits in the database
    ///
    /// Parameters
    /// ----------
    ///
    /// shmmr_pair : tuple
    ///     a shimmer pair used for query
    ///
    /// Returns
    /// -------
    ///
    /// int
    ///     number of hits
    #[pyo3(text_signature = "($self, shmmr_pair)")]
    pub fn get_shmmr_pair_count(&self, shmmr_pair: (u64, u64)) -> usize {
        let shmmr_to_frags = self.get_shmmr_map_internal();
        if shmmr_to_frags.contains_key(&shmmr_pair) {
            shmmr_to_frags.get(&shmmr_pair).unwrap().len()
        } else {
            0
        }
    }

    /// count the number of shimmer hits paritioned by the source file in the database
    ///
    /// Parameters
    /// ----------
    ///
    /// shmmr_pair : tuple
    ///     a shimmer pair used for query
    ///
    /// max_unique_count : int
    ///     a interger to filter out shimmer pairs with count that are greater
    ///     than the `max_unique_count`  
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     a list of the tuple (soure_name : string, count : int)
    ///
    #[pyo3(text_signature = "($self, shmmr_pair, max_unique_count)")]
    #[args(max_unique_count = "1")]
    pub fn get_shmmr_pair_source_count(
        &self,
        shmmr_pair: (u64, u64),
        max_unique_count: Option<usize>,
    ) -> Vec<(String, usize)> {
        let shmmr_to_frags = self.get_shmmr_map_internal();
        let mut count = FxHashMap::<String, usize>::default();
        if shmmr_to_frags.contains_key(&shmmr_pair) {
            shmmr_to_frags
                .get(&shmmr_pair)
                .unwrap()
                .iter()
                .for_each(|v| {
                    let sid = v.1;
                    let source = self
                        .seq_info
                        .as_ref()
                        .unwrap()
                        .get(&sid)
                        .unwrap()
                        .1
                        .as_ref()
                        .unwrap_or(&"".to_string())
                        .clone();
                    *count.entry(source).or_insert(0) += 1;
                });
            count
                .into_par_iter()
                .filter(|(_k, v)| {
                    if let Some(muc) = max_unique_count {
                        if *v > muc {
                            false
                        } else {
                            true
                        }
                    } else {
                        true
                    }
                })
                .map(|(k, v)| (k, v))
                .collect::<Vec<(String, usize)>>()
        } else {
            vec![]
        }
    }

    /// Output the specific of the shimmer used to build the index
    ///
    /// Returns
    /// -------
    ///
    /// tuple
    ///     (window_size, k_mer_size, reduction_factor, min_space, use_sketch)
    ///
    pub fn get_shmmr_spec(&self) -> PyResult<Option<(u32, u32, u32, u32, bool)>> {
        if let Some(spec) = self.shmmr_spec.as_ref() {
            Ok(Some((spec.w, spec.k, spec.r, spec.min_span, spec.sketch)))
        } else {
            Ok(None)
        }
    }

    /// get the ``shmmer_pair`` to ``fragment_id`` map in Python
    ///
    /// this can be very expensive to generate the Python objects of a large hashmap in Rust
    ///
    /// Returns
    /// -------
    ///
    /// dict
    ///     the ``shmmer_pair`` to ``fragments`` map
    ///
    ///     fragments: a list of ``FragmentSignature``: (frg_id, seq_id, bgn, end,
    ///     orientation(to the shimmer pair)) defined as::
    ///
    ///         pub type FragmentSignature = (u32, u32, u32, u32, u8);
    ///
    pub fn get_shmmr_map(&self) -> PyResult<PyObject> {
        // very expansive as the Rust FxHashMap will be converted to Python's dictionary
        // maybe limit the size that can be converted to avoid OOM
        let shmmr_to_frags = self.get_shmmr_map_internal();
        pyo3::Python::with_gil(|py| Ok(shmmr_to_frags.to_object(py)))
    }

    /// get the ``shmmer_pair`` to ``fragment_id`` map in Python as a list
    ///
    /// this can be very expensive to generate the Python objects of a large hashmap in Rust
    ///
    /// Returns
    /// -------
    ///
    /// list
    ///     list of the tuple (shmmr0, shmmr1, seq_id, position0, position1, orientation)
    ///   
    pub fn get_shmmr_pair_list(&mut self) -> PyResult<Vec<(u64, u64, u32, u32, u32, u8)>> {
        let shmmr_to_frags = self.get_shmmr_map_internal();
        let py_out = shmmr_to_frags
            .par_iter()
            .flat_map(|v| {
                v.1.iter()
                    .map(|vv| (v.0 .0, v.0 .1, vv.1, vv.2, vv.3, vv.4))
                    .collect::<Vec<(u64, u64, u32, u32, u32, u8)>>()
            })
            .collect::<Vec<(u64, u64, u32, u32, u32, u8)>>();
        Ok(py_out)
    }

    /// fetch a contiguous sub-sequence
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    /// bgn : int
    ///     the starting coordinate (0-based)
    /// end : int
    ///     the ending coordinate (exclusive)  
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(text_signature = "($self, sample_name, ctg_name, bgn, end)")]
    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> PyResult<Vec<u8>> {
        if self.agc_db.is_some() {
            Ok(self
                .agc_db
                .as_ref()
                .unwrap()
                .0
                .get_sub_seq(sample_name, ctg_name, bgn, end))
        } else {
            let &(sid, _) = self
                .seq_index
                .as_ref()
                .unwrap()
                .get(&(ctg_name, Some(sample_name)))
                .unwrap();
            let seq = self.seq_db.as_ref().unwrap().get_seq_by_id(sid);
            Ok(seq[bgn..end].to_vec())
        }
    }

    /// fetch a sequence
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[pyo3(text_signature = "($self, sample_name, ctg_name)")]
    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> PyResult<Vec<u8>> {
        if self.agc_db.is_some() {
            Ok(self
                .agc_db
                .as_ref()
                .unwrap()
                .0
                .get_seq(sample_name, ctg_name))
        } else {
            let &(sid, _) = self
                .seq_index
                .as_ref()
                .unwrap()
                .get(&(ctg_name, Some(sample_name)))
                .unwrap();
            Ok(self.seq_db.as_ref().unwrap().get_seq_by_id(sid))
        }
    }

    /// Get adjecent list of the shimmer graph shimmer_pair -> shimmer_pair
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     the minimum number of times a pair of shimmers must be observed to be included in the graph
    ///
    /// Returns
    /// -------
    /// list
    ///     list of pairs of shimmer pairs ((h00, h01, orientation0),(h10, h11, orientation1))  
    ///
    pub fn get_smp_adj_list(&self, min_count: usize) -> Vec<(u32, (u64, u64, u8), (u64, u64, u8))> {
        let frag_map = self.get_shmmr_map_internal();
        seq_db::frag_map_to_adj_list(frag_map, min_count)
    }

    /// Sort the adjecent list of the shimmer graph
    ///
    /// Parameters
    /// ----------
    /// adj_list : list
    ///     the list from the output from ``get_smp_adj_list``
    ///
    /// sort_by : (u64, u64)
    ///     the starting node signature
    ///
    /// Returns
    /// -------
    /// list
    ///     list of node (node.0, node.1, node_count, is_leaf)
    ///
    ///

    pub fn sort_adj_list_by_weighted_dfs(
        &self,
        adj_list: Vec<(u32, (u64, u64, u8), (u64, u64, u8))>,
        start: (u64, u64, u8),
    ) -> Vec<(
        (u64, u64, u8),
        Option<(u64, u64, u8)>,
        u32,
        bool,
        u32,
        u32,
        u32,
    )> {
        let frag_map = self.get_shmmr_map_internal();
        seq_db::sort_adj_list_by_weighted_dfs(&frag_map, &adj_list, start)
    }

    /// Convert the adjecent list of the shimmer graph shimmer_pair -> GFA
    ///
    /// Parameters
    /// ----------
    /// min_count : int
    ///     the minimum number of times a pair of shimmers must be observed to be included in the graph
    ///
    /// filenpath : string
    ///     the path to the output file
    ///
    /// Returns
    /// -------
    ///
    ///     the adj_list is written to a file in GFA v.1 format
    ///
    pub fn generate_smp_gfa(&self, min_count: usize, filepath: &str) -> PyResult<()> {
        let frag_map = self.get_shmmr_map_internal();
        let adj_list = seq_db::frag_map_to_adj_list(frag_map, min_count);
        let mut overlaps =
            FxHashMap::<((u64, u64, u8), (u64, u64, u8)), Vec<(u32, u8, u8)>>::default();
        let mut frag_id = FxHashMap::<(u64, u64), usize>::default();
        let mut id = 0_usize;
        adj_list.iter().for_each(|(k, v, w)| {
            if v.0 <= w.0 {
                let key = (*v, *w);
                let val = (*k, v.2, w.2);
                overlaps.entry(key).or_insert_with(|| vec![]).push(val);
                frag_id.entry((v.0, v.1)).or_insert_with(|| {
                    let c_id = id;
                    id += 1;
                    c_id
                });
                frag_id.entry((w.0, w.1)).or_insert_with(|| {
                    let c_id = id;
                    id += 1;
                    c_id
                });
            }
        });

        let mut out_file = BufWriter::new(File::create(filepath).unwrap());

        out_file.write("H\tVN:Z:1.0\tCM:Z:Sparse Genome Graph Generated By pgr-tk\n".as_bytes())?;
        (&frag_id)
            .into_iter()
            .try_for_each(|(smp, id)| -> PyResult<()> {
                let hits = frag_map.get(&smp).unwrap();
                let ave_len =
                    hits.iter().fold(0_u32, |len_sum, &s| len_sum + s.3 - s.2) / hits.len() as u32;
                let seg_line = format!(
                    "S\t{}\t*\tLN:i:{}\tSN:Z:{:016x}_{:016x}\n",
                    id, ave_len, smp.0, smp.1
                );
                out_file.write(seg_line.as_bytes())?;
                Ok(())
            })?;

        let kmer_size = self.shmmr_spec.as_ref().unwrap().k;

        overlaps
            .into_iter()
            .try_for_each(|(op, vs)| -> PyResult<()> {
                let o1 = if op.0 .2 == 0 { "+" } else { "-" };
                let o2 = if op.1 .2 == 0 { "+" } else { "-" };
                let id0 = frag_id.get(&(op.0 .0, op.0 .1)).unwrap();
                let id1 = frag_id.get(&(op.1 .0, op.1 .1)).unwrap();
                let overlap_line = format!(
                    "L\t{}\t{}\t{}\t{}\t{}M\tSC:i:{}\n",
                    id0,
                    o1,
                    id1,
                    o2,
                    kmer_size,
                    vs.len()
                );
                out_file.write(overlap_line.as_bytes())?;
                Ok(())
            })?;

        Ok(())
    }

    fn write_midx_to_text_file(&self, filepath: &str) -> Result<(), std::io::Error> {
        let mut writer = BufWriter::new(File::create(filepath)?);

        self.seq_info.as_ref().unwrap().iter().try_for_each(
            |(k, v)| -> Result<(), std::io::Error> {
                let line = format!(
                    "C\t{}\t{}\t{}\t{}\n",
                    k,
                    v.0,
                    v.1.clone().unwrap_or("NA".to_string()),
                    v.2
                );
                writer.write(line.as_bytes())?;
                Ok(())
            },
        )?;

        if let Some(shmmr_spec) = self.shmmr_spec.clone() {
            writer.write(
                format!(
                    "K\t{}\t{}\t{}\t{}\t{}\n",
                    shmmr_spec.w,
                    shmmr_spec.k,
                    shmmr_spec.r,
                    shmmr_spec.min_span,
                    shmmr_spec.sketch
                )
                .as_bytes(),
            )?;
        }
        let frag_map = self.get_shmmr_map_internal();
        frag_map
            .iter()
            .try_for_each(|v| -> Result<(), std::io::Error> {
                v.1.iter()
                    .try_for_each(|vv| -> Result<(), std::io::Error> {
                        writer.write(
                            format!(
                                "F\t{:016x}\t{:016x}\t{}\t{}\t{}\t{}\t{}\n",
                                v.0 .0, v.0 .1, vv.0, vv.1, vv.2, vv.3, vv.4
                            )
                            .as_bytes(),
                        )?;
                        Ok(())
                    })?;
                Ok(())
            })?;
        Ok(())
    }
}

impl SeqIndexDB {
    // depending on the storage type, return the corresponding index
    fn get_shmmr_map_internal(&self) -> &seq_db::ShmmrToFrags {
        let shmmr_to_frags;
        if self.agc_db.is_some() {
            shmmr_to_frags = &self.agc_db.as_ref().unwrap().1;
        } else {
            shmmr_to_frags = &self.seq_db.as_ref().unwrap().frag_map;
        }
        shmmr_to_frags
    }
}
/// A PyO3 class wrapping an exising AGC file for reading
///
/// Examaple::
///
///      >>> agc_file = AGCFile("/path/to/genomes.agc")
///  
#[pyclass(unsendable)] // lock in one thread (see https://github.com/PyO3/pyo3/blob/main/guide/src/class.md)
struct AGCFile {
    /// internal agc_file handle
    agc_file: agc_io::AGCFile,

    /// A hashmap mapping (source, ctg_name) to sequence length.
    /// It is more efficient to make a copy in Python to access it to avoid
    /// expensive Rust struct ot Python object conversion than using it directly. For example::
    ///
    ///     >>> agc_file = AGCFile("genomes.agc")
    ///     >>> agc_ctg_lens = agc_file.ctg_lens.copy()
    ///  
    #[pyo3(get)]
    pub ctg_lens: FxHashMap<(String, String), usize>,
}

#[pymethods]
impl AGCFile {
    /// constructor
    ///
    /// Parameters
    /// ----------
    /// filepath: string
    ///     the path to a AGC file
    #[args(filepath)]
    #[new]
    pub fn new(filepath: String) -> PyResult<Self> {
        let agc_file = agc_io::AGCFile::new(filepath)?;
        let mut ctg_lens = FxHashMap::<(String, String), usize>::default();
        agc_file.ctg_lens.iter().for_each(|(k, v)| {
            ctg_lens.insert((k.0.clone(), k.1.clone()), *v);
        });
        Ok(AGCFile { agc_file, ctg_lens })
    }

    /// fetch a contiguous sub-sequence from an AGC file
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    /// bgn : int
    ///     the starting coordinate (0-based)
    /// end : int
    ///     the ending coordinate (exclusive)  
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[args(sample_name, ctg_name, bgn, end)]
    #[pyo3(text_signature = "($self, sample_name, ctg_name, bgn, end)")]
    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> PyResult<Vec<u8>> {
        Ok(self.agc_file.get_sub_seq(sample_name, ctg_name, bgn, end))
    }

    /// fetch a full contig sequence from an AGC file
    ///
    /// Parameters
    /// ----------
    /// sample_name : string
    ///     the sample name stored in the AGC file
    /// ctg_name : string
    ///     the contig name stored in the AGC file
    ///
    /// Returns
    /// -------
    /// list
    ///     a list of bytes representing the sequence
    #[args(sample_name, ctg_name)]
    #[pyo3(text_signature = "($self, sample_name, ctg_name)")]
    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> PyResult<Vec<u8>> {
        Ok(self.agc_file.get_seq(sample_name, ctg_name))
    }
}

/// Perform sparse dynamic programming to identify alignment between sequence
/// using matched shimmer pairs
///
/// Parameters
/// ----------
/// sp_hits : list
///     a list of tuple of ``HitPair`` defined as
///     ``pub type HitPair = ((u32, u32, u8), (u32, u32, u8))``
///     This represents the hits as the position of matched Shimmer Pairs from
///     the two sequence. For example, if there two shimmers at positiong 2342
///     and 4322 of the query sequence that matches the shimmers at positions
///     6125465 and 6127445, them the HitPair will be ``(2342, 4322, 0)`` and
///     ``(6125465, 6127445, 0)``. The third number would be 1 if shimmer paired
///     are reversed matched to the sequence orientation and 0 otherwise.
///       
/// max_span : int
///     For a give hit, the max_span defines how many other following hits are
///     considered for the next aligned position. This will limit the search
///     space for the best alignment. If the two sequences are very repetitive,
///     then one needs to use larger ``max_span`` to ensure capturing the right
///     alignment path
///
/// penality : float
///     this parameter will determine when to break the alignment if there are big
///     gaps between the alignment segement. One can set it to zero to catch large
///     chunk alignment ignoring the gaps. Typically, a number between 0.1 to 0.5 should
///     be used.
///
#[pyfunction(sp_hits, max_span, penality)]
#[pyo3(text_signature = "($self, sp_hits, max_span, penality)")]
pub fn sparse_aln(
    sp_hits: Vec<HitPair>,
    max_span: u32,
    penality: f32,
) -> PyResult<Vec<(f32, Vec<HitPair>)>> {
    let mut hp = sp_hits.clone();
    Ok(aln::sparse_aln(&mut hp, max_span, penality))
}

/// Generate a list of shimmer pair from a sequence
///
/// Parameters
/// ----------
/// w : int
///     window size, default to 80, max allowed is 128
///
/// k : int
///     k-mer size, default to 56, max allowed is 56
///
/// r : int
///     reduction factor for generate sparse hierarchical minimiers (shimmer),
///     default to 4, max allowed is 12
///
/// min_span : int
///     - a parameter to remove close-by shimmer pairs
///     - if not zero, shimmer pairs whic the distance between them are less
///       than ``the min_span`` will be removed
///     - default to 16
///
/// padding : bool
///     - for short fragment that segemented by using shimmer, set ``padding`` to true
///       to preserve the first and last shimmers
///     - default to false
///
/// Returns
/// -------
/// list of tuple
///     a list fo tuple of ``(shimmr0, shimmr1, position0, position1, orientation)``  
///
#[pyfunction(w = "80", k = "56", r = "4", min_span = "16", padding = "false")]
#[pyo3(text_signature = "($self, w, k, r, min_span, padding)")]
fn get_shmmr_pairs_from_seq(
    seq: Vec<u8>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    padding: bool,
) -> PyResult<Vec<(u64, u64, u32, u32, u8)>> {
    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let shmmrs = sequence_to_shmmrs(0, &seq, &shmmr_spec, padding);
    let res = seq_db::pair_shmmrs(&shmmrs)
        .par_iter()
        .map(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.x >> 8;
            let s1 = s1.x >> 8;
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            }
        })
        .collect::<Vec<(u64, u64, u32, u32, u8)>>();
    Ok(res)
}

/// Generate a list of shimmer matches for creating a dot plot between two sequences
///
/// Parameters
/// ----------
///
/// seq0 : list
///     a list of bytes representing the first sequences
///
/// seq1 : list
///     a list of bytes representing the second sequences
///
/// w : int
///     window size, default to 80, max allowed is 128
///
/// k : int
///     k-mer size, default to 56, max allowed is 56
///
/// r : int
///     reduction factor for generate sparse hierarchical minimiers (shimmer),
///     default to 4, max allowed is 12
///
/// min_span : int
///     -  a parameter to remove close-by shimmer pairs
///     -  if not zero, shimmer pairs whic the distance between them are less
///        than ``the min_span`` will be removed
///
/// Returns
/// -------
/// tuple of two lists
///     ``(x, y)``:
///
///     -  ``x``: the matched shimmer positions in sequence 0
///     -  ``y``: the matched shimmer positions in sequence 1
///
#[pyfunction(w = "80", k = "56", r = "4", min_span = "16")]
#[pyo3(text_signature = "($self, seq0, seq1, w, k, r, min_span)")]
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

    let shmmr0 = sequence_to_shmmrs(0, &seq0, &shmmr_spec, false);
    let shmmr1 = sequence_to_shmmrs(1, &seq1, &shmmr_spec, false);
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

/// A wrapper class to represent alignment segement for python
///
/// This wraps the Rust struct ``seq2variants::AlnSegment`` mapping
/// the enum ``seq2variants::AlnSegType`` to intergers
///
#[pyclass]
#[derive(Clone)]
struct AlnSegment {
    /// alignment type: value =  ``b'M'``, ``b'I'``, ``b'D'``, ``b'X'``, or ``b'?'``
    #[pyo3(get, set)]
    t: u8,
    /// segment coordinate in the reference: (begin, end, length)
    #[pyo3(get, set)]
    ref_loc: (u32, u32, u32),
    /// segment coordinate in the target: (begin, end, length)
    #[pyo3(get, set)]
    tgt_loc: (u32, u32, u32),
}

/// A  class to represent alignment mapping between the two sequences
#[pyclass]
#[derive(Clone)]
pub struct AlnMap {
    /// mapped location between two sequence as a list of paired coordinate
    #[pyo3(get, set)]
    pmap: Vec<(u32, u32)>,
    /// alignment string of the reference
    #[pyo3(get, set)]
    ref_a_seq: Vec<u8>,
    /// alignment string of the target
    #[pyo3(get, set)]
    tgt_a_seq: Vec<u8>,
    /// alignment string of alignment symbols
    #[pyo3(get, set)]
    aln_seq: Vec<u8>,
}

/// Generate the CIGAR string from two sequences with WFA
///
/// Parameters
/// ----------
///
/// seq0 : string
///     a string representing the first sequence
///
/// seq1 : string
///     a string representing the second sequence
///
/// Returns
/// -------
/// tuple
///     tuple of (alignment_score, CIGAR_string, CIGAR_list)
///
#[pyfunction(seq0, seq1)]
#[pyo3(text_signature = "($self, seq0, seq1)")]
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

/// Get alignement segments from two sequences
///
/// Parameters
/// ----------
/// ref_id : int
///     a interger id for the reference sequence
///
/// ref_seq : string
///     a python string of the reference sequence
///
/// tgt_id : int
///     a interger id for the target sequnece
///
/// tgt_seq : string
///     a python string of the target sequence
///
/// Returns
/// -------
/// list
///     a list of ``AlnSegement``
///
///     the ``AlnSegment`` is a Rust struct defined as::
///
///         pub struct SeqLocus {
///             pub id: u32,
///             pub bgn: u32,
///             pub len: u32,
///         }
///
///         pub enum AlnSegType {
///             Match,
///             Mismatch,
///             Insertion,
///             Deletion,
///             Unspecified,
///         }
///
///         pub struct AlnSegment {
///             pub ref_loc: SeqLocus,
///             pub tgt_loc: SeqLocus,
///             pub t: AlnSegType,
///         }
///
#[pyfunction(ref_id, ref_seq, tgt_id, tgt_seq)]
#[pyo3(text_signature = "($self, ref_id, ref_seq, tgt_id, tgt_seq)")]
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
            .par_iter()
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

/// Get alignement map from a list of alignment segments
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the ``AlnSegment``
///
/// s0: string
///     a python string of the reference sequence
///
/// s1: int
///     a interger id for the target sequnece
///
/// Returns
/// -------
/// list
///     a list of ``AlnSegement``
#[pyfunction(aln_segs, s0, s1)]
#[pyo3(text_signature = "($self, aln_segs, s0, s1)")]
fn get_aln_map(
    aln_segs: Vec<AlnSegment>,
    ref_seq: &PyString,
    tgt_seq: &PyString,
) -> PyResult<AlnMap> {
    let s0 = ref_seq.to_string();
    let s1 = tgt_seq.to_string();
    let aln_segs = aln_segs
        .par_iter()
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

/// Perform a navie de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// kmer_size : int  
///     the size of kmers used for constructing the de Bruijn graph
///
/// min_cov : int
///     to keep hyplotype specific consensus, if a kmer has coverage more or equal to min_cov, it will be kept
///
/// Returns
/// -------
/// list
///     a list of bytes representing the consensus sequence
///
#[pyfunction(seqs, kmer_size = 33, min_cov = 2)]
#[pyo3(text_signature = "($self, seqs, kmer_size, min_cov)")]
pub fn naive_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    kmer_size: usize,
    min_cov: usize,
) -> PyResult<Vec<u8>> {
    let consensus = pgr_db::ec::naive_dbg_consensus(seqs, kmer_size, min_cov);
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// Perform a shimmer de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// k, w, r, min_span : int  
///     specification of the shimmers for construting graph
///
/// Returns
/// -------
/// list
///     a list of a set of bytes representing the consensus sequences of all branches in the graph
///
#[pyfunction(seqs, w = 33, k = 33, r = 1, min_span = 0)]
#[pyo3(text_signature = "($self, seqs, w, k, r, min_span)")]
pub fn shmmr_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
) -> PyResult<Vec<(Vec<u8>, Vec<u32>)>> {
    let spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let consensus = pgr_db::ec::shmmr_dbg_consensus(seqs, &Some(spec));
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// Perform a guided shimmer de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// k, w, r, min_span : int  
///     specification of the shimmers for construting graph
///
/// min_cov : int
///     to keep hyplotype specific consensus, if a kmer has coverage more or equal to min_cov, it will be kept
///
/// Returns
/// -------
/// list
///     a list of a set of bytes representing the consensus sequences of all branches in the graph
///
#[pyfunction(seqs, w = 33, k = 33, r = 1, min_span = 0, min_cov = 2)]
#[pyo3(text_signature = "($self, seqs, w, k, r, min_span, min_cov)")]
pub fn guided_shmmr_dbg_consensus(
    seqs: Vec<Vec<u8>>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    min_cov: u32,
) -> PyResult<(Vec<u8>, Vec<u32>)> {
    let spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let consensus = pgr_db::ec::guided_shmmr_dbg_consensus(seqs, &Some(spec), min_cov);
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// Perform a shimmer de Bruijn graph consensus
///
/// Parameters
/// ----------
/// aln_segs : list
///     a list of the list of bytes representing the bases of each sequence
///
/// k, w, r, min_span : int  
///     specification of the shimmers for construting graph
///
/// Returns
/// -------
/// list
///     a list of a set of bytes representing the consensus sequences of all branches in the graph
///
#[pyfunction(seqs, w = 33, k = 33, r = 1, min_span = 0, min_cov=2)]
#[pyo3(text_signature = "($self, seqs, w, k, r, min_span, min_cov)")]
pub fn shmmr_sparse_aln_consensus(
    seqs: Vec<Vec<u8>>,
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    min_cov: u32,
) -> PyResult<Vec<(Vec<u8>, Vec<u32>)>> {
    let spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch: false,
    };
    let consensus = pgr_db::ec::shmmr_sparse_aln_consensus(seqs, &Some(spec), min_cov);
    match consensus {
        Ok(seq) => Ok(seq),
        Err(_) => Err(exceptions::PyException::new_err(
            "consensus failed, trying bigger kmer size",
        )),
    }
}

/// The internal `pgrtk` modules implemented with Rust.
/// These classes and fucntion are re-exported as `pgrtk.*`
/// so `import pgrtk` will bring these classes and function
/// into `pgrtk.*` scope to avoid using the verbose
/// `pgrtk.pgrtk.*`.
#[pymodule]
fn pgrtk(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<SeqIndexDB>()?;
    m.add_class::<AGCFile>()?;
    m.add_function(wrap_pyfunction!(sparse_aln, m)?)?;
    m.add_function(wrap_pyfunction!(get_shmmr_dots, m)?)?;
    m.add_function(wrap_pyfunction!(get_cigar, m)?)?;
    m.add_function(wrap_pyfunction!(get_aln_segements, m)?)?;
    m.add_function(wrap_pyfunction!(get_aln_map, m)?)?;
    m.add_function(wrap_pyfunction!(pgr_lib_version, m)?)?;
    m.add_function(wrap_pyfunction!(get_shmmr_pairs_from_seq, m)?)?;
    m.add_function(wrap_pyfunction!(naive_dbg_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(shmmr_dbg_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(guided_shmmr_dbg_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(shmmr_sparse_aln_consensus, m)?)?;
    Ok(())
}
