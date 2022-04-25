// src/lib.rs
use pgr_db::agc_io;
use pgr_db::seq_db::{query_fragment, read_shmr_map_file, FragmentSignature, ShmmrToFrags};
use pgr_utils::fasta_io;
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::PyString;
use pyo3::wrap_pyfunction;
use pyo3::Python;
use rustc_hash::FxHashMap;

#[pyclass]
#[derive(Clone)]
struct ShmmrFragMap {
    pub shmmr_to_frags: ShmmrToFrags,
}

#[pymethods]
impl ShmmrFragMap {
    #[new]
    pub fn new() -> Self {
        let shmmr_to_frags = ShmmrToFrags::default();
        ShmmrFragMap { shmmr_to_frags }
    }

    pub fn load_from_file(&mut self, filename: String) -> () {
        let new_map = read_shmr_map_file(filename).unwrap();
        self.shmmr_to_frags = new_map;
    }

    pub fn query_fragment(
        &self,
        seq: Vec<u8>,
    ) -> PyResult<Vec<((u64, u64), (u32, u32, u8), Vec<FragmentSignature>)>> {
        let res: Vec<((u64, u64), (u32, u32, u8), Vec<FragmentSignature>)> =
            query_fragment(&self.shmmr_to_frags, &seq);
        Ok(res)
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

#[pymodule]
fn pgrlite(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<ShmmrFragMap>()?;
    m.add_class::<AGCFile>()?;
    Ok(())
}
