// src/lib.rs
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::PyString;
use pyo3::wrap_pyfunction;
use pyo3::Python;

use pgr_utils::fasta_io::FastaReader;



#[pymodule]
fn pgrlite(_: Python, m: &PyModule) -> PyResult<()> {
    Ok(())
}