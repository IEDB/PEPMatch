use pyo3::prelude::*;

mod preprocess;
#[path = "match.rs"]
mod matching;

#[pyfunction]
fn rs_version() -> &'static str {
    "0.1.0"
}

#[pyfunction]
fn rs_preprocess(fasta_path: &str, k: usize, output_path: &str) {
    preprocess::run(fasta_path, k, output_path);
}

#[pyfunction]
fn rs_match(pepidx_path: &str, peptides: Vec<(String, String)>, k: usize, max_mismatches: usize) -> Vec<Vec<String>> {
    matching::run(pepidx_path, peptides, k, max_mismatches)
}

#[pymodule]
fn _rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rs_version, m)?)?;
    m.add_function(wrap_pyfunction!(rs_preprocess, m)?)?;
    m.add_function(wrap_pyfunction!(rs_match, m)?)?;
    Ok(())
}
