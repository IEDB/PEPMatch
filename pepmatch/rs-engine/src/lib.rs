use pyo3::prelude::*;

mod preprocess;

#[pyfunction]
fn rs_version() -> &'static str {
    "0.1.0"
}

#[pyfunction]
fn rs_preprocess(fasta_path: &str, k: usize, output_path: &str) {
    preprocess::run(fasta_path, k, output_path);
}

#[pymodule]
fn _rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rs_version, m)?)?;
    m.add_function(wrap_pyfunction!(rs_preprocess, m)?)?;
    Ok(())
}
