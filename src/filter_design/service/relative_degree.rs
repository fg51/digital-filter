use num_complex::Complex64 as Complex;

//    Return relative degree of transfer function from zeros and poles

//def _relative_degree(z, p):
pub fn relative_degree(z: &[f64], p: &[Complex]) -> usize {
    let degree = p.len() as isize - z.len() as isize;
    if degree < 0 {
        //        raise ValueError("Improper transfer function. "
        //                         "Must have at least as many poles as zeros.")
        todo!();
    } else {
        degree as usize
    }
}
