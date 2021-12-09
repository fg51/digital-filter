use super::relative_degree;

//    r"""
//    Transform a lowpass filter prototype to a different frequency.
//
//    Return an analog low-pass filter with cutoff frequency `wo`
//    from an analog low-pass filter prototype with unity cutoff frequency,
//    using zeros, poles, and gain ('zpk') representation.
//
//    Parameters
//    ----------
//    z : array_like
//        Zeros of the analog filter transfer function.
//    p : array_like
//        Poles of the analog filter transfer function.
//    k : float
//        System gain of the analog filter transfer function.
//    wo : float
//        Desired cutoff, as angular frequency (e.g., rad/s).
//        Defaults to no change.
//
//    Returns
//    -------
//    z : ndarray
//        Zeros of the transformed low-pass filter transfer function.
//    p : ndarray
//        Poles of the transformed low-pass filter transfer function.
//    k : float
//        System gain of the transformed low-pass filter.
//
//    See Also
//    --------
//    lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
//    lp2lp
//
//    Notes
//    -----
//    This is derived from the s-plane substitution
//
//    .. math:: s \rightarrow \frac{s}{\omega_0}
//
//    .. versionadded:: 1.1.0
//
//    """

use num_complex::Complex64 as Complex;

//def lp2lp_zpk(z, p, k, wo=1.0):
pub fn lp2lp_zpk(
    z: &[f64],
    p: &[Complex],
    k: f64,
    wo: Option<f64>,
) -> (Vec<f64>, Vec<Complex>, f64) {
    //z = atleast_1d(z)
    //p = atleast_1d(p)
    let wo = match wo {
        //wo = float(wo)  # Avoid int wraparound
        Some(wo) => wo,
        None => 1.,
    };

    let degree = relative_degree(z, p);

    // Scale all points radially from origin to shift cutoff frequency
    let z_lp: Vec<f64> = z.iter().map(|x| wo * x).collect();
    let p_lp: Vec<Complex> = p.iter().map(|x| wo * x).collect();

    // Each shifted pole decreases gain by wo, each shifted zero increases it.
    // Cancel out the net change to keep overall gain the same
    let k_lp = k * (wo.powf(degree as f64));

    return (z_lp, p_lp, k_lp);
}

// //    """
// //    Return relative degree of transfer function from zeros and poles
// //    """
// //def _relative_degree(z, p):
// fn relative_degree(z: &[f64], p: &[Complex]) -> f64 {
//     let degree = p.len() as isize - z.len() as isize;
//     if degree < 0 {
//         //        raise ValueError("Improper transfer function. "
//         //                         "Must have at least as many poles as zeros.")
//         todo!();
//     } else {
//         degree as f64
//     }
// }
