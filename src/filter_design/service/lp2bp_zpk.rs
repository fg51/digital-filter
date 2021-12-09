//r"""
//Transform a lowpass filter prototype to a bandpass filter.

//Return an analog band-pass filter with center frequency `wo` and
//bandwidth `bw` from an analog low-pass filter prototype with unity
//cutoff frequency, using zeros, poles, and gain ('zpk') representation.

//Parameters
//----------
//z : array_like
//    Zeros of the analog filter transfer function.
//p : array_like
//    Poles of the analog filter transfer function.
//k : float
//    System gain of the analog filter transfer function.
//wo : float
//    Desired passband center, as angular frequency (e.g., rad/s).
//    Defaults to no change.
//bw : float
//    Desired passband width, as angular frequency (e.g., rad/s).
//    Defaults to 1.

//Returns
//-------
//z : ndarray
//    Zeros of the transformed band-pass filter transfer function.
//p : ndarray
//    Poles of the transformed band-pass filter transfer function.
//k : float
//    System gain of the transformed band-pass filter.

//See Also
//--------
//lp2lp_zpk, lp2hp_zpk, lp2bs_zpk, bilinear
//lp2bp

//Notes
//-----
//This is derived from the s-plane substitution

//.. math:: s \rightarrow \frac{s^2 + {\omega_0}^2}{s \cdot \mathrm{BW}}

//This is the "wideband" transformation, producing a passband with
//geometric (log frequency) symmetry about `wo`.

//.. versionadded:: 1.1.0

//"""
use super::relative_degree;
use num_complex::Complex64 as Complex;

//def lp2bp_zpk(z, p, k, wo=1.0, bw=1.0):
pub fn lp2bp_zpk(
    z: &[f64],
    p: &[Complex],
    k: f64,
    wo: Option<f64>,
    bw: Option<f64>,
) -> (Vec<Complex>, Vec<Complex>, f64) {
    //    z = atleast_1d(z)
    //    p = atleast_1d(p)
    let wo = match wo {
        Some(wo) => wo,
        None => 1.,
    };
    let bw = match bw {
        Some(bw) => bw,
        None => 1.,
    };

    let degree = relative_degree(z, p);

    // Scale poles and zeros to desired bandwidth
    let z_lp: Vec<Complex> = z.iter().map(|z| Complex::new(z * bw / 2., 0.)).collect();
    let p_lp: Vec<Complex> = p.iter().map(|p| p * bw / 2.).collect();

    // Duplicate poles and zeros and shift from baseband to +wo and -wo
    //    z_bp = concatenate((z_lp + sqrt(z_lp**2 - wo**2), z_lp - sqrt(z_lp**2 - wo**2)))
    let mut z_bp: Vec<Complex> = z_lp
        .iter()
        .map(|z| z + ((z * z - wo * wo).sqrt()))
        .collect();
    z_bp.extend(
        z_lp.iter()
            .map(|z| z - ((z * z - wo * wo).sqrt()))
            .collect::<Vec<Complex>>(),
    );
    //    p_bp = concatenate((p_lp + sqrt(p_lp**2 - wo**2),
    //                        p_lp - sqrt(p_lp**2 - wo**2)))
    let mut p_bp: Vec<Complex> = p_lp
        .iter()
        .map(|p| p + ((p * p - wo * wo).sqrt()))
        .collect();
    p_bp.extend(
        p_lp.iter()
            .map(|p| p - ((p * p - wo * wo).sqrt()))
            .collect::<Vec<Complex>>(),
    );

    // Move degree zeros to origin, leaving degree zeros at infinity for BPF
    z_bp.extend(vec![Complex::new(0., 0.); degree]);

    // Cancel out gain change from frequency scaling
    let k_bp = k * (bw.powf(degree as f64));

    return (z_bp, p_bp, k_bp);
}
