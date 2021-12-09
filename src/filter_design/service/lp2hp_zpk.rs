//r"""
//Transform a lowpass filter prototype to a highpass filter.

//Return an analog high-pass filter with cutoff frequency `wo`
//from an analog low-pass filter prototype with unity cutoff frequency,
//using zeros, poles, and gain ('zpk') representation.

//Parameters
//----------
//z : array_like
//    Zeros of the analog filter transfer function.
//p : array_like
//    Poles of the analog filter transfer function.
//k : float
//    System gain of the analog filter transfer function.
//wo : float
//    Desired cutoff, as angular frequency (e.g., rad/s).
//    Defaults to no change.

//Returns
//-------
//z : ndarray
//    Zeros of the transformed high-pass filter transfer function.
//p : ndarray
//    Poles of the transformed high-pass filter transfer function.
//k : float
//    System gain of the transformed high-pass filter.

//See Also
//--------
//lp2lp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
//lp2hp

//Notes
//-----
//This is derived from the s-plane substitution

//.. math:: s \rightarrow \frac{\omega_0}{s}

//This maintains symmetry of the lowpass and highpass responses on a
//logarithmic scale.

//.. versionadded:: 1.1.0

use super::relative_degree;
use num_complex::Complex64 as Complex;

//def lp2hp_zpk(z, p, k, wo=1.0):
pub fn lp2hp_zpk(
    z: &[f64],
    p: &[Complex],
    k: f64,
    wo: Option<f64>,
) -> (Vec<f64>, Vec<Complex>, f64) {
    let wo = match wo {
        //wo = float(wo)  # Avoid int wraparound
        Some(wo) => wo,
        None => 1.,
    };
    let degree = relative_degree(z, p);

    // Invert positions radially about unit circle to convert LPF to HPF
    // Scale all points radially from origin to shift cutoff frequency
    let mut z_hp: Vec<f64> = z.iter().map(|x| wo / x).collect();
    let p_hp: Vec<Complex> = p.iter().map(|x| wo / x).collect();

    // If lowpass had zeros at infinity, inverting moves them to origin.
    z_hp.extend(vec![0.; degree]);

    // Cancel out gain change caused by inversion
    let prod_z = {
        let mut prod = 1.;
        for i in z {
            prod *= -i
        }
        prod
    };
    let prod_p = {
        let mut prod = Complex::new(1., 0.);
        for i in p {
            prod *= -i;
        }
        prod
    };
    //k_hp = k * real(prod(-z) / prod(-p))
    let k_hp = k * ((Complex::new(prod_z, 0.) / (prod_p)).re);

    return (z_hp, p_hp, k_hp);
}
