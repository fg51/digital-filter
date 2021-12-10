//    IIR digital and analog filter design given order and critical points.
//
//    Design an Nth-order digital or analog filter and return the filter
//    coefficients.
//
//    Parameters
//    ----------
//    N : int
//        The order of the filter.
//    Wn : array_like
//        A scalar or length-2 sequence giving the critical frequencies.
//
//        For digital filters, `Wn` are in the same units as `fs`. By default,
//        `fs` is 2 half-cycles/sample, so these are normalized from 0 to 1,
//        where 1 is the Nyquist frequency. (`Wn` is thus in
//        half-cycles / sample.)
//
//        For analog filters, `Wn` is an angular frequency (e.g., rad/s).
//    rp : float, optional
//        For Chebyshev and elliptic filters, provides the maximum ripple
//        in the passband. (dB)
//    rs : float, optional
//        For Chebyshev and elliptic filters, provides the minimum attenuation
//        in the stop band. (dB)
//    btype : {'bandpass', 'lowpass', 'highpass', 'bandstop'}, optional
//        The type of filter.  Default is 'bandpass'.
//    analog : bool, optional
//        When True, return an analog filter, otherwise a digital filter is
//        returned.
//    ftype : str, optional
//        The type of IIR filter to design:
//
//            - Butterworth   : 'butter'
//            - Chebyshev I   : 'cheby1'
//            - Chebyshev II  : 'cheby2'
//            - Cauer/elliptic: 'ellip'
//            - Bessel/Thomson: 'bessel'
//
//    output : {'ba', 'zpk', 'sos'}, optional
//        Filter form of the output:
//
//            - second-order sections (recommended): 'sos'
//            - numerator/denominator (default)    : 'ba'
//            - pole-zero                          : 'zpk'
//
//        In general the second-order sections ('sos') form  is
//        recommended because inferring the coefficients for the
//        numerator/denominator form ('ba') suffers from numerical
//        instabilities. For reasons of backward compatibility the default
//        form is the numerator/denominator form ('ba'), where the 'b'
//        and the 'a' in 'ba' refer to the commonly used names of the
//        coefficients used.
//
//        Note: Using the second-order sections form ('sos') is sometimes
//        associated with additional computational costs: for
//        data-intense use cases it is therefore recommended to also
//        investigate the numerator/denominator form ('ba').
//
//    fs : float, optional
//        The sampling frequency of the digital system.
//
//        .. versionadded:: 1.2.0
//
//    Returns
//    -------
//    b, a : ndarray, ndarray
//        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
//        Only returned if ``output='ba'``.
//    z, p, k : ndarray, ndarray, float
//        Zeros, poles, and system gain of the IIR filter transfer
//        function.  Only returned if ``output='zpk'``.
//    sos : ndarray
//        Second-order sections representation of the IIR filter.
//        Only returned if ``output=='sos'``.
//
//    See Also
//    --------
//    butter : Filter design using order and critical points
//    cheby1, cheby2, ellip, bessel
//    buttord : Find order and critical points from passband and stopband spec
//    cheb1ord, cheb2ord, ellipord
//    iirdesign : General filter design using passband and stopband spec
//
//    Notes
//    -----
//    The ``'sos'`` output parameter was added in 0.16.0.
//
//    Examples
//    --------
//    Generate a 17th-order Chebyshev II analog bandpass filter from 50 Hz to
//    200 Hz and plot the frequency response:
//
//    >>> from scipy import signal
//    >>> import matplotlib.pyplot as plt
//
//    >>> b, a = signal.iirfilter(17, [2*np.pi*50, 2*np.pi*200], rs=60,
//    ...                         btype='band', analog=True, ftype='cheby2')
//    >>> w, h = signal.freqs(b, a, 1000)
//    >>> fig = plt.figure()
//    >>> ax = fig.add_subplot(1, 1, 1)
//    >>> ax.semilogx(w / (2*np.pi), 20 * np.log10(np.maximum(abs(h), 1e-5)))
//    >>> ax.set_title('Chebyshev Type II bandpass frequency response')
//    >>> ax.set_xlabel('Frequency [Hz]')
//    >>> ax.set_ylabel('Amplitude [dB]')
//    >>> ax.axis((10, 1000, -100, 10))
//    >>> ax.grid(which='both', axis='both')
//    >>> plt.show()
//
//    Create a digital filter with the same properties, in a system with
//    sampling rate of 2000 Hz, and plot the frequency response. (Second-order
//    sections implementation is required to ensure stability of a filter of
//    this order):
//
//    >>> sos = signal.iirfilter(17, [50, 200], rs=60, btype='band',
//    ...                        analog=False, ftype='cheby2', fs=2000,
//    ...                        output='sos')
//    >>> w, h = signal.sosfreqz(sos, 2000, fs=2000)
//    >>> fig = plt.figure()
//    >>> ax = fig.add_subplot(1, 1, 1)
//    >>> ax.semilogx(w, 20 * np.log10(np.maximum(abs(h), 1e-5)))
//    >>> ax.set_title('Chebyshev Type II bandpass frequency response')
//    >>> ax.set_xlabel('Frequency [Hz]')
//    >>> ax.set_ylabel('Amplitude [dB]')
//    >>> ax.axis((10, 1000, -100, 10))
//    >>> ax.grid(which='both', axis='both')
//    >>> plt.show()

use super::service::lp2bp_zpk;
use super::service::lp2hp_zpk;
use super::service::lp2lp_zpk;
use super::service::zpk2tf;
use super::values::{FilterForm, FilterKind, IIRFilterKind};
use crate::errors::{ErrorKind, Result};

//def iirfilter(N, Wn, rp=None, rs=None, btype='band', analog=False, ftype='butter', output='ba', fs=None):
pub fn iirfilter(
    num: usize,
    wn: &[f64],
    rp: Option<f64>, // ripple
    rs: Option<f64>,
    btype: Option<FilterKind>,
    analog: Option<bool>,
    ftype: Option<IIRFilterKind>,
    output: Option<FilterForm>,
    fs: Option<f64>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let analog = match analog {
        Some(analog) => analog,
        None => false,
    };

    let wn: Vec<f64> = match fs {
        Some(fs) => {
            if analog {
                return Err(ErrorKind::ValueError(
                    "fs cannot be specified for an analog filter".to_string(),
                ));
            } else {
                wn.iter().map(|w| 2. * w / fs).collect()
            }
        }
        None => wn.iter().map(|w| *w).collect(),
    };

    if let Some(rp) = rp {
        if rp < 0. {
            return Err(ErrorKind::ValueError(
                "passband ripple (rp) must be positive".to_string(),
            ));
        }
    }

    if let Some(rs) = rs {
        if rs < 0. {
            return Err(ErrorKind::ValueError(
                "stopband attenuation (rs) must be positive".to_string(),
            ));
        }
    }

    // Get analog lowpass prototype
    let (z, p, k) = match ftype {
        Some(IIRFilterKind::Butterworth) | None => {
            //if typefunc == buttap:
            buttap(num)
        }
        Some(IIRFilterKind::Bessel) => {
            //elif typefunc == besselap:
            //    z, p, k = typefunc(N, norm=bessel_norms[ftype])
            // z,p,k = besselap(num, rs); //  z, p, k = typefunc(N, rs)
            //  bessel_norms = 'phase'(default), delay', 'mag'}
            todo!()
        }
        Some(IIRFilterKind::Chebyshev1) => {
            //elif typefunc == cheb1ap:
            if rp.is_none() {
                return Err(ErrorKind::ValueError(
                    "passband ripple (rp) must be provided to design a Chebyshev I filter."
                        .to_string(),
                ));
            }
            cheb1ap(num, rp.unwrap())
        }
        Some(IIRFilterKind::Chebyshev2) => {
            //elif typefunc == cheb2ap:
            if rs.is_none() {
                return Err(ErrorKind::ValueError(
                    "stopband attenuation (rs) must be provided to design an Chebyshev II filter."
                        .to_string(),
                ));
            }
            //z, p, k = cheb2ap(num, rs.unwrap());

            todo!()
        }
        Some(IIRFilterKind::Elliptic) => {
            //elif typefunc == ellipap:
            if rs.is_none() | rp.is_none() {
                return Err(ErrorKind::ValueError(
                    "Both rp and rs must be provided to design an elliptic filter.".to_string(),
                ));
            }
            //    z, p, k = typefunc(N, rp, rs)
            // z, p, k = ellipap(num, rp.unwrap(), rs.unwrap());
            todo!()
        }
    };

    // Pre-warp frequencies for digital filter design
    let (warped, fs) = if analog == false {
        pre_warp(&wn, fs)?
        //if is_ranged(&wn) == false {
        //    if let Some(fs) = fs {
        //        return Err(ErrorKind::ValueError(format!(
        //            "Digital filter critical frequencies must be 0 < Wn < fs/2 (fs={} -> fs/2={})",
        //            fs,
        //            fs / 2.
        //        )));
        //    }
        //    return Err(ErrorKind::ValueError(
        //        "Digital filter critical frequencies must be 0 < Wn < 1".to_string(),
        //    ));
        //}
        //let fs = 2.0;
        //let warped = wn.iter().map(|x| 2. * fs * ((PI * x / fs).tan())).collect();
        //(warped, Some(fs))
    } else {
        (wn, fs)
    };

    // transform to lowpass, bandpass, highpass, or bandstop
    let (z, p, k) = match btype {
        Some(FilterKind::LowPass) => {
            //if btype in ('lowpass', 'highpass'):
            //    if numpy.size(Wn) != 1:
            //        raise ValueError('Must specify a single critical frequency Wn for lowpass or highpass filter')
            //    if btype == 'lowpass':
            lp2lp_zpk(&z, &p, k, Some(warped[0]))
        }
        Some(FilterKind::HighPass) => {
            //if btype in ('lowpass', 'highpass'):
            //    if numpy.size(Wn) != 1:
            //        raise ValueError('Must specify a single critical frequency Wn for lowpass or highpass filter')
            //    elif btype == 'highpass':
            lp2hp_zpk(&z, &p, k, Some(warped[0]))
        }
        Some(FilterKind::BandPass) | None => {
            //elif btype in ('bandpass', 'bandstop'):
            if warped.len() < 2 {
                return Err(ErrorKind::ValueError(
                    "Wn must specify start and stop frequencies for bandpass or bandstop "
                        .to_string(),
                ));
            }

            let _bw = warped[1] - warped[0];
            let _wo = (warped[0] * warped[1]).sqrt();
            //    if btype == 'bandpass':
            // z, p, klp2bp_zpk(&z, &p, k, Some(wo), Some(bw));
            todo!();
        } //    elif btype == 'bandstop':
        Some(FilterKind::BandStop) => {
            //elif btype in ('bandpass', 'bandstop'):
            //    try:
            let _bw = warped[1] - warped[0];
            let _wo = (warped[0] * warped[1]).sqrt();
            //  z, p, k = lp2bs_zpk(&z,&p, k, Some(wo), Some(bw));
            todo!()
        } //else:
          //    raise NotImplementedError("'%s' not implemented in iirfilter." % btype)
    };

    // Find discrete equivalent if necessary
    let (z, p, k) = if analog == false {
        bilinear_zpk(&z, &p, k, fs.unwrap())
    } else {
        (z, p, k)
    };

    // Transform to proper out type (pole-zero, state-space, numer-denom)
    match output {
        Some(FilterForm::Zpk) => {
            //    return z, p, k
            todo!()
        }
        Some(FilterForm::Sos) => {
            //    return zpk2sos(z, p, k)
            todo!()
        }
        Some(FilterForm::Ba) | None => return Ok(zpk2tf(&z, &p, k)),
    };
}

use num_complex::Complex64 as Complex;
use std::f64::consts::{E, PI};

// Return (z,p,k) for analog prototype of Nth-order Butterworth filter.
// The filter will have an angular (e.g., rad/s) cutoff frequency of 1.
fn buttap(num: usize) -> (Vec<f64>, Vec<Complex>, f64) {
    //if abs(int(N)) != N:
    //    raise ValueError("Filter order must be a nonnegative integer")
    let z = vec![]; //numpy.array([])
    let m = arange(-(num as f64) + 1., num as f64, 2.);
    // Middle value is 0 to ensure an exactly real pole
    //p = -numpy.exp(1j * PI * m / (2 * N));
    let p = m
        .into_iter()
        .map(|x| Complex::new(0., PI * x / (2. * (num as f64))).expf(E) * -1.)
        .collect();
    let k = 1.;
    return (z, p, k);
}

fn arange(start: f64, stop: f64, step: f64) -> Vec<f64> {
    let mut x = start;
    let mut xs = vec![];
    while x < stop {
        xs.push(x);
        x += step;
    }
    return xs;
}

// Return (z,p,k) for Nth-order Chebyshev type I analog lowpass filter.
// The returned filter prototype has `rp` decibels of ripple in the passband.
//
// The filter's angular (e.g. rad/s) cutoff frequency is normalized to 1,
// defined as the point at which the gain first drops below ``-rp``.
fn cheb1ap(num: usize, rp: f64) -> (Vec<f64>, Vec<Complex>, f64) {
    //if abs(int(N)) != N:
    //    raise ValueError("Filter order must be a nonnegative integer")
    if num == 0 {
        //    # Avoid divide-by-zero error
        //    # Even order filters have DC gain of -rp dB
        return (vec![], vec![], 10f64.powf(-rp / 20.));
    }
    let z: Vec<f64> = vec![]; // numpy.array([])

    // Ripple factor (epsilon)
    let eps = (10f64.powf(0.1 * rp) - 1.0).sqrt();
    let mu = 1.0 / (num as f64) * ((1. / eps).asinh());

    // Arrange poles in an ellipse on the left half of the S-plane
    let m = arange(-(num as f64) + 1., num as f64, 2.); //m = numpy.arange(-N+1, N, 2)
    let theta: Vec<f64> = m.iter().map(|x| PI * x / (2. * num as f64)).collect(); // theta = pi * m / (2*N)
    let p: Vec<Complex> = theta
        .iter()
        .map(|t| -1. * Complex::new(mu, *t).sinh())
        .collect(); //p = -sinh(mu + 1j*theta)

    let k = {
        let mut kc = Complex::new(1., 0.); // k = numpy.prod(-p, axis=0).real
        for i in &p {
            kc *= -i;
        }
        let k = kc.re;
        if (num & 1) == 0 {
            // if N % 2 == 0:
            k / ((1. + eps * eps).sqrt())
        } else {
            k
        }
    };
    return (z, p, k);
}

use super::service::relative_degree; //::relative_degree;
                                     //Return a digital IIR filter from an analog one using a bilinear transform.

//Transform a set of poles and zeros from the analog s-plane to the digital
//z-plane using Tustin's method, which substitutes ``(z-1) / (z+1)`` for
//``s``, maintaining the shape of the frequency response.

//Parameters
//----------
//z : array_like
//    Zeros of the analog filter transfer function.
//p : array_like
//    Poles of the analog filter transfer function.
//k : float
//    System gain of the analog filter transfer function.
//fs : float
//    Sample rate, as ordinary frequency (e.g., hertz). No prewarping is
//    done in this function.

//Returns
//-------
//z : ndarray
//    Zeros of the transformed digital filter transfer function.
//p : ndarray
//    Poles of the transformed digital filter transfer function.
//k : float
//    System gain of the transformed digital filter.

//See Also
//--------
//lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk
//bilinear

//Notes
//-----
//.. versionadded:: 1.1.0

//Examples
//--------
//>>> from scipy import signal
//>>> import matplotlib.pyplot as plt

//>>> fs = 100
//>>> bf = 2 * np.pi * np.array([7, 13])
//>>> filts = signal.lti(*signal.butter(4, bf, btype='bandpass', analog=True,
//...                                   output='zpk'))
//>>> filtz = signal.lti(*signal.bilinear_zpk(filts.zeros, filts.poles,
//...                                         filts.gain, fs))
//>>> wz, hz = signal.freqz_zpk(filtz.zeros, filtz.poles, filtz.gain)
//>>> ws, hs = signal.freqs_zpk(filts.zeros, filts.poles, filts.gain,
//...                           worN=fs*wz)
//>>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hz).clip(1e-15)),
//...              label=r'$|H_z(e^{j \omega})|$')
//>>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hs).clip(1e-15)),
//...              label=r'$|H(j \omega)|$')
//>>> plt.legend()
//>>> plt.xlabel('Frequency [Hz]')
//>>> plt.ylabel('Magnitude [dB]')
//>>> plt.grid()
//"""
//def bilinear_zpk(z, p, k, fs):
fn bilinear_zpk(z: &[f64], p: &[Complex], k: f64, fs: f64) -> (Vec<f64>, Vec<Complex>, f64) {
    //z = atleast_1d(z)
    //p = atleast_1d(p)

    let degree = relative_degree(z, p);

    let fs2 = 2.0 * fs;

    // Bilinear transform the poles and zeros
    let mut z_z: Vec<f64> = z.iter().map(|z| (fs2 + z) / (fs2 - z)).collect();
    let p_z: Vec<Complex> = p.iter().map(|p| (fs2 + p) / (fs2 - p)).collect();

    // Any zeros that were at infinity get moved to the Nyquist frequency
    //z_z = append(z_z, -ones(degree))
    z_z.extend(vec![-1.; degree]);

    // Compensate for gain change
    let prod_fsz = {
        let mut fsz = 1.;
        for i in z.iter() {
            fsz *= fs2 - i;
        }
        fsz
    };
    let prod_fsp = {
        let mut fsp = Complex::new(1., 0.);
        for i in p.iter() {
            fsp *= Complex::new(fs2, 0.) - i;
        }
        fsp
    };
    let k_z = k * ((Complex::new(prod_fsz, 0.) / prod_fsp).re);
    //let k_z = k * real(prod(fs2 - z) / prod(fs2 - p))

    return (z_z, p_z, k_z);
}

fn is_ranged(wn: &[f64]) -> bool {
    let (min, max) = (0., 1.);
    for &i in wn {
        if i <= min {
            return false;
        }
        if i >= max {
            return false;
        }
    }
    return true;
}

// Pre-warp frequencies for digital filter design
fn pre_warp(wn: &[f64], fs: Option<f64>) -> Result<(Vec<f64>, Option<f64>)> {
    if is_ranged(&wn) == false {
        if let Some(fs) = fs {
            return Err(ErrorKind::ValueError(format!(
                "Digital filter critical frequencies must be 0 < Wn < fs/2 (fs={} -> fs/2={})",
                fs,
                fs / 2.
            )));
        }
        return Err(ErrorKind::ValueError(
            "Digital filter critical frequencies must be 0 < Wn < 1".to_string(),
        ));
    }
    let fs = 2.0;
    let warped = wn.iter().map(|x| 2. * fs * ((PI * x / fs).tan())).collect();
    Ok((warped, Some(fs)))
}
