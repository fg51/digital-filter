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

use super::service::lp2lp_zpk;
use super::service::zpk2tf;
use super::values::{FilterForm, FilterKind, IIRFilterKind};

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
) -> (Vec<f64>, Vec<f64>) {
    let analog = match analog {
        Some(analog) => analog,
        None => false,
    };

    //ftype, btype, output = [x.lower() for x in (ftype, btype, output)]
    //Wn = asarray(Wn)
    let wn: Vec<f64> = match fs {
        Some(fs) => {
            if analog {
                // raise ValueError("fs cannot be specified for an analog filter")
                todo!();
            } else {
                //    Wn = 2*Wn/fs
                wn.iter().map(|w| 2. * w / fs).collect()
            }
        }
        None => wn.iter().map(|w| *w).collect(),
    };

    //typefunc = filter_dict[ftype][0]

    //if rp is not None and rp < 0:
    //    raise ValueError("passband ripple (rp) must be positive")

    //if rs is not None and rs < 0:
    //    raise ValueError("stopband attenuation (rs) must be positive")

    // Get analog lowpass prototype
    let (z, p, k) = match ftype {
        Some(IIRFilterKind::Butterworth) | None => {
            //if typefunc == buttap:
            buttap(num)
            //    z, p, k = typefunc(N)
            //todo!()
        }
        Some(IIRFilterKind::Bessel) => {
            //elif typefunc == besselap:
            //    z, p, k = typefunc(N, norm=bessel_norms[ftype])
            todo!()
        }
        Some(IIRFilterKind::Chebyshev1) => {
            //elif typefunc == cheb1ap:
            //    if rp is None:
            //        raise ValueError("passband ripple (rp) must be provided to "
            //                         "design a Chebyshev I filter.")
            cheb1ap(num, rp.unwrap())
        }
        Some(IIRFilterKind::Chebyshev2) => {
            //elif typefunc == cheb2ap:
            //    if rs is None:
            //        raise ValueError("stopband attenuation (rs) must be provided to "
            //                         "design an Chebyshev II filter.")
            //    z, p, k = typefunc(N, rs)
            todo!()
        }
        Some(IIRFilterKind::Elliptic) => {
            //elif typefunc == ellipap:
            //    if rs is None or rp is None:
            //        raise ValueError("Both rp and rs must be provided to design an "
            //                         "elliptic filter.")
            //    z, p, k = typefunc(N, rp, rs)
            todo!()
        } //else:
          //    raise NotImplementedError("'%s' not implemented in iirfilter." % ftype)
    };

    // Pre-warp frequencies for digital filter design
    let warped = if analog {
        //    if numpy.any(Wn <= 0) or numpy.any(Wn >= 1):
        //        if fs is not None:
        //            raise ValueError("Digital filter critical frequencies "
        //                             "must be 0 < Wn < fs/2 (fs={} -> fs/2={})".format(fs, fs/2))
        //        raise ValueError("Digital filter critical frequencies "
        //                         "must be 0 < Wn < 1")
        //    fs = 2.0
        //    warped = 2 * fs * tan(pi * Wn / fs)
        todo!();
    } else {
        wn
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
            //        z, p, k = lp2hp_zpk(z, p, k, wo=warped)
            todo!()
        }
        //elif btype in ('bandpass', 'bandstop'):
        //    try:
        //        bw = warped[1] - warped[0]
        //        wo = sqrt(warped[0] * warped[1])
        //    except IndexError as e:
        //        raise ValueError('Wn must specify start and stop frequencies for bandpass or bandstop '
        //                         'filter') from e
        Some(FilterKind::BandPass) | None => {
            //elif btype in ('bandpass', 'bandstop'):
            //    try:
            //        bw = warped[1] - warped[0]
            //        wo = sqrt(warped[0] * warped[1])
            //    if btype == 'bandpass':
            //        z, p, k = lp2bp_zpk(z, p, k, wo=wo, bw=bw)
            todo!()
        } //    elif btype == 'bandstop':
        Some(FilterKind::BandStop) => {
            //elif btype in ('bandpass', 'bandstop'):
            //    try:
            //        bw = warped[1] - warped[0]
            //        wo = sqrt(warped[0] * warped[1])
            //        z, p, k = lp2bs_zpk(z, p, k, wo=wo, bw=bw)
            todo!()
        } //else:
          //    raise NotImplementedError("'%s' not implemented in iirfilter." % btype)
    };

    // Find discrete equivalent if necessary
    let (z, p, k) = if analog {
        // bilinear_zpk(z, p, k, fs=fs)
        todo!();
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
        Some(FilterForm::Ba) | None => return zpk2tf(&z, &p, k),
    };
}

use num_complex::Complex64 as Complex;
use std::f64::consts::{E, PI};

// """Return (z,p,k) for analog prototype of Nth-order Butterworth filter.
//
// The filter will have an angular (e.g., rad/s) cutoff frequency of 1.
//
// See Also
// --------
// butter : Filter design function using this prototype
//def buttap(N):
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

//    """
//    Return (z,p,k) for Nth-order Chebyshev type I analog lowpass filter.
//
//    The returned filter prototype has `rp` decibels of ripple in the passband.
//
//    The filter's angular (e.g. rad/s) cutoff frequency is normalized to 1,
//    defined as the point at which the gain first drops below ``-rp``.
//
//    See Also
//    --------
//    cheby1 : Filter design function using this prototype
//
//    """
//def cheb1ap(N, rp):
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
