//    r"""Return a Kaiser window.
//
//    The Kaiser window is a taper formed by using a Bessel function.
//
//    Parameters
//    ----------
//    M : int
//        Number of points in the output window. If zero or less, an empty
//        array is returned.
//    beta : float
//        Shape parameter, determines trade-off between main-lobe width and
//        side lobe level. As beta gets large, the window narrows.
//    sym : bool, optional
//        When True (default), generates a symmetric window, for use in filter
//        design.
//        When False, generates a periodic window, for use in spectral analysis.
//
//    Returns
//    -------
//    w : ndarray
//        The window, with the maximum value normalized to 1 (though the value 1
//        does not appear if `M` is even and `sym` is True).
//
//    Notes
//    -----
//    The Kaiser window is defined as
//
//    .. math::  w(n) = I_0\left( \beta \sqrt{1-\frac{4n^2}{(M-1)^2}}
//               \right)/I_0(\beta)
//
//    with
//
//    .. math:: \quad -\frac{M-1}{2} \leq n \leq \frac{M-1}{2},
//
//    where :math:`I_0` is the modified zeroth-order Bessel function.
//
//    The Kaiser was named for Jim Kaiser, who discovered a simple approximation
//    to the DPSS window based on Bessel functions.
//    The Kaiser window is a very good approximation to the Digital Prolate
//    Spheroidal Sequence, or Slepian window, which is the transform which
//    maximizes the energy in the main lobe of the window relative to total
//    energy.
//
//    The Kaiser can approximate other windows by varying the beta parameter.
//    (Some literature uses alpha = beta/pi.) [4]_
//
//    ====  =======================
//    beta  Window shape
//    ====  =======================
//    0     Rectangular
//    5     Similar to a Hamming
//    6     Similar to a Hann
//    8.6   Similar to a Blackman
//    ====  =======================
//
//    A beta value of 14 is probably a good starting point. Note that as beta
//    gets large, the window narrows, and so the number of samples needs to be
//    large enough to sample the increasingly narrow spike, otherwise NaNs will
//    be returned.
//
//    Most references to the Kaiser window come from the signal processing
//    literature, where it is used as one of many windowing functions for
//    smoothing values.  It is also known as an apodization (which means
//    "removing the foot", i.e. smoothing discontinuities at the beginning
//    and end of the sampled signal) or tapering function.
//
//    References
//    ----------
//    .. [1] J. F. Kaiser, "Digital Filters" - Ch 7 in "Systems analysis by
//           digital computer", Editors: F.F. Kuo and J.F. Kaiser, p 218-285.
//           John Wiley and Sons, New York, (1966).
//    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The
//           University of Alberta Press, 1975, pp. 177-178.
//    .. [3] Wikipedia, "Window function",
//           https://en.wikipedia.org/wiki/Window_function
//    .. [4] F. J. Harris, "On the use of windows for harmonic analysis with the
//           discrete Fourier transform," Proceedings of the IEEE, vol. 66,
//           no. 1, pp. 51-83, Jan. 1978. :doi:`10.1109/PROC.1978.10837`.
//
//    Examples
//    --------
//    Plot the window and its frequency response:
//
//    >>> from scipy import signal
//    >>> from scipy.fft import fft, fftshift
//    >>> import matplotlib.pyplot as plt
//
//    >>> window = signal.windows.kaiser(51, beta=14)
//    >>> plt.plot(window)
//    >>> plt.title(r"Kaiser window ($\beta$=14)")
//    >>> plt.ylabel("Amplitude")
//    >>> plt.xlabel("Sample")
//
//    >>> plt.figure()
//    >>> A = fft(window, 2048) / (len(window)/2.0)
//    >>> freq = np.linspace(-0.5, 0.5, len(A))
//    >>> response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
//    >>> plt.plot(freq, response)
//    >>> plt.axis([-0.5, 0.5, -120, 0])
//    >>> plt.title(r"Frequency response of the Kaiser window ($\beta$=14)")
//    >>> plt.ylabel("Normalized magnitude [dB]")
//    >>> plt.xlabel("Normalized frequency [cycles per sample]")
//
//    """

//def kaiser(M, beta, sym=True):
pub fn kaiser(m: usize, beta: f64, sym: Option<bool>) -> Vec<f64> {
    //if _len_guards(M):
    //    return np.ones(M)
    //M, needs_trunc = _extend(M, sym)
    let (m, needs_trunc) = match sym {
        Some(true) => (m, false),
        Some(false) | None => (m + 1, true),
    };

    // np.arange(0, M)
    let alpha = (m as f64 - 1.) / 2.0;
    //w = (special.i0(beta * np.sqrt(1 - ((n - alpha) / alpha) ** 2.0)) /
    //     special.i0(beta))
    let denom = i0(beta);
    let mut w = Vec::with_capacity(m);
    w.resize(m, 0.);
    for i in 0..m {
        let ia = (i as f64 - alpha) / alpha;
        let ia2 = (1. - ia * ia).sqrt();
        w[i] = i0(beta * ia2) / i0(beta);
    }

    return truncate(w, needs_trunc);
}

// Truncate window by 1 sample if needed for DFT-even symmetry
//def _truncate(w, needed):
fn truncate(w: Vec<f64>, needed: bool) -> Vec<f64> {
    if needed {
        return w.into_iter().rev().collect(); //return w[:-1];
    } else {
        return w;
    }
}

//double i0(x) double x; {
fn i0(x: f64) -> f64 {
    //    double y;

    let x = if x < 0. { -x } else { x };
    if x <= 8.0 {
        let y = (x / 2.0) - 2.0;
        return x.exp() * chbevl(y, &A, 30);
    }
    return x.exp() * chbevl(32.0 / x - 2.0, &B, 25) / (x.sqrt());
}

//double chbevl(double x, double array[], int n)
fn chbevl(x: f64, array: &[f64], n: usize) -> f64 {
    //double b0, b1, b2, *p;
    //int i;

    //p = array;
    let mut b0 = array[0]; // b0 = *p++;
                           //let array = array[1..];
    let mut b1 = 0.0;
    let mut i = n - 1;

    let mut b2 = 0.;
    for p in &array[1..] {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + p; //b0 = x * b1 - b2 + *p++;
                              //while (--i);
        i -= 1;
        if i == 0 {
            break;
        }
    }

    return 0.5 * (b0 - b2);
}

const A: [f64; 30] = [
    -4.41534164647933937950E-18,
    3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
    1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
    7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
    2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
    9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
    2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
    6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
    1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
    1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
    1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
    1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
    1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
    4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
    1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
    6.76795274409476084995E-1,
];

const B: [f64; 25] = [
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
    4.46562142029675999901E-17,
    3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
    1.77256013305652638360E-15,
    3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
    1.54008621752140982691E-14,
    3.85277838274214270114E-13,
    7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
    1.18891471078464383424E-11,
    4.94060238822496958910E-10,
    3.39623202570838634515E-9,
    2.26666899049817806459E-8,
    2.04891858946906374183E-7,
    2.89137052083475648297E-6,
    6.88975834691682398426E-5,
    3.36911647825569408990E-3,
    8.04490411014108831608E-1,
];
