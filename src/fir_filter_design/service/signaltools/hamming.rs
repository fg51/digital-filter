use std::f64::consts::PI;

use super::super::linspace;

//    r"""Return a Hamming window.
//
//    The Hamming window is a taper formed by using a raised cosine with
//    non-zero endpoints, optimized to minimize the nearest side lobe.
//
//    Parameters
//    ----------
//    M : int
//        Number of points in the output window. If zero or less, an empty
//        array is returned.
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
//    The Hamming window is defined as
//
//    .. math::  w(n) = 0.54 - 0.46 \cos\left(\frac{2\pi{n}}{M-1}\right)
//               \qquad 0 \leq n \leq M-1
//
//    The Hamming was named for R. W. Hamming, an associate of J. W. Tukey and
//    is described in Blackman and Tukey. It was recommended for smoothing the
//    truncated autocovariance function in the time domain.
//    Most references to the Hamming window come from the signal processing
//    literature, where it is used as one of many windowing functions for
//    smoothing values.  It is also known as an apodization (which means
//    "removing the foot", i.e. smoothing discontinuities at the beginning
//    and end of the sampled signal) or tapering function.
//
//    References
//    ----------
//    .. [1] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power
//           spectra, Dover Publications, New York.
//    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The
//           University of Alberta Press, 1975, pp. 109-110.
//    .. [3] Wikipedia, "Window function",
//           https://en.wikipedia.org/wiki/Window_function
//    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
//           "Numerical Recipes", Cambridge University Press, 1986, page 425.
//
//    Examples
//    --------
//    Plot the window and its frequency response:
//
//    >>> from scipy import signal
//    >>> from scipy.fft import fft, fftshift
//    >>> import matplotlib.pyplot as plt
//
//    >>> window = signal.windows.hamming(51)
//    >>> plt.plot(window)
//    >>> plt.title("Hamming window")
//    >>> plt.ylabel("Amplitude")
//    >>> plt.xlabel("Sample")
//
//    >>> plt.figure()
//    >>> A = fft(window, 2048) / (len(window)/2.0)
//    >>> freq = np.linspace(-0.5, 0.5, len(A))
//    >>> response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
//    >>> plt.plot(freq, response)
//    >>> plt.axis([-0.5, 0.5, -120, 0])
//    >>> plt.title("Frequency response of the Hamming window")
//    >>> plt.ylabel("Normalized magnitude [dB]")
//    >>> plt.xlabel("Normalized frequency [cycles per sample]")

//def hamming(M, sym=True):
pub fn hamming(m: usize, sym: Option<bool>) -> Vec<f64> {
    return general_hamming(m, 0.54, sym);
}

//    """Return a generalized Hamming window.
//
//    The generalized Hamming window is constructed by multiplying a rectangular
//    window by one period of a cosine function [1]_.
//
//    Parameters
//    ----------
//    M : int
//        Number of points in the output window. If zero or less, an empty
//        array is returned.
//    alpha : float
//        The window coefficient, :math:`\alpha`
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
//    The generalized Hamming window is defined as
//
//    .. math:: w(n) = \alpha - \left(1 - \alpha\right) \cos\left(\frac{2\pi{n}}{M-1}\right)
//              \qquad 0 \leq n \leq M-1
//
//    Both the common Hamming window and Hann window are special cases of the
//    generalized Hamming window with :math:`\alpha` = 0.54 and :math:`\alpha` =
//    0.5, respectively [2]_.
//
//    See Also
//    --------
//    hamming, hann
//
//    Examples
//    --------
//    The Sentinel-1A/B Instrument Processing Facility uses generalized Hamming
//    windows in the processing of spaceborne Synthetic Aperture Radar (SAR)
//    data [3]_. The facility uses various values for the :math:`\alpha`
//    parameter based on operating mode of the SAR instrument. Some common
//    :math:`\alpha` values include 0.75, 0.7 and 0.52 [4]_. As an example, we
//    plot these different windows.
//
//    >>> from scipy.signal.windows import general_hamming
//    >>> from scipy.fft import fft, fftshift
//    >>> import matplotlib.pyplot as plt
//
//    >>> fig1, spatial_plot = plt.subplots()
//    >>> spatial_plot.set_title("Generalized Hamming Windows")
//    >>> spatial_plot.set_ylabel("Amplitude")
//    >>> spatial_plot.set_xlabel("Sample")
//
//    >>> fig2, freq_plot = plt.subplots()
//    >>> freq_plot.set_title("Frequency Responses")
//    >>> freq_plot.set_ylabel("Normalized magnitude [dB]")
//    >>> freq_plot.set_xlabel("Normalized frequency [cycles per sample]")
//
//    >>> for alpha in [0.75, 0.7, 0.52]:
//    ...     window = general_hamming(41, alpha)
//    ...     spatial_plot.plot(window, label="{:.2f}".format(alpha))
//    ...     A = fft(window, 2048) / (len(window)/2.0)
//    ...     freq = np.linspace(-0.5, 0.5, len(A))
//    ...     response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
//    ...     freq_plot.plot(freq, response, label="{:.2f}".format(alpha))
//    >>> freq_plot.legend(loc="upper right")
//    >>> spatial_plot.legend(loc="upper right")
//
//    References
//    ----------
//    .. [1] DSPRelated, "Generalized Hamming Window Family",
//           https://www.dsprelated.com/freebooks/sasp/Generalized_Hamming_Window_Family.html
//    .. [2] Wikipedia, "Window function",
//           https://en.wikipedia.org/wiki/Window_function
//    .. [3] Riccardo Piantanida ESA, "Sentinel-1 Level 1 Detailed Algorithm
//           Definition",
//           https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Level-1-Detailed-Algorithm-Definition
//    .. [4] Matthieu Bourbigot ESA, "Sentinel-1 Product Definition",
//           https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Definition

//def general_hamming(M, alpha, sym=True):
fn general_hamming(m: usize, alpha: f64, sym: Option<bool>) -> Vec<f64> {
    return general_cosine(m, &[alpha, 1. - alpha], sym);
}

//    r"""
//    Generic weighted sum of cosine terms window
//
//    Parameters
//    ----------
//    M : int
//        Number of points in the output window
//    a : array_like
//        Sequence of weighting coefficients. This uses the convention of being
//        centered on the origin, so these will typically all be positive
//        numbers, not alternating sign.
//    sym : bool, optional
//        When True (default), generates a symmetric window, for use in filter
//        design.
//        When False, generates a periodic window, for use in spectral analysis.
//
//    References
//    ----------
//    .. [1] A. Nuttall, "Some windows with very good sidelobe behavior," IEEE
//           Transactions on Acoustics, Speech, and Signal Processing, vol. 29,
//           no. 1, pp. 84-91, Feb 1981. :doi:`10.1109/TASSP.1981.1163506`.
//    .. [2] Heinzel G. et al., "Spectrum and spectral density estimation by the
//           Discrete Fourier transform (DFT), including a comprehensive list of
//           window functions and some new flat-top windows", February 15, 2002
//           https://holometer.fnal.gov/GH_FFT.pdf
//
//    Examples
//    --------
//    Heinzel describes a flat-top window named "HFT90D" with formula: [2]_
//
//    .. math::  w_j = 1 - 1.942604 \cos(z) + 1.340318 \cos(2z)
//               - 0.440811 \cos(3z) + 0.043097 \cos(4z)
//
//    where
//
//    .. math::  z = \frac{2 \pi j}{N}, j = 0...N - 1
//
//    Since this uses the convention of starting at the origin, to reproduce the
//    window, we need to convert every other coefficient to a positive number:
//
//    >>> HFT90D = [1, 1.942604, 1.340318, 0.440811, 0.043097]
//
//    The paper states that the highest sidelobe is at -90.2 dB.  Reproduce
//    Figure 42 by plotting the window and its frequency response, and confirm
//    the sidelobe level in red:
//
//    >>> from scipy.signal.windows import general_cosine
//    >>> from scipy.fft import fft, fftshift
//    >>> import matplotlib.pyplot as plt
//
//    >>> window = general_cosine(1000, HFT90D, sym=False)
//    >>> plt.plot(window)
//    >>> plt.title("HFT90D window")
//    >>> plt.ylabel("Amplitude")
//    >>> plt.xlabel("Sample")
//
//    >>> plt.figure()
//    >>> A = fft(window, 10000) / (len(window)/2.0)
//    >>> freq = np.linspace(-0.5, 0.5, len(A))
//    >>> response = np.abs(fftshift(A / abs(A).max()))
//    >>> response = 20 * np.log10(np.maximum(response, 1e-10))
//    >>> plt.plot(freq, response)
//    >>> plt.axis([-50/1000, 50/1000, -140, 0])
//    >>> plt.title("Frequency response of the HFT90D window")
//    >>> plt.ylabel("Normalized magnitude [dB]")
//    >>> plt.xlabel("Normalized frequency [cycles per sample]")
//    >>> plt.axhline(-90.2, color='red')
//    >>> plt.show()
//    """
//def general_cosine(M, a, sym=True):
fn general_cosine(m: usize, a: &[f64], sym: Option<bool>) -> Vec<f64> {
    //    if _len_guards(M):
    //        return np.ones(M)
    // m, needs_trunc = _extend(m, sym);
    let (m, needs_trunc) = match sym {
        Some(true) => (m, false),
        Some(false) | None => (m + 1, true),
    };

    //    fac = np.linspace(-np.pi, np.pi, M)
    let fac = linspace(-PI, PI, m, None);
    let mut w = vec![0.; m]; //w = np.zeros(M)
                             //    for k in range(len(a)):
                             //        w += a[k] * np.cos(k * fac)
                             //
    for k in 0..a.len() {
        for i in 0..w.len() {
            w[i] += a[k] * ((k as f64 * fac[i]).cos()); //w += a[k] * np.cos(k * fac)
        }
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
