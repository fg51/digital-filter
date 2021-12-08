use super::hamming::hamming;
use super::kaiser::kaiser;

//    """
//    Return a window of a given length and type.
//
//    Parameters
//    ----------
//    window : string, float, or tuple
//        The type of window to create. See below for more details.
//    Nx : int
//        The number of samples in the window.
//    fftbins : bool, optional
//        If True (default), create a "periodic" window, ready to use with
//        `ifftshift` and be multiplied by the result of an FFT (see also
//        :func:`~scipy.fft.fftfreq`).
//        If False, create a "symmetric" window, for use in filter design.
//
//    Returns
//    -------
//    get_window : ndarray
//        Returns a window of length `Nx` and type `window`
//
//    Notes
//    -----
//    Window types:
//
//    - `~scipy.signal.windows.boxcar`
//    - `~scipy.signal.windows.triang`
//    - `~scipy.signal.windows.blackman`
//    - `~scipy.signal.windows.hamming`
//    - `~scipy.signal.windows.hann`
//    - `~scipy.signal.windows.bartlett`
//    - `~scipy.signal.windows.flattop`
//    - `~scipy.signal.windows.parzen`
//    - `~scipy.signal.windows.bohman`
//    - `~scipy.signal.windows.blackmanharris`
//    - `~scipy.signal.windows.nuttall`
//    - `~scipy.signal.windows.barthann`
//    - `~scipy.signal.windows.kaiser` (needs beta)
//    - `~scipy.signal.windows.gaussian` (needs standard deviation)
//    - `~scipy.signal.windows.general_gaussian` (needs power, width)
//    - `~scipy.signal.windows.dpss` (needs normalized half-bandwidth)
//    - `~scipy.signal.windows.chebwin` (needs attenuation)
//    - `~scipy.signal.windows.exponential` (needs center, decay scale)
//    - `~scipy.signal.windows.tukey` (needs taper fraction)
//    - `~scipy.signal.windows.taylor` (needs number of constant sidelobes,
//      sidelobe level)
//
//    If the window requires no parameters, then `window` can be a string.
//
//    If the window requires parameters, then `window` must be a tuple
//    with the first argument the string name of the window, and the next
//    arguments the needed parameters.
//
//    If `window` is a floating point number, it is interpreted as the beta
//    parameter of the `~scipy.signal.windows.kaiser` window.
//
//    Each of the window types listed above is also the name of
//    a function that can be called directly to create a window of
//    that type.
//
//    Examples
//    --------
//    >>> from scipy import signal
//    >>> signal.get_window('triang', 7)
//    array([ 0.125,  0.375,  0.625,  0.875,  0.875,  0.625,  0.375])
//    >>> signal.get_window(('kaiser', 4.0), 9)
//    array([ 0.08848053,  0.29425961,  0.56437221,  0.82160913,  0.97885093,
//            0.97885093,  0.82160913,  0.56437221,  0.29425961])
//    >>> signal.get_window(('exponential', None, 1.), 9)
//    array([ 0.011109  ,  0.03019738,  0.082085  ,  0.22313016,  0.60653066,
//            0.60653066,  0.22313016,  0.082085  ,  0.03019738])
//    >>> signal.get_window(4.0, 9)
//    array([ 0.08848053,  0.29425961,  0.56437221,  0.82160913,  0.97885093,
//            0.97885093,  0.82160913,  0.56437221,  0.29425961])
//
//    """
use super::super::super::firwin::WindowKind;

//def get_window(window, Nx, fftbins=True):
pub fn get_window(window: WindowKind, nx: usize, fftbins: Option<bool>) -> Vec<f64> {
    let fftbins = match fftbins {
        Some(fftbins) => fftbins,
        None => true,
    };

    let sym = Some(!fftbins);
    //try:
    //    beta = float(window)
    //except (TypeError, ValueError) as e:
    //    args = ()
    //    if isinstance(window, tuple):
    //        winstr = window[0]
    //        if len(window) > 1:
    //            args = window[1:]
    //    elif isinstance(window, str):
    //        if window in _needs_param:
    //            raise ValueError("The '" + window + "' window needs one or "
    //                             "more parameters -- pass a tuple.") from e
    //        else:
    //            winstr = window
    //    else:
    //        raise ValueError("%s as window type is not supported." %
    //                         str(type(window))) from e

    //    try:
    //        winfunc = _win_equiv[winstr]
    //    except KeyError as e:
    //        raise ValueError("Unknown window type.") from e

    //    if winfunc is dpss:
    //        params = (Nx,) + args + (None, sym)
    //    else:
    //        params = (Nx,) + args + (sym,)
    //else:
    //    winfunc = kaiser
    //    params = (Nx, beta, sym)
    match window {
        WindowKind::Value(beta) => {
            //let winfunc = kaiser;
            //let params = (Nx, beta, sym);
            todo!()
        }
        WindowKind::Kaiser(arg) => {
            //        winfunc = _win_equiv["kaiser"]
            //        if len(window) > 1:
            //            args = window[1:]

            //    if winfunc is dpss:
            //        params = (Nx,) + args + (None, sym)
            //    else:
            //        params = (Nx,) + args + (sym,)
            return kaiser(nx, arg, sym);
        }
        WindowKind::Hamming => {
            //        params = (Nx,) + args + (sym,)
            return hamming(nx, sym);
        }
        _ => {
            //    return winfunc(*params)
            todo!()
        }
    }
}

//_win_equiv_raw = {
//    ('barthann', 'brthan', 'bth'): (barthann, False),
//    ('bartlett', 'bart', 'brt'): (bartlett, False),
//    ('blackman', 'black', 'blk'): (blackman, False),
//    ('blackmanharris', 'blackharr', 'bkh'): (blackmanharris, False),
//    ('bohman', 'bman', 'bmn'): (bohman, False),
//    ('boxcar', 'box', 'ones',
//        'rect', 'rectangular'): (boxcar, False),
//    ('chebwin', 'cheb'): (chebwin, True),
//    ('cosine', 'halfcosine'): (cosine, False),
//    ('dpss',): (dpss, True),
//    ('exponential', 'poisson'): (exponential, True),
//    ('flattop', 'flat', 'flt'): (flattop, False),
//    ('gaussian', 'gauss', 'gss'): (gaussian, True),
//    ('general gaussian', 'general_gaussian',
//        'general gauss', 'general_gauss', 'ggs'): (general_gaussian, True),
//    ('hamming', 'hamm', 'ham'): (hamming, False),
//    ('hanning', 'hann', 'han'): (hann, False),
//    ('kaiser', 'ksr'): (kaiser, True),
//    ('nuttall', 'nutl', 'nut'): (nuttall, False),
//    ('parzen', 'parz', 'par'): (parzen, False),
//    ('taylor', 'taylorwin'): (taylor, False),
//    ('triangle', 'triang', 'tri'): (triang, False),
//    ('tukey', 'tuk'): (tukey, True),
//}
