//FIR filter design using the window method.
//
//This function computes the coefficients of a finite impulse response
//filter. The filter will have linear phase; it will be Type I if
//`numtaps` is odd and Type II if `numtaps` is even.
//
//Type II filters always have zero response at the Nyquist frequency, so a
//ValueError exception is raised if firwin is called with `numtaps` even and
//having a passband whose right end is at the Nyquist frequency.
//
//Parameters
//----------
//numtaps : int
//    Length of the filter (number of coefficients, i.e. the filter
//    order + 1).  `numtaps` must be odd if a passband includes the
//    Nyquist frequency.
//cutoff : float or 1-D array_like
//    Cutoff frequency of filter (expressed in the same units as `fs`)
//    OR an array of cutoff frequencies (that is, band edges). In the
//    latter case, the frequencies in `cutoff` should be positive and
//    monotonically increasing between 0 and `fs/2`. The values 0 and
//    `fs/2` must not be included in `cutoff`.
//width : float or None, optional
//    If `width` is not None, then assume it is the approximate width
//    of the transition region (expressed in the same units as `fs`)
//    for use in Kaiser FIR filter design. In this case, the `window`
//    argument is ignored.
//window : string or tuple of string and parameter values, optional
//    Desired window to use. See `scipy.signal.get_window` for a list
//    of windows and required parameters.
//pass_zero : {True, False, 'bandpass', 'lowpass', 'highpass', 'bandstop'}, optional
//    If True, the gain at the frequency 0 (i.e., the "DC gain") is 1.
//    If False, the DC gain is 0. Can also be a string argument for the
//    desired filter type (equivalent to ``btype`` in IIR design functions).
//
//    .. versionadded:: 1.3.0
//       Support for string arguments.
//scale : bool, optional
//    Set to True to scale the coefficients so that the frequency
//    response is exactly unity at a certain frequency.
//    That frequency is either:
//
//    - 0 (DC) if the first passband starts at 0 (i.e. pass_zero
//      is True)
//    - `fs/2` (the Nyquist frequency) if the first passband ends at
//      `fs/2` (i.e the filter is a single band highpass filter);
//      center of first passband otherwise
//
//nyq : float, optional
//    *Deprecated. Use `fs` instead.* This is the Nyquist frequency.
//    Each frequency in `cutoff` must be between 0 and `nyq`. Default
//    is 1.
//fs : float, optional
//    The sampling frequency of the signal. Each frequency in `cutoff`
//    must be between 0 and ``fs/2``.  Default is 2.
//
//Returns
//-------
//h : (numtaps,) ndarray
//    Coefficients of length `numtaps` FIR filter.
//
//Raises
//------
//ValueError
//    If any value in `cutoff` is less than or equal to 0 or greater
//    than or equal to ``fs/2``, if the values in `cutoff` are not strictly
//    monotonically increasing, or if `numtaps` is even but a passband
//    includes the Nyquist frequency.
//
//See Also
//--------
//firwin2
//firls
//minimum_phase
//remez

use crate::errors::{ErrorKind, Result};
use std::f64::consts::PI;

use super::service::get_fs;
use super::service::signaltools::window::get_window;
use super::values::WindowKind;

//def firwin(numtaps, cutoff, width=None, window='hamming', pass_zero=True,
//           scale=True, nyq=None, fs=None):
pub fn firwin(
    numtaps: usize,
    cutoff: &[f64],
    width: Option<usize>,
    window: Option<WindowKind>,
    pass_zero: Option<bool>,
    scale: Option<bool>,
    nyq: Option<f64>,
    fs: Option<f64>,
) -> Result<Vec<f64>> {
    let window = match window {
        Some(window) => window,
        None => WindowKind::Hamming,
    };
    let pass_zero = match pass_zero {
        Some(pass_zero) => pass_zero,
        None => true,
    };
    let scale = match scale {
        Some(scale) => scale,
        None => true,
    };
    let nyq = 0.5 * get_fs(fs, nyq);

    let cutoff: Vec<f64> = cutoff.iter().map(|x| *x / nyq).collect();

    // Check for invalid input.
    //if cutoff.ndim > 1:
    //    raise ValueError("The cutoff argument must be at most " "one-dimensional.")
    //if cutoff.size == 0:
    //    raise ValueError("At least one cutoff frequency must be given.")
    //if cutoff.min() <= 0 or cutoff.max() >= 1:
    if (cutoff.iter().fold(0. / 0., |m, v| v.min(m)) <= 0.)
        | (cutoff.iter().fold(0. / 0., |m, v| v.max(m)) >= 1.)
    {
        //    raise ValueError("Invalid cutoff frequency: frequencies must be "
        //                     "greater than 0 and less than fs/2.")
        return Err(ErrorKind::ValueError);
    }

    //if np.any(np.diff(cutoff) <= 0):
    //    raise ValueError("Invalid cutoff frequencies: the frequencies "
    //                     "must be strictly increasing.")

    // A width was given.  Find the beta parameter of the Kaiser window  and set `window`.  This overrides the value of `window` passed in.
    let window = match width {
        Some(width) => {
            //    atten = kaiser_atten(numtaps, float(width) / nyq)
            //    beta = kaiser_beta(atten)
            //    window = ('kaiser', beta)
            todo!()
        }
        None => window,
    };

    //if isinstance(pass_zero, str):
    //    if pass_zero in ('bandstop', 'lowpass'):
    //        if pass_zero == 'lowpass':
    //            if cutoff.size != 1:
    //                raise ValueError('cutoff must have one element if '
    //                                 'pass_zero=="lowpass", got %s'
    //                                 % (cutoff.shape,))
    //        elif cutoff.size <= 1:
    //            raise ValueError('cutoff must have at least two elements if '
    //                             'pass_zero=="bandstop", got %s'
    //                             % (cutoff.shape,))
    //        pass_zero = True
    //    elif pass_zero in ('bandpass', 'highpass'):
    //        if pass_zero == 'highpass':
    //            if cutoff.size != 1:
    //                raise ValueError('cutoff must have one element if '
    //                                 'pass_zero=="highpass", got %s'
    //                                 % (cutoff.shape,))
    //        elif cutoff.size <= 1:
    //            raise ValueError('cutoff must have at least two elements if '
    //                             'pass_zero=="bandpass", got %s'
    //                             % (cutoff.shape,))
    //        pass_zero = False
    //    else:
    //        raise ValueError('pass_zero must be True, False, "bandpass", '
    //                         '"lowpass", "highpass", or "bandstop", got '
    //                         '%s' % (pass_zero,))
    //pass_zero = bool(operator.index(pass_zero))  # ensure bool-like

    //pass_nyquist = bool(cutoff.size & 1) ^ pass_zero
    let pass_nyquist = ((cutoff.len() & 1) == 1) ^ pass_zero;
    //if pass_nyquist and (numtaps % 2) == 0:
    if pass_nyquist & ((numtaps & 1) == 0) {
        //    raise ValueError("A filter with an even number of coefficients must "
        //                     "have zero response at the Nyquist frequency.")
        return Err(ErrorKind::ValueError);
    }

    //# Insert 0 and/or 1 at the ends of cutoff so that the length of cutoff
    //# is even, and each pair in cutoff corresponds to passband.
    //cutoff = np.hstack(([0.0] * pass_zero, cutoff, [1.0] * pass_nyquist))
    let cutoff = {
        let mut xs = vec![];
        if pass_zero {
            xs.push(0.);
        }
        for i in cutoff {
            xs.push(i);
        }
        if pass_nyquist {
            xs.push(1.);
        }
        xs
    };

    // `bands` is a 2-D array; each row gives the left and right edges of a passband.
    //bands = cutoff.reshape(-1, 2)
    let bands = {
        let mut bands = vec![];
        for i in 0..cutoff.len() / 2 {
            bands.push((cutoff[2 * i], cutoff[2 * i + 1]));
        }
        bands
    };

    // Build up the coefficients.
    let alpha = 0.5 * ((numtaps - 1) as f64);
    //let m = np.arange(0, numtaps) - alpha
    let ms: Vec<f64> = (0..numtaps).map(|x| (x as f64) - alpha).collect();
    let mut h = vec![0.; ms.len()];
    //for left, right in bands:
    //    h += right * sinc(right * m)
    //    h -= left * sinc(left * m)

    for (left, right) in bands.iter() {
        for (i, m) in ms.iter().enumerate() {
            h[i] += right * sinc(right * m);
            h[i] -= left * sinc(left * m);
        }
    }

    // Get and apply the window function.
    //from .signaltools import get_window
    let win = get_window(window, numtaps, Some(false));

    for i in 0..h.len() {
        h[i] *= win[i]; //h *= win
    }

    //# Now handle scaling if desired.
    if scale {
        // Get the first passband.
        let (left, right) = bands[0];
        let scale_frequency = if left == 0. {
            0.0
        } else if right == 1. {
            1.0
        } else {
            0.5 * (left + right)
        };
        let cs: Vec<f64> = ms
            .iter()
            .map(|m| (PI * m * scale_frequency).cos())
            .collect();
        let mut hcs = vec![0.; h.len()];
        for i in 0..h.len() {
            hcs[i] = h[i] * cs[i];
        }
        let s: f64 = hcs.into_iter().sum();
        for i in 0..h.len() {
            h[i] /= s;
        }
    }

    return Ok(h);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_pass_from_0_to_f() {
        assert_eq!(2 + 2, 4);
        let numtaps = 3;
        let f = [0.1];
        //>>> signal.firwin(numtaps, f)
        let actual = firwin(numtaps, &f, None, None, None, None, None, None).unwrap();
        let expect = [0.06799017, 0.86401967, 0.06799017];

        for (i, ex) in expect.into_iter().enumerate() {
            assert!((ex - actual[i]).abs() < 1E-5, "Error at {}.", i);
        }

        //Use a specific window function:
        //
        //>>> signal.firwin(numtaps, f, window='nuttall')
        //array([  3.56607041e-04,   9.99286786e-01,   3.56607041e-04])
    }

    //High-pass ('stop' from 0 to f):
    //
    //>>> signal.firwin(numtaps, f, pass_zero=False)
    //array([-0.00859313,  0.98281375, -0.00859313])
    //
    //Band-pass:
    //
    //>>> f1, f2 = 0.1, 0.2
    //>>> signal.firwin(numtaps, [f1, f2], pass_zero=False)
    //array([ 0.06301614,  0.88770441,  0.06301614])
    //
    //Band-stop:
    //
    //>>> signal.firwin(numtaps, [f1, f2])
    //array([-0.00801395,  1.0160279 , -0.00801395])
    //
    //Multi-band (passbands are [0, f1], [f2, f3] and [f4, 1]):
    //
    //>>> f3, f4 = 0.3, 0.4
    //>>> signal.firwin(numtaps, [f1, f2, f3, f4])
    //array([-0.01376344,  1.02752689, -0.01376344])
    //
    //Multi-band (passbands are [f1, f2] and [f3,f4]):
    //
    //>>> signal.firwin(numtaps, [f1, f2, f3, f4], pass_zero=False)
    //array([ 0.04890915,  0.91284326,  0.04890915])
}

fn sinc(x: f64) -> f64 {
    let y = if x == 0. { PI * 1E-20 } else { PI * x };
    y.sin() / y
}
