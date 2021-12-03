//Compute the frequency response of a digital filter.
//
//Given the M-order numerator `b` and N-order denominator `a` of a digital
//filter, compute its frequency response::
//
//             jw                 -jw              -jwM
//    jw    B(e  )    b[0] + b[1]e    + ... + b[M]e
// H(e  ) = ------ = -----------------------------------
//             jw                 -jw              -jwN
//          A(e  )    a[0] + a[1]e    + ... + a[N]e
//
//Parameters
//----------
//b : array_like
//    Numerator of a linear filter. If `b` has dimension greater than 1,
//    it is assumed that the coefficients are stored in the first dimension,
//    and ``b.shape[1:]``, ``a.shape[1:]``, and the shape of the frequencies
//    array must be compatible for broadcasting.
//a : array_like
//    Denominator of a linear filter. If `b` has dimension greater than 1,
//    it is assumed that the coefficients are stored in the first dimension,
//    and ``b.shape[1:]``, ``a.shape[1:]``, and the shape of the frequencies
//    array must be compatible for broadcasting.
//worN : {None, int, array_like}, optional
//    If a single integer, then compute at that many frequencies (default is
//    N=512). This is a convenient alternative to::
//
//        np.linspace(0, fs if whole else fs/2, N, endpoint=include_nyquist)
//
//    Using a number that is fast for FFT computations can result in
//    faster computations (see Notes).
//
//    If an array_like, compute the response at the frequencies given.
//    These are in the same units as `fs`.
//whole : bool, optional
//    Normally, frequencies are computed from 0 to the Nyquist frequency,
//    fs/2 (upper-half of unit-circle). If `whole` is True, compute
//    frequencies from 0 to fs. Ignored if worN is array_like.
//plot : callable
//    A callable that takes two arguments. If given, the return parameters
//    `w` and `h` are passed to plot. Useful for plotting the frequency
//    response inside `freqz`.
//fs : float, optional
//    The sampling frequency of the digital system. Defaults to 2*pi
//    radians/sample (so w is from 0 to pi).
//
//    .. versionadded:: 1.2.0
//include_nyquist : bool, optional
//    If `whole` is False and `worN` is an integer, setting `include_nyquist` to True
//    will include the last frequency (Nyquist frequency) and is otherwise ignored.
//
//    .. versionadded:: 1.5.0
//
//Returns
//-------
//w : ndarray
//    The frequencies at which `h` was computed, in the same units as `fs`.
//    By default, `w` is normalized to the range [0, pi) (radians/sample).
//h : ndarray
//    The frequency response, as complex numbers.
//
//See Also
//--------
//freqz_zpk
//sosfreqz
//
//Notes
//-----
//Using Matplotlib's :func:`matplotlib.pyplot.plot` function as the callable
//for `plot` produces unexpected results, as this plots the real part of the
//complex transfer function, not the magnitude.
//Try ``lambda w, h: plot(w, np.abs(h))``.
//
//A direct computation via (R)FFT is used to compute the frequency response
//when the following conditions are met:
//
//1. An integer value is given for `worN`.
//2. `worN` is fast to compute via FFT (i.e.,
//   `next_fast_len(worN) <scipy.fft.next_fast_len>` equals `worN`).
//3. The denominator coefficients are a single value (``a.shape[0] == 1``).
//4. `worN` is at least as long as the numerator coefficients
//   (``worN >= b.shape[0]``).
//5. If ``b.ndim > 1``, then ``b.shape[-1] == 1``.
//
//For long FIR filters, the FFT approach can have lower error and be much
//faster than the equivalent direct polynomial calculation.
//
//Examples
//--------
//>>> from scipy import signal
//>>> b = signal.firwin(80, 0.5, window=('kaiser', 8))
//>>> w, h = signal.freqz(b)
//
//>>> import matplotlib.pyplot as plt
//>>> fig, ax1 = plt.subplots()
//>>> ax1.set_title('Digital filter frequency response')
//
//>>> ax1.plot(w, 20 * np.log10(abs(h)), 'b')
//>>> ax1.set_ylabel('Amplitude [dB]', color='b')
//>>> ax1.set_xlabel('Frequency [rad/sample]')
//
//>>> ax2 = ax1.twinx()
//>>> angles = np.unwrap(np.angle(h))
//>>> ax2.plot(w, angles, 'g')
//>>> ax2.set_ylabel('Angle (radians)', color='g')
//>>> ax2.grid()
//>>> ax2.axis('tight')
//>>> plt.show()
//
//Broadcasting Examples
//
//Suppose we have two FIR filters whose coefficients are stored in the
//rows of an array with shape (2, 25). For this demonstration, we'll
//use random data:
//
//>>> np.random.seed(42)
//>>> b = np.random.rand(2, 25)
//
//To compute the frequency response for these two filters with one call
//to `freqz`, we must pass in ``b.T``, because `freqz` expects the first
//axis to hold the coefficients. We must then extend the shape with a
//trivial dimension of length 1 to allow broadcasting with the array
//of frequencies.  That is, we pass in ``b.T[..., np.newaxis]``, which has
//shape (25, 2, 1):
//
//>>> w, h = signal.freqz(b.T[..., np.newaxis], worN=1024)
//>>> w.shape
//(1024,)
//>>> h.shape
//(2, 1024)
//
//Now, suppose we have two transfer functions, with the same numerator
//coefficients ``b = [0.5, 0.5]``. The coefficients for the two denominators
//are stored in the first dimension of the 2-D array  `a`::
//
//    a = [   1      1  ]
//        [ -0.25, -0.5 ]
//
//>>> b = np.array([0.5, 0.5])
//>>> a = np.array([[1, 1], [-0.25, -0.5]])
//
//Only `a` is more than 1-D. To make it compatible for
//broadcasting with the frequencies, we extend it with a trivial dimension
//in the call to `freqz`:
//
//>>> w, h = signal.freqz(b, a[..., np.newaxis], worN=1024)
//>>> w.shape
//(1024,)
//>>> h.shape
//(2, 1024)

use std::f64::consts::PI;

mod service;
use service::linspace;

mod stub;

pub enum WorN<'a> {
    W(&'a [f64]),
    N(usize),
}

//def freqz(b, a=1, worN=512, whole=False, plot=None, fs=2*pi, include_nyquist=False):
pub fn freqz(
    b: &[f64],
    a: Option<&[f64]>,
    w_or_n: Option<WorN>,
    whole: Option<bool>,
    fs: Option<f64>,
    include_nyquist: Option<()>,
) -> (Vec<f64>, Vec<f64>) {
    todo!();
    //    //let default_a = [1.];
    //    //b = atleast_1d(b)
    //    let a = match a {
    //        Some(a) => a, //a = atleast_1d(a)
    //        None => &[1.],
    //    };
    //
    //    //if worN is None:
    //    //    worN = 512
    //    // For backwards compatibility
    //    let w_or_n = match w_or_n {
    //        Some(w_or_n) => w_or_n,
    //        None => WorN::N(512),
    //    };
    //
    //    let fs = match fs {
    //        Some(fs) => fs,
    //        None => 2. * PI,
    //    };
    //
    //    //h = None
    //
    //    //if _is_int_type(worN):
    //    let (w, h): (Vec<f64>, Option<Vec<f64>>) = match w_or_n {
    //        WorN::N(number) => {
    //            //    N = operator.index(worN)
    //            //    del worN
    //            //    if N < 0:
    //            //        raise ValueError('worN must be nonnegative, got %s' % (N,))
    //            //let lastpoint = 2. * PI if whole else PI;
    //            let lastpoint = match whole {
    //                Some(true) => 2. * PI,
    //                Some(false) | None => PI,
    //            };
    //            // if include_nyquist is true and whole is false, w should include end point
    //            //let w = np.linspace(0, lastpoint, N, endpoint=include_nyquist and not whole);
    //            let w = linspace(
    //                0.,
    //                lastpoint,
    //                number,
    //                match include_nyquist {
    //                    None => None,
    //                    Some(_) => match whole {
    //                        None | Some(false) => Some(true),
    //                        Some(true) => None,
    //                    },
    //                },
    //            );
    //            if (a.len() == 1)
    //                    & (number >= b.len()) // & (number >= b.shape[0])
    //                    & (sp_fft::next_fast_len(number) == number)
    //            //& ((b.ndim == 1) | (b.shape[-1] == 1))
    //            {
    //                //    # if N is fast, 2 * N will be fast, too, so no need to check
    //                //    n_fft = N if whole else N * 2
    //                let n_fft = match whole {
    //                    Some(true) => number,
    //                    Some(false) | None => number * 2,
    //                };
    //                //    if np.isrealobj(b) and np.isrealobj(a):
    //                //        fft_func = sp_fft.rfft
    //                //    else:
    //                //        fft_func = sp_fft.fft
    //                //    h = fft_func(b, n=n_fft, axis=0)[:N]
    //                let h: Vec<f64> = vec![];
    //                let h: Vec<f64> = h.into_iter().map(|x| x / a[0]).collect();
    //                //    if fft_func is sp_fft.rfft and whole:
    //                let is_rfft = true;
    //                if is_rfft {
    //                    if let Some(whole) = whole {
    //                        if whole == true {
    //                            // exclude DC and maybe Nyquist (no need to use axis_reverse
    //                            // here because we can build reversal with the truncation)
    //                            let stop = if (n_fft & 0b1) == 1 { -1 } else { -2 };
    //                            let h_flip = [stop, 0, -1]; // slice(stop, 0, -1)
    //                                                        // # xs = [1,2,3,4,5]
    //                                                        // # xs[slice(-1,0,-1)] == [5,4,3,2]
    //                                                        // # xs[slice(-2,0,-1)] == [4,3,2]
    //                                                        //let h = np.concatenate((h, h[h_flip].conj()));
    //                            h.extend(h[h_flip].conj());
    //                            (w, Some(h))
    //                        }
    //                    }
    //                }
    //                //    if b.ndim > 1:
    //                //        # Last axis of h has length 1, so drop it.
    //                //        h = h[..., 0]
    //                //        # Rotate the first axis of h to the end.
    //                //        h = np.rollaxis(h, 0, h.ndim)
    //                todo!();
    //            } else {
    //                (w, None)
    //            }
    //        }
    //        WorN::W(w) => {
    //            //else:
    //            //    w = atleast_1d(worN)
    //            //    del worN
    //            let w = w.into_iter().map(|x| 2. * PI * x / fs).collect();
    //            (w, None)
    //        }
    //    };
    //
    //    let h = match h {
    //        Some(h) => h,
    //        None => {
    //            // still need to compute using freqs w
    //            //    zm1 = exp(-1j * w)
    //            //    h = (npp_polyval(zm1, b, tensor=False) /
    //            //         npp_polyval(zm1, a, tensor=False))
    //            todo!();
    //        }
    //    };
    //
    //    let w = w.into_iter().map(|x| x * fs / (2. * PI)).collect();
    //
    //    //if plot is not None:
    //    //    plot(w, h)
    //
    //    return (w, h);
}

mod sp_fft {
    pub fn next_fast_len(n: usize) -> usize {
        good_size_cmplx(n)
        //good_size_real(n)
    }

    /* returns the smallest composite of 2, 3, 5 which is >= n */
    pub fn good_size_real(n: usize) -> usize {
        if n <= 6 {
            return n;
        }

        let mut best_fac = 2 * n;
        let mut f5 = 1;
        while f5 < best_fac {
            let mut x = f5; //  size_t x = f5;
            while x < n {
                x <<= 1; //x *= 2;
            }
            loop {
                if x < n {
                    x *= 3;
                } else if x > n {
                    if x < best_fac {
                        best_fac = x;
                    }
                    if (x & 1) == 0b1 {
                        break;
                    };
                    x >>= 1;
                } else {
                    return n;
                }
            }
            f5 *= 5;
        }
        return best_fac;
    }

    // returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n
    fn good_size_cmplx(n: usize) -> usize {
        if n <= 12 {
            return n;
        }

        let mut best_fac: usize = 2 * n;
        let mut f11 = 1;
        while f11 < best_fac {
            let mut f11_7 = f11;
            while f11_7 < best_fac {
                let mut f11_7_5 = f11_7;
                while f11_7_5 < best_fac {
                    let mut x = f11_7_5;
                    while x < n {
                        x *= 2;
                    }
                    loop {
                        if x < n {
                            x *= 3;
                        } else if n < x {
                            if x < best_fac {
                                best_fac = x;
                            }
                            if (x & 0b1) == 0b1 {
                                break;
                            }
                            x >>= 1;
                        } else {
                            return n;
                        }
                    }
                    f11_7_5 *= 5;
                }
                f11_7 *= 7;
            }
            f11 *= 11;
        }
        return best_fac;
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn it_works() {
            let (n, expect) = (512, 512);
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (1944, 1944); // (2.pow(3)) * (3.pow(5))
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (513, 525);
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (514, 525);
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (526, 528);
            assert_eq!(next_fast_len(n), expect);
        }
    }
}
