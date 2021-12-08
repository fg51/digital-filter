use std::f64::consts::E;
use std::f64::consts::PI;

use num_complex::Complex64 as Complex;
use rustfft::FftPlanner;

use super::service::{linspace, polyval, sp_fft::next_fast_len};

//def freqz(b, a=1, worN=512, whole=False, plot=None, fs=2*pi, include_nyquist=False):
//Returns
//-------
//w : ndarray
//    The frequencies at which `h` was computed, in the same units as `fs`.
//    By default, `w` is normalized to the range [0, pi) (radians/sample).
//h : ndarray
//    The frequency response, as complex numbers.
//Examples
//--------
//>>> from scipy import signal
//>>> b = signal.firwin(80, 0.5, window=('kaiser', 8))
//>>> w, h = signal.freqz(b)

//>>> import matplotlib.pyplot as plt
//>>> fig, ax1 = plt.subplots()
//>>> ax1.set_title('Digital filter frequency response')

//>>> ax1.plot(w, 20 * np.log10(abs(h)), 'b')
//>>> ax1.set_ylabel('Amplitude [dB]', color='b')
//>>> ax1.set_xlabel('Frequency [rad/sample]')

//>>> ax2 = ax1.twinx()
//>>> angles = np.unwrap(np.angle(h))
//>>> ax2.plot(w, angles, 'g')
//>>> ax2.set_ylabel('Angle (radians)', color='g')
//>>> ax2.grid()
//>>> ax2.axis('tight')
//>>> plt.show()

//Broadcasting Examples

//Suppose we have two FIR filters whose coefficients are stored in the
//rows of an array with shape (2, 25). For this demonstration, we'll
//use random data:

//>>> np.random.seed(42)
//>>> b = np.random.rand(2, 25)

//To compute the frequency response for these two filters with one call
//to `freqz`, we must pass in ``b.T``, because `freqz` expects the first
//axis to hold the coefficients. We must then extend the shape with a
//trivial dimension of length 1 to allow broadcasting with the array
//of frequencies.  That is, we pass in ``b.T[..., np.newaxis]``, which has
//shape (25, 2, 1):

//>>> w, h = signal.freqz(b.T[..., np.newaxis], worN=1024)
//>>> w.shape
//(1024,)
//>>> h.shape
//(2, 1024)

//Now, suppose we have two transfer functions, with the same numerator
//coefficients ``b = [0.5, 0.5]``. The coefficients for the two denominators
//are stored in the first dimension of the 2-D array  `a`::

//    a = [   1      1  ]
//        [ -0.25, -0.5 ]

//>>> b = np.array([0.5, 0.5])
//>>> a = np.array([[1, 1], [-0.25, -0.5]])

//Only `a` is more than 1-D. To make it compatible for
//broadcasting with the frequencies, we extend it with a trivial dimension
//in the call to `freqz`:

//>>> w, h = signal.freqz(b, a[..., np.newaxis], worN=1024)
//>>> w.shape
//(1024,)
//>>> h.shape
//(2, 1024)

pub fn freqz(
    b: &[f64],
    a: &[f64],
    num: Option<usize>,
    whole: Option<bool>,
    fs: Option<f64>,
    include_nyquist: Option<bool>,
) -> (Vec<f64>, Vec<Complex>) {
    //if worN is None:
    //    # For backwards compatibility
    //    worN = 512
    let num = match num {
        Some(num) => num,
        None => 512,
    };
    let whole = match whole {
        Some(whole) => whole,
        None => false,
    };

    let w_or_n = 512;

    let fs = match fs {
        Some(fs) => fs,
        None => 2. * PI,
    };

    let _include_nyquist = match include_nyquist {
        Some(include_nyquist) => include_nyquist,
        None => false,
    };

    //h = None

    let (w, h): (Vec<f64>, Option<Vec<Complex>>) = if is_int_type(w_or_n) {
        //    N = operator.index(worN)
        //    del worN
        //    if N < 0:
        //        raise ValueError('worN must be nonnegative, got %s' % (N,))
        let lastpoint = if whole { 2. * PI } else { PI };
        // if include_nyquist is true and whole is false, w should include end point
        let endpoint = Some(false);
        let w = linspace(0., lastpoint, num, endpoint);
        if is_plan_fft(a.len(), b.len(), num) == true {
            // if N is fast, 2 * N will be fast, too, so no need to check
            let n_fft = if whole { num } else { num * 2 };
            //        if np.isrealobj(b) and np.isrealobj(a):
            //            fft_func = sp_fft.rfft
            //        else:
            //            fft_func = sp_fft.fft
            // h = fft_func(b, n_fft, axis = 0)[0..N];

            let mut planner = FftPlanner::<f64>::new();
            let fft = planner.plan_fft_forward(n_fft);
            let mut h = vec![Complex::new(0., 0.); n_fft];
            for (i, v) in b.iter().enumerate() {
                h[i].re = *v;
            }
            fft.process(&mut h);
            let mut h: Vec<Complex> = (0..num).map(|i| h[i]).collect();

            let denom = Complex::new(a[0], 0.);
            for i in 0..h.len() {
                h[i] = h[i] / denom;
            }
            // if fft_func is sp_fft.rfft and whole:
            //     # exclude DC and maybe Nyquist (no need to use axis_reverse
            //     # here because we can build reversal with the truncation)
            //     stop = -1 if n_fft % 2 == 1 else -2
            //     h_flip = slice(stop, 0, -1)
            //     h = np.concatenate((h, h[h_flip].conj()))
            (w, Some(h))
        } else {
            (w, None)
        }
    } else {
        //    w = atleast_1d(worN)
        //    del worN
        //    w = 2*pi*w/fs
        todo!();
    };

    let h = match h {
        Some(h) => h,
        None => {
            let zm1: Vec<Complex> = w.iter().map(|x| Complex::new(0., -(*x)).expf(E)).collect(); //zm1 = exp(-1j * w)
            let h_numer = polyval(&zm1, b, Some(false)); //    h = (npp_polyval(zm1, b, tensor=False) /
            let h_denom = polyval(&zm1, a, Some(false)); //         npp_polyval(zm1, a, tensor=False))

            let mut h = vec![];
            for i in 0..h_numer.len() {
                h.push(h_numer[i] / h_denom[i]);
            }
            h
        }
    };

    return (normalized_omega(&w, fs), h);
}

// range: [0,pi), unit: [radians/sample]
fn normalized_omega(w: &[f64], fs: f64) -> Vec<f64> {
    w.into_iter().map(|x| x * fs / (2. * PI)).collect()
}

fn is_int_type(_w_or_n: usize) -> bool {
    true
}

fn is_plan_fft(num_of_a: usize, num_of_b: usize, num: usize) -> bool {
    if (num_of_a == 1) == false {
        return false;
    }
    if (num >= num_of_b) == false {
        return false;
    }
    if (next_fast_len(num) == num) == false {
        return false;
    }
    //    //& ((b.ndim == 1) | (b.shape[-1] == 1))
    return true;
}
