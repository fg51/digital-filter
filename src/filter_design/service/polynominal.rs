//"""
//Find the coefficients of a polynomial with the given sequence of roots.

//.. note::
//   This forms part of the old polynomial API. Since version 1.4, the
//   new polynomial API defined in `numpy.polynomial` is preferred.
//   A summary of the differences can be found in the
//   :doc:`transition guide </reference/routines.polynomials>`.

//Returns the coefficients of the polynomial whose leading coefficient
//is one for the given sequence of zeros (multiple roots must be included
//in the sequence as many times as their multiplicity; see Examples).
//A square matrix (or array, which will be treated as a matrix) can also
//be given, in which case the coefficients of the characteristic polynomial
//of the matrix are returned.

//Parameters
//----------
//seq_of_zeros : array_like, shape (N,) or (N, N)
//    A sequence of polynomial roots, or a square array or matrix object.

//Returns
//-------
//c : ndarray
//    1D array of polynomial coefficients from highest to lowest degree:

//    ``c[0] * x**(N) + c[1] * x**(N-1) + ... + c[N-1] * x + c[N]``
//    where c[0] always equals 1.

//Raises
//------
//ValueError
//    If input is the wrong shape (the input must be a 1-D or square
//    2-D array).

//See Also
//--------
//polyval : Compute polynomial values.
//roots : Return the roots of a polynomial.
//polyfit : Least squares polynomial fit.
//poly1d : A one-dimensional polynomial class.

//Notes
//-----
//Specifying the roots of a polynomial still leaves one degree of
//freedom, typically represented by an undetermined leading
//coefficient. [1]_ In the case of this function, that coefficient -
//the first one in the returned array - is always taken as one. (If
//for some reason you have one other point, the only automatic way
//presently to leverage that information is to use ``polyfit``.)

//The characteristic polynomial, :math:`p_a(t)`, of an `n`-by-`n`
//matrix **A** is given by

//    :math:`p_a(t) = \\mathrm{det}(t\\, \\mathbf{I} - \\mathbf{A})`,

//where **I** is the `n`-by-`n` identity matrix. [2]_

//References
//----------
//.. [1] M. Sullivan and M. Sullivan, III, "Algebra and Trignometry,
//   Enhanced With Graphing Utilities," Prentice-Hall, pg. 318, 1996.

//.. [2] G. Strang, "Linear Algebra and Its Applications, 2nd Edition,"
//   Academic Press, pg. 182, 1980.

//Examples
//--------
//Given a sequence of a polynomial's zeros:

//>>> np.poly((0, 0, 0)) # Multiple root example
//array([1., 0., 0., 0.])

//The line above represents z**3 + 0*z**2 + 0*z + 0.

//>>> np.poly((-1./2, 0, 1./2))
//array([ 1.  ,  0.  , -0.25,  0.  ])

//The line above represents z**3 - z/4

//>>> np.poly((np.random.random(1)[0], 0, np.random.random(1)[0]))
//array([ 1.        , -0.77086955,  0.08618131,  0.        ]) # random

//Given a square array object:

//>>> P = np.array([[0, 1./3], [-1./2, 0]])
//>>> np.poly(P)
//array([1.        , 0.        , 0.16666667])

//Note how in all cases the leading coefficient is always 1.

//def poly(seq_of_zeros):
//seq_of_zeros = atleast_1d(seq_of_zeros)
//sh = seq_of_zeros.shape

//if len(sh) == 2 and sh[0] == sh[1] and sh[0] != 0:
//    seq_of_zeros = eigvals(seq_of_zeros)
//elif len(sh) == 1:
//    dt = seq_of_zeros.dtype
//    # Let object arrays slip through, e.g. for arbitrary precision
//    if dt != object:
//        seq_of_zeros = seq_of_zeros.astype(mintypecode(dt.char))
//else:
//    raise ValueError("input must be 1d or non-empty square 2d array.")

//if len(seq_of_zeros) == 0:
//    return 1.0
//dt = seq_of_zeros.dtype
//a = ones((1,), dtype=dt)
//for zero in seq_of_zeros:
//    a = NX.convolve(a, array([1, -zero], dtype=dt), mode='full')

//if issubclass(a.dtype.type, NX.complexfloating):
//    # if complex roots are all complex conjugates, the roots are real.
//    roots = NX.asarray(seq_of_zeros, complex)
//    if NX.all(NX.sort(roots) == NX.sort(roots.conjugate())):
//        a = a.real.copy()

//return a

use num_complex::Complex64 as Complex;

pub trait Poly {
    type Value;
    // poly([-1,1,1,10]) = [1, -11, 9, 11, -10]
    // (x-1)(x+1)(x+1)(x+10) = 1*x^4 -11 * x^3 + 9 * x^2 + 11 * x^1 -10
    fn poly(xs: &[Self::Value]) -> Vec<Self::Value>;

    //"""
    //Returns the discrete, linear convolution of two one-dimensional sequences.

    //The convolution operator is often seen in signal processing, where it
    //models the effect of a linear time-invariant system on a signal [1]_.  In
    //probability theory, the sum of two independent random variables is
    //distributed according to the convolution of their individual
    //distributions.

    //If `v` is longer than `a`, the arrays are swapped before computation.

    //Parameters
    //----------
    //a : (N,) array_like
    //    First one-dimensional input array.
    //v : (M,) array_like
    //    Second one-dimensional input array.
    //mode : {'full', 'valid', 'same'}, optional
    //    'full':
    //      By default, mode is 'full'.  This returns the convolution
    //      at each point of overlap, with an output shape of (N+M-1,). At
    //      the end-points of the convolution, the signals do not overlap
    //      completely, and boundary effects may be seen.

    //    'same':
    //      Mode 'same' returns output of length ``max(M, N)``.  Boundary
    //      effects are still visible.

    //    'valid':
    //      Mode 'valid' returns output of length
    //      ``max(M, N) - min(M, N) + 1``.  The convolution product is only given
    //      for points where the signals overlap completely.  Values outside
    //      the signal boundary have no effect.

    //Returns
    //-------
    //out : ndarray
    //    Discrete, linear convolution of `a` and `v`.

    //See Also
    //--------
    //scipy.signal.fftconvolve : Convolve two arrays using the Fast Fourier
    //                           Transform.
    //scipy.linalg.toeplitz : Used to construct the convolution operator.
    //polymul : Polynomial multiplication. Same output as convolve, but also
    //          accepts poly1d objects as input.

    //Notes
    //-----
    //The discrete convolution operation is defined as

    //.. math:: (a * v)[n] = \\sum_{m = -\\infty}^{\\infty} a[m] v[n - m]

    //It can be shown that a convolution :math:`x(t) * y(t)` in time/space
    //is equivalent to the multiplication :math:`X(f) Y(f)` in the Fourier
    //domain, after appropriate padding (padding is necessary to prevent
    //circular convolution).  Since multiplication is more efficient (faster)
    //than convolution, the function `scipy.signal.fftconvolve` exploits the
    //FFT to calculate the convolution of large data-sets.

    //References
    //----------
    //.. [1] Wikipedia, "Convolution",
    //    https://en.wikipedia.org/wiki/Convolution

    //Examples
    //--------
    //Note how the convolution operator flips the second array
    //before "sliding" the two across one another:

    //>>> np.convolve([1, 2, 3], [0, 1, 0.5])
    //array([0. , 1. , 2.5, 4. , 1.5])

    //Only return the middle values of the convolution.
    //Contains boundary effects, where zeros are taken
    //into account:

    //>>> np.convolve([1,2,3],[0,1,0.5], 'same')
    //array([1. ,  2.5,  4. ])

    //The two arrays are of the same length, so there
    //is only one position where they completely overlap:

    //>>> np.convolve([1,2,3],[0,1,0.5], 'valid')
    //array([2.5])

    //"""
    //def convolve(a, v, mode='full'):
    fn convolve(a: &[Self::Value], v: &[Self::Value]) -> Vec<Self::Value>;

    fn one() -> Self::Value;
}

impl Poly for f64 {
    type Value = f64;
    fn poly(seq_of_zeros: &[Self::Value]) -> Vec<Self::Value> {
        //sh = seq_of_zeros.shape

        //if len(sh) == 2 and sh[0] == sh[1] and sh[0] != 0:
        //    seq_of_zeros = eigvals(seq_of_zeros)
        //elif len(sh) == 1:
        //    dt = seq_of_zeros.dtype
        //    # Let object arrays slip through, e.g. for arbitrary precision
        //    if dt != object:
        //        seq_of_zeros = seq_of_zeros.astype(mintypecode(dt.char))
        //else:
        //    raise ValueError("input must be 1d or non-empty square 2d array.")

        if seq_of_zeros.len() == 0 {
            return vec![Self::one()];
        }
        //dt = seq_of_zeros.dtype
        let mut a = vec![Self::one()]; //a = ones((1,), dtype=dt)
        for zero in seq_of_zeros {
            a = Self::convolve(&a, &[Self::one(), -zero]); // NX.convolve(a, [1., -zero], mode='full')
        }

        //if issubclass(a.dtype.type, NX.complexfloating):
        //    # if complex roots are all complex conjugates, the roots are real.
        //    roots = NX.asarray(seq_of_zeros, complex)
        //    if NX.all(NX.sort(roots) == NX.sort(roots.conjugate())):
        //        a = a.real.copy()

        return a;
    }

    fn convolve(a: &[Self::Value], v: &[Self::Value]) -> Vec<Self::Value> {
        //    a, v = array(a, copy=False, ndmin=1), array(v, copy=False, ndmin=1)
        let (a, v) = if v.len() > a.len() { (v, a) } else { (a, v) };

        if a.len() == 0 {
            todo!();
            //return Err(ErrorKind::ValueError("a cannot be empty".to_string()));
        }
        if v.len() == 0 {
            todo!();
            //return Err(ErrorKind::ValueError("v cannot be empty".to_string()));
        }
        if v.len() == 1 {
            let v = v[0];
            a.iter().map(|x| x * v).collect()
        } else {
            let a = {
                let mut a1 = vec![Default::default(); (v.len() - 1) + a.len() + (v.len() - 1)];
                for (i, &a0) in a.iter().enumerate() {
                    a1[i + 1] = a0;
                }
                a1
            };
            let mut xs = vec![];
            for i in 0..(a.len() - v.len() + 1) {
                let mut sum = Default::default();
                for (j, v1) in v.iter().rev().enumerate() {
                    sum += a[i + j] * v1;
                }
                xs.push(sum);
            }
            return xs;
        }
    }

    fn one() -> Self::Value {
        1.
    }
}

impl Poly for Complex {
    type Value = Complex;
    fn poly(seq_of_zeros: &[Self::Value]) -> Vec<Self::Value> {
        //sh = seq_of_zeros.shape

        //if len(sh) == 2 and sh[0] == sh[1] and sh[0] != 0:
        //    seq_of_zeros = eigvals(seq_of_zeros)
        //elif len(sh) == 1:
        //    dt = seq_of_zeros.dtype
        //    # Let object arrays slip through, e.g. for arbitrary precision
        //    if dt != object:
        //        seq_of_zeros = seq_of_zeros.astype(mintypecode(dt.char))
        //else:
        //    raise ValueError("input must be 1d or non-empty square 2d array.")

        if seq_of_zeros.len() == 0 {
            return vec![Self::one()];
        }
        //dt = seq_of_zeros.dtype
        let mut a = vec![Self::one()]; //a = ones((1,), dtype=dt)
        for zero in seq_of_zeros {
            a = Self::convolve(&a, &[Self::one(), -zero]); // NX.convolve(a, [1., -zero], mode='full')
        }

        //if issubclass(a.dtype.type, NX.complexfloating):
        //    # if complex roots are all complex conjugates, the roots are real.
        //    roots = NX.asarray(seq_of_zeros, complex)
        //    if NX.all(NX.sort(roots) == NX.sort(roots.conjugate())):
        //        a = a.real.copy()

        return a;
    }

    fn convolve(a: &[Self::Value], v: &[Self::Value]) -> Vec<Self::Value> {
        let (a, v) = if v.len() > a.len() { (v, a) } else { (a, v) };

        if a.len() == 0 {
            todo!();
            //return Err(ErrorKind::ValueError("a cannot be empty".to_string()));
        }
        if v.len() == 0 {
            todo!();
            //return Err(ErrorKind::ValueError("v cannot be empty".to_string()));
        }
        if v.len() == 1 {
            a.iter().map(|x| x * v[0]).collect()
        } else {
            let a = {
                let mut a1 = vec![Default::default(); (v.len() - 1) + a.len() + (v.len() - 1)];
                for (i, &a0) in a.iter().enumerate() {
                    a1[i + 1] = a0;
                }
                a1
            };
            let mut xs = vec![];
            for i in 0..(a.len() - v.len() + 1) {
                let mut sum = Default::default();
                for (j, v1) in v.iter().rev().enumerate() {
                    sum += a[i + j] * v1;
                }
                xs.push(sum);
            }
            return xs;
        }
    }

    fn one() -> Self::Value {
        Complex::new(1., 0.)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_f64_works() {
        let expect = [1., -11., 9., 11., -10.];
        let actual = f64::poly(&[-1., 1., 1., 10.]);
        for (e, a) in expect.into_iter().zip(actual) {
            assert_eq!(e, a);
        }
    }

    #[test]
    fn convolve_f64_works() {
        let expect = [0., 0.2, 1.2, 2.2, 3.2, 4.2, 4.];
        let actual = f64::convolve(&[0., 1., 2., 3., 4., 5.], &[0.2, 0.8]);
        for (e, a) in expect.into_iter().zip(actual) {
            assert!((e - a).abs() < 1E-5, "expect:{}, actual:{}", e, a);
        }
    }
}
