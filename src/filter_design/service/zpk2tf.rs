//Return polynomial transfer function representation from zeros and poles

//Parameters
//----------
//z : array_like
//    Zeros of the transfer function.
//p : array_like
//    Poles of the transfer function.
//k : float
//    System gain.

//Returns
//-------
//b : ndarray
//    Numerator polynomial coefficients.
//a : ndarray
//    Denominator polynomial coefficients.

//"""

use super::polynominal::Poly;
use num_complex::Complex64 as Complex;

//def zpk2tf(z, p, k):
pub fn zpk2tf(z: &[f64], p: &[Complex], k: f64) -> (Vec<f64>, Vec<f64>) {
    //println("{:?}", z);

    //z = atleast_1d(z)
    //k = atleast_1d(k)
    //if len(z.shape) > 1:
    let b: Vec<f64> = if false {
        //    temp = poly(z[0])
        //    b = np.empty((z.shape[0], z.shape[1] + 1), temp.dtype.char)
        //    if len(k) == 1:
        //        k = [k[0]] * z.shape[0]
        //    for i in range(z.shape[0]):
        //        b[i] = k[i] * poly(z[i])
        todo!()
    } else {
        f64::poly(z).into_iter().map(|x| k * x).collect()
    };
    let a = Complex::poly(p); // atleast_1d(poly(p))

    // Use real output if possible. Copied from numpy.poly, since
    // we can't depend on a specific version of numpy.
    //if issubclass(b.dtype.type, numpy.complexfloating):
    //    # if complex roots are all complex conjugates, the roots are real.
    //    roots = numpy.asarray(z, complex)
    //    pos_roots = numpy.compress(roots.imag > 0, roots)
    //    neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
    //    if len(pos_roots) == len(neg_roots):
    //        if numpy.all(numpy.sort_complex(neg_roots) ==
    //                     numpy.sort_complex(pos_roots)):
    //            b = b.real.copy()

    //if issubclass(a.dtype.type, numpy.complexfloating):
    //    # if complex roots are all complex conjugates, the roots are real.
    //    roots = numpy.asarray(p, complex)
    //    pos_roots = numpy.compress(roots.imag > 0, roots)
    //    neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
    //    if len(pos_roots) == len(neg_roots):
    //        if numpy.all(numpy.sort_complex(neg_roots) ==
    //                     numpy.sort_complex(pos_roots)):
    //            a = a.real.copy()

    let a = a.into_iter().map(|x| x.re).collect();
    return (b, a);
}
