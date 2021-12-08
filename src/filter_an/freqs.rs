pub type Result<T> = std::result::Result<T, String>;

use crate::polyval::polyval_complex;
use crate::values::Complex;

// int DSPL_API freqs(double* b, double* a, int ord, double* w, int n, complex_t *h)
pub fn freqs(
    b: &[f64],
    a: &[f64],
    order: usize,
    w: &[f64],
    n: usize,
    h: &mut [Complex],
) -> Result<()> {
    //    complex_t jw;
    //    complex_t *bc = NULL;
    //    complex_t *ac = NULL;
    //    complex_t num, den;
    //    double mag;
    //    int k;
    //    int res;

    //    if(!b || !a || !w || !h)
    //        return ERROR_PTR;
    //    if(ord<0)
    //        return ERROR_FILTER_ORD;
    //    if(n<1)
    //        return ERROR_SIZE;

    //let bc = vec![Complex::new(0., 0.); order + 1];
    //    res = re2cmplx(b, ord+1, bc);
    let bc: Vec<Complex> = (0..order + 1).map(|i| Complex::new(b[i], 0.)).collect();

    //    if( res!=RES_OK )
    //        goto exit_label;

    //    ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    //    res = re2cmplx(a, ord+1, ac);
    //    if( res!=RES_OK )
    //        goto exit_label;
    let ac: Vec<Complex> = (0..order + 1).map(|i| Complex::new(a[i], 0.)).collect();

    //let jw = Complex::new(0., 0.);
    for k in 0..n {
        //jw.im = w[k];
        let jw = [Complex::new(0., w[k])];

        let numer = {
            let mut numer = [Complex::new(0., 0.)];
            polyval_complex(&bc, order, &jw, 1, &mut numer)?;
            numer[0]
        };

        let denom = {
            let mut denom = [Complex::new(0., 0.)];
            polyval_complex(&ac, order, &jw, 1, &mut denom)?;
            denom[0]
        };
        let mag = denom.norm_sqr(); // ABSSQR(denom);
                                    //if mag == 0.0 {
                                    //    //            res = ERROR_DIV_ZERO;
                                    //    return Err("ERROR DIV ZERO");
                                    //}
        let mag = 1.0 / mag;
        //        RE(h[k]) = CMCONJRE(num, den) * mag;
        //        IM(h[k]) = CMCONJIM(num, den) * mag;
        h[k].re = (numer.re * denom.re + numer.im * denom.im) * mag;
        h[k].im = (numer.im * denom.re - numer.im * denom.im) * mag;
    }
    //    res = RES_OK;
    //exit_label:
    //    if(bc)
    //        free(bc);
    //    if(ac)
    //        free(ac);
    //    return res;
    todo!();
}
