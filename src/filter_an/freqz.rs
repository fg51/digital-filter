use crate::errors::Result;
use crate::polyval::polyval_complex;
use crate::values::Complex;

// int DSPL_API freqz(double* b, double* a, int ord, double* w, int n, complex_t *h)
pub fn freqz(
    b: &[f64],
    a: Option<&[f64]>,
    order: usize,
    w: &[f64],
    n: usize,
    h: &mut [Complex],
) -> Result<()> {
    //    if(!b || !w || !h)
    //        return ERROR_PTR;
    //    if(ord<0)
    //        return ERROR_FILTER_ORD;
    //    if(n<1)
    //        return ERROR_SIZE;

    let bc: Vec<Complex> = (0..order + 1).map(|i| Complex::new(b[i], 0.)).collect();

    match a {
        Some(a) => {
            // IIR filter
            let ac: Vec<Complex> = (0..order + 1).map(|i| Complex::new(a[i], 0.)).collect();
            for i in 0..n {
                let jw = [Complex::new(w[i].cos(), -(w[i].sin()))];
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
                let mag = {
                    let mag = denom.norm_sqr(); // ABSSQR(den);
                    1.0 / mag // ERROR_DIV_ZERO
                };
                h[i] = Complex::new(
                    (numer.re * denom.re + numer.im * denom.im) * mag,
                    (numer.im * denom.re - numer.re * denom.im) * mag,
                );
            }
        }
        None => {
            // FIR filter
            for i in 0..n {
                let jw = [Complex::new(w[i].cos(), -(w[i].sin()))];
                polyval_complex(&bc, order, &jw, 1, &mut h[i..])?;
                //if(res != RES_OK)
                //    goto exit_label;
            }
        }
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
