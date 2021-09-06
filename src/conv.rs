use crate::errors::Result;
use crate::service::complex::{complex_mul_im, complex_mul_re};
use crate::values::Complex;

// int DSPL_API conv_cmplx(complex_t* a, int na, complex_t* b, int nb, complex_t* c)
pub fn conv_cmplx(
    a: &[Complex],
    num_of_a: usize,
    b: &[Complex],
    num_of_b: usize,
) -> Result<Vec<Complex>> {
    //int k;
    //int n;

    //complex_t *t;
    //size_t bufsize;

    //if(!a || !b || !c)
    //    return ERROR_PTR;
    //if(na < 1 || nb < 1)
    //    return ERROR_SIZE;

    let bufsize = num_of_a + num_of_b - 1;

    //if((a != c) && (b != c))
    //    t = c;
    //else
    //    t = (complex_t*)malloc(bufsize);

    let mut ts = vec![Complex::new(0., 0.); bufsize];

    for k in 0..num_of_a {
        for n in 0..num_of_b {
            ts[k + n].re += complex_mul_re(&a[k], &b[n]);
            ts[k + n].im += complex_mul_im(&a[k], &b[n]);
        }
    }

    //if(t!=c)
    //{
    //    memcpy(c, t, bufsize);
    //    free(t);
    //}

    //return RES_OK;
    return Ok(ts);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conv_cmplx_works() {
        //    void* handle;           /* DSPL handle        */
        //    handle = dspl_load();   /* Load DSPL function */
        let ac: Vec<Complex> = [(0.0, 1.0), (1.0, 1.0), (2.0, 2.0)]
            .iter()
            .map(|x| Complex::new(x.0, x.1))
            .collect();
        let bc: Vec<Complex> = [(3.0, 3.0), (4.0, 4.0), (5.0, 5.0), (6.0, 6.0)]
            .iter()
            .map(|x| Complex::new(x.0, x.1))
            .collect();
        //    complex_t cc[6];

        //    double ar[3] = {1.0, 2.0, 3.0};
        //    double br[4] = {3.0, -1.0, 2.0, 4.0};
        //    double cr[6];
        //
        //    int n;
        //
        //    printf("\nconv\n--------------------------------\n");
        //    conv(ar, 3, br, 4, cr);
        //    for(n = 0; n < 6; n++)
        //        printf("cr[%d] = %5.1f\n", n, cr[n]);

        //    printf("\nconv_cmplx\n--------------------------------\n");
        let cc = conv_cmplx(&ac, ac.len(), &bc, bc.len()).unwrap();
        //    for(n = 0; n < 6; n++)
        //        printf("cc[%d] = %5.1f%+5.1fj\n", n, RE(cc[n]),IM(cc[n]));
        //
        //    dspl_free(handle);      // free dspl handle
        //    return 0;
        //}

        for (i, expect) in vec![
            Complex::new(-3.0, 3.0),
            Complex::new(-4.0, 10.0),
            Complex::new(-5.0, 25.0),
            Complex::new(-6.0, 32.0),
            Complex::new(0.0, 32.0),
            Complex::new(0.0, 24.0),
        ]
        .iter()
        .enumerate()
        {
            assert!((cc[i].re - expect.re).abs() < 1E-1);
            assert!((cc[i].im - expect.im).abs() < 1E-1);
        }
    }
}
