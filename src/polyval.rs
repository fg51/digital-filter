use crate::conv::conv_cmplx;
use crate::service::complex::{complex_mul_im, complex_mul_re};
use crate::values::Complex;

pub type Result<T> = std::result::Result<T, String>;

pub fn poly_z2a_cmplx(z: &[Complex], nz: usize, ord: usize) -> Result<Vec<Complex>> {
    //  if(!z || !a)
    //      return ERROR_PTR;
    //  if(nz < 0)
    //      return ERROR_SIZE;
    //  if(nz > ord || ord < 1)
    //      return ERROR_POLY_ORD;

    //let mut xs = [Complex::new(0., 0.), Complex::new(1.0, 0.)];

    let mut a = vec![Complex::new(0., 0.); ord + 1];
    a[0].re = 1.0;

    let mut index = 1;
    for k in 0..nz {
        let xs = [-z[k], Complex::new(1.0, 0.)];
        let c = conv_cmplx(&a, index, &xs, xs.len())?;
        for i in 0..c.len() {
            a[i] = c[i];
        }
        index += 1;
    }

    //  return RES_OK;
    return Ok(a);
}

//int DSPL_API polyval_cmplx(complex_t* a, int ord, complex_t* x, int n, complex_t* y)
pub fn polyval_complex(
    a: &[Complex],
    order: usize,
    x: &[Complex],
    num: usize,
    y: &mut [Complex],
) -> Result<()> {
    //int k, m;

    //if(!a || !x || !y)
    //    return ERROR_PTR;
    //if(ord<0)
    //    return ERROR_POLY_ORD;
    //if(n<1)
    //    return ERROR_SIZE;

    for i in 0..num {
        y[i] = a[order];
        let mut j = (order - 1) as isize;
        while j > -1 {
            let t = Complex::new(complex_mul_re(&y[i], &x[i]), complex_mul_im(&y[i], &x[i]));
            y[i] = t + a[j as usize];
            j -= 1;
        }
    }
    //return RES_OK;
    Ok(())
}
