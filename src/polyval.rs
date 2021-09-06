use crate::conv::conv_cmplx;
use crate::values::Complex;

use crate::errors::Result;

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
