pub type Result<T> = std::result::Result<T, String>;

use crate::polyval::poly_z2a_cmplx;
use crate::values::Complex;

pub fn filter_zp2ab(
    zeros: &[Complex],
    nz: usize,
    poles: &[Complex],
    np: usize,
    ord: usize,
    b: &mut [f64],
    a: &mut [f64],
) -> Result<()> {
    //    if(!z || !p || !b || !a)
    //        return ERROR_PTR;
    //    if(nz < 0 || np < 0)
    //        return ERROR_SIZE;
    //    if(nz > ord || np > ord)
    //        return ERROR_POLY_ORD;

    let bcc = poly_z2a_cmplx(zeros, nz, ord)?;
    for i in 0..ord + 1 {
        b[i] = bcc[i].re;
    }

    let acc = poly_z2a_cmplx(poles, np, ord)?;
    for i in 0..ord + 1 {
        a[i] = acc[i].re;
    }
    return Ok(());
}
