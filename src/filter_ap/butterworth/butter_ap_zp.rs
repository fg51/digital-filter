use crate::values::Complex;
use crate::values::PI;

pub type Result<T> = std::result::Result<T, String>;

//int DSPL_API butter_ap_zp(int ord, double rp, complex_t* z, int* nz,
pub fn butter_ap_zp(
    ord: usize,
    ripple: f64,
    _zeros: &mut [Complex],
    poles: &mut [Complex],
) -> Result<(usize, usize)> {
    //    if(rp < 0 || rp == 0)
    //        return ERROR_FILTER_RP;
    //    if(ord < 1)
    //        return ERROR_FILTER_ORD;
    //    if(!z || !p || !nz || !np)
    //        return ERROR_PTR;

    //let ep = sqrt(pow(10.0, ripple * 0.1) - 1.0);
    let ep = (10.0f64.powf(ripple * 0.1) - 1.0).sqrt();
    let r = ord % 2;
    let length = (ord - r) / 2;

    let alpha = ep.powf(-1.0 / ord as f64);
    let mut index = 0;
    if r > 0 {
        poles[index] = Complex::new(-alpha, 0.);
        index += 1;
    }
    for k in 0..length {
        let theta = PI * (2 * k + 1) as f64 / (2 * ord) as f64;
        poles[index] = Complex::new(-alpha * theta.sin(), alpha * theta.cos());
        poles[index + 1] = Complex::new(-alpha * theta.sin(), -alpha * theta.cos());
        index += 2;
    }
    return Ok((0, ord));
}

#[cfg(test)]
mod tests {
    use super::*;
    //use crate::values::Complex;

    #[test]
    fn butter_ap_zp_works() {
        let ord = 7;

        let mut zeros = vec![Complex::new(0., 0.); ord + 1]; // H(s) zeros vector
        let mut poles = vec![Complex::new(0., 0.); ord + 1]; // H(s) poles vector
        let ripple = 1.; // Magnitude ripple from 0 to 1 rad/s

        //      int res, k, nz, np;

        // Zeros and poles vectors calculation
        let (num_of_zeros, num_of_poles) =
            butter_ap_zp(ord, ripple, &mut zeros, &mut poles).unwrap();

        assert_eq!(num_of_zeros, 0);
        assert_eq!(num_of_poles, 7);
        //      if(res != RES_OK)
        //          printf("error code = 0x%8x\n", res);

        //      /* print H(s) zeros values */
        //      printf("Butterworth filter zeros: %d\n", nz);
        //      for(k = 0; k < nz; k++)
        //          printf("z[%2d] = %9.3f %+9.3f j\n", k, RE(z[k]), IM(z[k]));

        //      /* print H(s) poles values */
        //      printf("Butterworth filter poles: %d\n", np);
        //      for(k = 0; k < np; k++)
        //          printf("p[%2d] = %9.3f %+9.3f j\n", k, RE(p[k]), IM(p[k]));

        let excepts = excepts();
        for i in 0..7 {
            assert!((poles[i].re - excepts[i].re).abs() < 1E-3);
            assert!((poles[i].im - excepts[i].im).abs() < 1E-3);
        }
    }

    fn excepts() -> [Complex; 7] {
        [
            Complex::new(-1.101, 0.000),
            Complex::new(-0.245, 1.074),
            Complex::new(-0.245, -1.074),
            Complex::new(-0.687, 0.861),
            Complex::new(-0.687, -0.861),
            Complex::new(-0.992, 0.478),
            Complex::new(-0.992, -0.478),
        ]
    }
}
