use crate::values::Complex;

use crate::errors::Result;

pub fn complex_mul_re(x: &Complex, y: &Complex) -> f64 {
    x.re * y.re - x.im * y.im
}

pub fn complex_mul_im(x: &Complex, y: &Complex) -> f64 {
    x.re * y.im + x.im * y.re
}

#[allow(dead_code)]
fn cmplx2re(
    xs: &[Complex],
    num: usize,
    re: Option<&mut [f64]>,
    im: Option<&mut [f64]>,
) -> Result<()> {
    //if(!x)
    //    return ERROR_PTR;
    //if(n < 1)
    //    return ERROR_SIZE;

    if let Some(re) = re {
        for i in 0..num {
            re[i] = xs[i].re;
        }
    }
    if let Some(im) = im {
        for i in 0..num {
            im[i] = xs[i].im;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_mul_re_should_multiplies_two_complex_numbers() {
        let z1 = Complex::new(2.0, 3.0);
        let z2 = Complex::new(4.0, 5.0);

        let z3 = Complex::new(complex_mul_re(&z1, &z2), 0.0);
        assert_eq!(z3.re, -7.0);
    }

    //START_TEST (test_complex_mul_re_multiplies_two_more_complex_f_numbers)
    //{
    //    complex_f z1 = { 2.0f, -3.0f };
    //    complex_f z2 = { 4.0f, 5.0f };
    //    complex_f z3;
    //
    //    z3.re = complex_mul_re(z1.re, z1.im, z2.re, z2.im);
    //    z3.im = 0.0f;
    //
    //    ck_assert_flt_eq(z3.re, 23.0f);
    //}
    //END_TEST

    //START_TEST (test_complex_mul_re_multiplies_complex_f_and_complex_d_numbers)
    //{
    //    complex_f z1 = { 2.0f, -3.0f };
    //    complex_d z2 = { 4.0, 5.0 };
    //    complex_d z3;
    //
    //    z3.re = complex_mul_re(z1.re, z1.im, z2.re, z2.im);
    //    z3.im = 0.0;
    //
    //    ck_assert_dbl_eq(z3.re, 23.0);
    //}
    //END_TEST

    #[test]
    fn complex_mul_im_should_multiplies_two_complex_f_numbers() {
        let z1 = Complex::new(2.0, 3.0);
        let z2 = Complex::new(4.0, 5.0);
        //    complex_f z3;

        let z3 = Complex::new(0.0, complex_mul_im(&z1, &z2));
        assert_eq!(z3.im, 22.0);
    }

    //START_TEST (test_complex_mul_im_multiplies_two_more_complex_f_numbers)
    //{
    //    complex_f z1 = { 2.0f, -3.0f };
    //    complex_f z2 = { 4.0f, 5.0f };
    //    complex_f z3;
    //
    //    z3.re = 0.0f;
    //    z3.im = complex_mul_im(z1.re, z1.im, z2.re, z2.im);
    //
    //    ck_assert_flt_eq(z3.im, -2.0f);
    //}
    //END_TEST
    //
    //START_TEST (test_complex_mul_im_multiplies_complex_f_and_complex_d_numbers)
    //{
    //    complex_f z1 = { 2.0f, -3.0f };
    //    complex_d z2 = { 4.0, 5.0 };
    //    complex_d z3;
    //
    //    z3.re = 0.0;
    //    z3.im = complex_mul_im(z1.re, z1.im, z2.re, z2.im);
    //
    //    ck_assert_dbl_eq(z3.im, -2.0);
    //}
}
