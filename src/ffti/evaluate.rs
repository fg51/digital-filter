use crate::values::{Complex, Direction, PI};

/*
 * In-place FFT butterfly algorithm
 *
 * input:
 *     A[] = array of N shuffled complex values where N is a power of 2
 * output:
 *     A[] = the DFT of input A[]
 *
 * for r = 1 to log2(N)
 *     m = 2^r
 *     Wm = exp(−j2π/m)
 *     for n = 0 to N-1 by m
 *         Wmk = 1
 *         for k = 0 to m/2 - 1
 *             u = A[n + k]
 *             t = Wmk * A[n + k + m/2]
 *             A[n + k]       = u + t
 *             A[n + k + m/2] = u - t
 *             Wmk = Wmk * Wm
 *
 * For inverse FFT, use Wm = exp(+j2π/m)
 */
pub fn ffti_evaluate(xs: &mut [Complex], log2n: usize, direction: Direction) {
    let num = 1 << log2n;
    let theta_2pi = match direction {
        Direction::Forward => -PI,
        Direction::Inverse => PI,
    } * 2.;

    for i in 1..=log2n {
        let m = 1 << i;
        let md2 = m >> 1;
        let theta = theta_2pi / m as f64;
        let wm = Complex::new(theta.cos(), theta.sin());
        //    for (n = 0; n < N; n += m)
        let mut j = 0;
        while j < num {
            let mut wmk = Complex::new(1., 0.);
            for k in 0..md2 {
                let ie = j + k;
                let io = ie + md2;
                let u = xs[ie].clone();
                let t = Complex::new(complex_mul_re(&wmk, &xs[io]), complex_mul_im(&wmk, &xs[io]));
                xs[ie] = u + t;
                xs[io] = u - t;
                wmk = Complex::new(complex_mul_re(&wmk, &wm), complex_mul_im(&wmk, &wm));
            }
            j += m;
        }
    }
}

fn complex_mul_re(x: &Complex, y: &Complex) -> f64 {
    x.re * y.re - x.im * y.im
}

fn complex_mul_im(x: &Complex, y: &Complex) -> f64 {
    x.re * y.im + x.im * y.re
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::copy_shuffle::ffti_copy_shuffle_f;

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

    #[test]
    fn ffti_evaluate_performs_8pt_dft() {
        let ts = [
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ];
        let expected_ws = [
            Complex::new(5.0, 0.000000),
            Complex::new(0.0, -2.414214),
            Complex::new(1.0, 0.000000),
            Complex::new(0.0, -0.414214),
            Complex::new(1.0, 0.000000),
            Complex::new(0.0, 0.414214),
            Complex::new(1.0, 0.000000),
            Complex::new(0.0, 2.414214),
        ];
        //    complex_f data[8];
        //    int i;

        let mut xs = vec![Complex::new(0., 0.); 8];
        ffti_copy_shuffle_f(&ts, &mut xs, 3); // 3 = log2(8)

        ffti_evaluate(&mut xs, 3, Direction::Forward); // 3 = log2(8)

        for i in 0..8 {
            assert!(
                (xs[i].re - expected_ws[i].re).abs() < 1E-5,
                "[{}]re:: expect: {}, actual: {}",
                i,
                expected_ws[i].re,
                xs[i].re
            );
            assert!(
                (xs[i].im - expected_ws[i].im).abs() < 1E-5,
                "[{}]im:: expect: {}, actual: {}",
                i,
                expected_ws[i].im,
                xs[i].im,
            );
        }
    }
}
//START_TEST (ffti_evaluate_f_performs_8pt_inverse_dft)
//{
//    complex_f F_w[8] = {
//        { 5.0f , 0.000000f },
//        { 0.0f , -2.414214f },
//        { 1.0f , 0.000000f },
//        { 0.0f , -0.414214f },
//        { 1.0f , 0.000000f },
//        { 0.0f , 0.414214f },
//        { 1.0f , 0.000000f },
//        { 0.0f , 2.414214f }
//    };
//    complex_f expected_f_t_N[8] = {
//        { 8.0f , 0.0f },
//        { 8.0f , 0.0f },
//        { 8.0f , 0.0f },
//        { 8.0f , 0.0f },
//        { 8.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f }
//    };
//    complex_f data[8];
//    int i;
//
//    ffti_copy_shuffle_f(F_w, data, 3);  /* 3 = log2(8) */
//
//    ffti_evaluate_f(data, 3, FFT_INVERSE);  /* 3 = log2(8) */
//
//    for (i = 0; i < 8; i++)
//    {
//        ck_assert_flt_eq_eps(data[i].re, expected_f_t_N[i].re, FLOAT_EQ_TOLERANCE);
//        ck_assert_flt_eq_eps(data[i].im, expected_f_t_N[i].im, FLOAT_EQ_TOLERANCE);
//    }
//}
//END_TEST
//
//START_TEST (ffti_evaluate_f_performs_16pt_dft)
//{
//    complex_f f_t[16] = {
//        { 1.0f , 0.0f },
//        { 1.0f , 0.0f },
//        { 1.0f , 0.0f },
//        { 1.0f , 0.0f },
//        { 1.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f }
//    };
//    complex_f expected_F_w[16] = {
//        { 5.000000f , 0.000000f },
//        { 3.013670f , -3.013670f },
//        { 0.000000f , -2.414214f },
//        { -0.248303f , -0.248303f },
//        { 1.000000f , 0.000000f },
//        { 0.834089f , -0.834089f },
//        { 0.000000f , -0.414214f },
//        { 0.400544f , 0.400544f },
//        { 1.000000f , 0.000000f },
//        { 0.400544f , -0.400544f },
//        { 0.000000f , 0.414214f },
//        { 0.834089f , 0.834089f },
//        { 1.000000f , 0.000000f },
//        { -0.248303f , 0.248303f },
//        { 0.000000f , 2.414214f },
//        { 3.013670f , 3.013670f }
//    };
//    complex_f data[16];
//    int i;
//
//    ffti_copy_shuffle_f(f_t, data, 4);  /* 4 = log2(16) */
//
//    ffti_evaluate_f(data, 4, FFT_FORWARD);  /* 4 = log2(16) */
//
//    for (i = 0; i < 16; i++)
//    {
//        ck_assert_flt_eq_eps(data[i].re, expected_F_w[i].re, FLOAT_EQ_TOLERANCE);
//        ck_assert_flt_eq_eps(data[i].im, expected_F_w[i].im, FLOAT_EQ_TOLERANCE);
//    }
//}
//END_TEST
//
//START_TEST (ffti_evaluate_f_performs_16pt_inverse_dft)
//{
//    complex_f F_w[16] = {
//        { 5.000000f , 0.000000f },
//        { 3.013670f , -3.013670f },
//        { 0.000000f , -2.414214f },
//        { -0.248303f , -0.248303f },
//        { 1.000000f , 0.000000f },
//        { 0.834089f , -0.834089f },
//        { 0.000000f , -0.414214f },
//        { 0.400544f , 0.400544f },
//        { 1.000000f , 0.000000f },
//        { 0.400544f , -0.400544f },
//        { 0.000000f , 0.414214f },
//        { 0.834089f , 0.834089f },
//        { 1.000000f , 0.000000f },
//        { -0.248303f , 0.248303f },
//        { 0.000000f , 2.414214f },
//        { 3.013670f , 3.013670f }
//    };
//    complex_f expected_f_t_N[16] = {
//        { 16.0f , 0.0f },
//        { 16.0f , 0.0f },
//        { 16.0f , 0.0f },
//        { 16.0f , 0.0f },
//        { 16.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f },
//        { 0.0f , 0.0f }
//    };
//    complex_f data[16];
//    int i;
//
//    ffti_copy_shuffle_f(F_w, data, 4);  /* 4 = log2(16) */
//
//    ffti_evaluate_f(data, 4, FFT_INVERSE);  /* 4 = log2(16) */
//
//    for (i = 0; i < 16; i++)
//    {
//        ck_assert_flt_eq_eps(data[i].re, expected_f_t_N[i].re, FLOAT_EQ_TOLERANCE);
//        ck_assert_flt_eq_eps(data[i].im, expected_f_t_N[i].im, FLOAT_EQ_TOLERANCE);
//    }
//}
