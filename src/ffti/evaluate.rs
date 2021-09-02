use crate::service::complex::{complex_mul_im, complex_mul_re};
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

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::shuffle::ffti_copy_shuffle;

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

        let mut xs = ffti_copy_shuffle(&ts, 3); // 3 = log2(8)

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
