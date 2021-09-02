use crate::values::{Complex, Direction};

mod evaluate;
use evaluate::ffti_evaluate;

mod shuffle;
use shuffle::ffti_shuffle;

mod copy_shuffle;

pub fn ffti(xs: &mut [Complex], log2n: usize, direction: Direction) {
    ffti_shuffle(xs, log2n);
    ffti_evaluate(xs, log2n, direction);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ffti_should_performs_4pt_inplace_dft() {
        let mut xs = [
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        ];
        let expected_fw = [
            Complex::new(10.0, 0.0),
            Complex::new(-2.0, 2.0),
            Complex::new(-2.0, 0.0),
            Complex::new(-2.0, -2.0),
        ];

        ffti(&mut xs, 2, Direction::Forward); // 2 = log2(4)

        for i in 0..4 {
            assert!((xs[i].re - expected_fw[i].re).abs() < 1E-5);
            assert!((xs[i].im - expected_fw[i].im).abs() < 1E-5);
        }
    }

    //START_TEST (ffti_f_performs_4pt_inplace_inverse_DFT)
    //{
    //    complex_f data[4] = {
    //        { 10.0f , 0.0f },
    //        { -2.0f , 2.0f },
    //        { -2.0f , 0.0f },
    //        { -2.0f , -2.0f }
    //    };
    //    complex_f expected_f_t_N[4] = {
    //        { 4.0f , 0.0f },
    //        { 8.0f , 0.0f },
    //        { 12.0f , 0.0f },
    //        { 16.0f , 0.0f }
    //    };
    //    int i;
    //
    //    ffti_f(data, 2, FFT_INVERSE);  /* 2 = log2(4) */
    //
    //    for (i = 0; i < 4; i++)
    //    {
    //        ck_assert_flt_eq_eps(data[i].re, expected_f_t_N[i].re, FLOAT_EQ_TOLERANCE);
    //        ck_assert_flt_eq_eps(data[i].im, expected_f_t_N[i].im, FLOAT_EQ_TOLERANCE);
    //    }
    //}
    //END_TEST
    //
    //START_TEST (ffti_f_performs_8pt_inplace_DFT)
    //{
    //    complex_f data[8] = {
    //        { 1.0f , 0.0f },
    //        { 1.0f , 0.0f },
    //        { 1.0f , 0.0f },
    //        { 1.0f , 0.0f },
    //        { 1.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f }
    //    };
    //    complex_f expected_F_w[8] = {
    //        { 5.0f , 0.000000f },
    //        { 0.0f , -2.414214f },
    //        { 1.0f , 0.000000f },
    //        { 0.0f , -0.414214f },
    //        { 1.0f , 0.000000f },
    //        { 0.0f , 0.414214f },
    //        { 1.0f , 0.000000f },
    //        { 0.0f , 2.414214f }
    //    };
    //    int i;
    //
    //    ffti_f(data, 3, FFT_FORWARD);  /* 3 = log2(8) */
    //
    //    for (i = 0; i < 8; i++)
    //    {
    //        ck_assert_flt_eq_eps(data[i].re, expected_F_w[i].re, FLOAT_EQ_TOLERANCE);
    //        ck_assert_flt_eq_eps(data[i].im, expected_F_w[i].im, FLOAT_EQ_TOLERANCE);
    //    }
    //}
    //END_TEST
    //
    //START_TEST (ffti_f_performs_8pt_inplace_inverse_DFT)
    //{
    //    complex_f data[8] = {
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
    //    int i;
    //
    //    ffti_f(data, 3, FFT_INVERSE);  /* 3 = log2(8) */
    //
    //    for (i = 0; i < 8; i++)
    //    {
    //        ck_assert_flt_eq_eps(data[i].re, expected_f_t_N[i].re, FLOAT_EQ_TOLERANCE);
    //        ck_assert_flt_eq_eps(data[i].im, expected_f_t_N[i].im, FLOAT_EQ_TOLERANCE);
    //    }
    //}
    //END_TEST
    //
    //START_TEST (ffti_f_performs_32pt_inplace_DFT)
    //{
    //    complex_f data[32] = {
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
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f }
    //    };
    //    complex_f expected_F_w[32] = {
    //        { 5.000000f , 0.000000f },
    //        { 4.443241f , -1.840451f },
    //        { 3.013670f , -3.013670f },
    //        { 1.311956f , -3.167342f },
    //        { 0.000000f , -2.414214f },
    //        { -0.515005f , -1.243333f },
    //        { -0.248303f , -0.248303f },
    //        { 0.422747f , 0.175108f },
    //        { 1.000000f , 0.000000f },
    //        { 1.143707f , -0.473739f },
    //        { 0.834089f , -0.834089f },
    //        { 0.335425f , -0.809787f },
    //        { 0.000000f , -0.414214f },
    //        { 0.039197f , 0.094631f },
    //        { 0.400544f , 0.400544f },
    //        { 0.818731f , 0.339130f },
    //        { 1.000000f , 0.000000f },
    //        { 0.818731f , -0.339130f },
    //        { 0.400544f , -0.400544f },
    //        { 0.039197f , -0.094631f },
    //        { 0.000000f , 0.414214f },
    //        { 0.335425f , 0.809787f },
    //        { 0.834089f , 0.834089f },
    //        { 1.143707f , 0.473739f },
    //        { 1.000000f , 0.000000f },
    //        { 0.422747f , -0.175108f },
    //        { -0.248303f , 0.248303f },
    //        { -0.515005f , 1.243333f },
    //        { 0.000000f , 2.414214f },
    //        { 1.311956f , 3.167342f },
    //        { 3.013670f , 3.013670f },
    //        { 4.443241f , 1.840451f }
    //    };
    //    int i;
    //
    //    ffti_f(data, 5, FFT_FORWARD);  /* 5 = log2(32) */
    //
    //    for (i = 0; i < 32; i++)
    //    {
    //        ck_assert_flt_eq_eps(data[i].re, expected_F_w[i].re, FLOAT_EQ_TOLERANCE);
    //        ck_assert_flt_eq_eps(data[i].im, expected_F_w[i].im, FLOAT_EQ_TOLERANCE);
    //    }
    //}
    //END_TEST
    //
    //START_TEST (ffti_f_performs_32pt_inplace_inverse_DFT)
    //{
    //    complex_f data[32] = {
    //        { 5.000000f , 0.000000f },
    //        { 4.443241f , -1.840451f },
    //        { 3.013670f , -3.013670f },
    //        { 1.311956f , -3.167342f },
    //        { 0.000000f , -2.414214f },
    //        { -0.515005f , -1.243333f },
    //        { -0.248303f , -0.248303f },
    //        { 0.422747f , 0.175108f },
    //        { 1.000000f , 0.000000f },
    //        { 1.143707f , -0.473739f },
    //        { 0.834089f , -0.834089f },
    //        { 0.335425f , -0.809787f },
    //        { 0.000000f , -0.414214f },
    //        { 0.039197f , 0.094631f },
    //        { 0.400544f , 0.400544f },
    //        { 0.818731f , 0.339130f },
    //        { 1.000000f , 0.000000f },
    //        { 0.818731f , -0.339130f },
    //        { 0.400544f , -0.400544f },
    //        { 0.039197f , -0.094631f },
    //        { 0.000000f , 0.414214f },
    //        { 0.335425f , 0.809787f },
    //        { 0.834089f , 0.834089f },
    //        { 1.143707f , 0.473739f },
    //        { 1.000000f , 0.000000f },
    //        { 0.422747f , -0.175108f },
    //        { -0.248303f , 0.248303f },
    //        { -0.515005f , 1.243333f },
    //        { 0.000000f , 2.414214f },
    //        { 1.311956f , 3.167342f },
    //        { 3.013670f , 3.013670f },
    //        { 4.443241f , 1.840451f }
    //    };
    //    complex_f expected_f_t_N[32] = {
    //        { 32.0f , 0.0f },
    //        { 32.0f , 0.0f },
    //        { 32.0f , 0.0f },
    //        { 32.0f , 0.0f },
    //        { 32.0f , 0.0f },
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
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f },
    //        { 0.0f , 0.0f }
    //    };
    //    int i;
    //
    //    ffti_f(data, 5, FFT_INVERSE);  /* 5 = log2(32) */
    //
    //    for (i = 0; i < 32; i++)
    //    {
    //        ck_assert_flt_eq_eps(data[i].re, expected_f_t_N[i].re, FLOAT_EQ_TOLERANCE);
    //        ck_assert_flt_eq_eps(data[i].im, expected_f_t_N[i].im, FLOAT_EQ_TOLERANCE);
    //    }
    //}
}
