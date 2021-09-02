use crate::values::Complex;
/*
 * Basic Bit-Reversal Scheme:
 *
 * The incrementing pattern operations used here correspond
 * to the logic operations of a synchronous counter.
 *
 * Incrementing a binary number simply flips a sequence of
 * least-significant bits, for example from 0111 to 1000.
 * So in order to compute the next bit-reversed index, we
 * have to flip a sequence of most-significant bits.
 */
pub fn ffti_copy_shuffle_f(src: &[Complex], dst: &mut [Complex], log2n: usize) {
    let num = 1 << log2n;
    let nd2 = num >> 1; /* N/2 = number range midpoint */
    let nm1 = num - 1; /* N-1 = digit mask */

    let mut j = 0; // index for next destination element
    for i in 0..num {
        dst[j] = src[i];

        // Find least significant zero bit
        let lszb = (!i) & (i + 1);

        /*
         * Use division to bit-reverse the single bit so that we now have
         * the most significant zero bit
         *
         * N = 2^r = 2^(m+1)
         * Nd2 = N/2 = 2^m
         * if lszb = 2^k, where k is within the range of 0...m, then
         *     mszb = Nd2 / lszb
         *          = 2^m / 2^k
         *          = 2^(m-k)
         *          = bit-reversed value of lszb
         */
        let mszb = nd2 / lszb;

        /*
         * Toggle bits with bit-reverse mask
         */
        //let bits = nm1 & !(if mszb == 0 { usize::MAX } else { mszb - 1 });
        let bits = nm1 & !(mszb.wrapping_sub(1));
        j ^= bits;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ffti_copy_shuffle_f_copies_four_shuffled_values() {
        /*
         * 0: 00 -> 00 : 0
         * 1: 01 -> 10 : 2
         * 2: 10 -> 01 : 1
         * 3: 11 -> 11 : 3
         */

        let src = [
            Complex::new(1., -1.),
            Complex::new(2., -2.),
            Complex::new(3., -3.),
            Complex::new(4., -4.),
        ];

        let mut dst = src.clone();

        ffti_copy_shuffle_f(&src, &mut dst, 2); // 2 = log2(4)

        assert_eq!(dst[0].re, 1.);
        assert_eq!(dst[0].im, -1.);

        assert_eq!(dst[1].re, 3.);
        assert_eq!(dst[1].im, -3.);

        assert_eq!(dst[2].re, 2.);
        assert_eq!(dst[2].im, -2.);

        assert_eq!(dst[3].re, 4.);
        assert_eq!(dst[3].im, -4.);
    }
}

//
//START_TEST (ffti_copy_shuffle_f_copies_eight_shuffled_values)
//{
//    /*
//     * 0: 000 -> 000 : 0
//     * 1: 001 -> 100 : 4
//     * 2: 010 -> 010 : 2
//     * 3: 011 -> 110 : 6
//     * 4: 100 -> 001 : 1
//     * 5: 101 -> 101 : 5
//     * 6: 110 -> 011 : 3
//     * 7: 111 -> 111 : 7
//     */
//
//    complex_f src[8] = {
//        { 1.f, -1.f},
//        { 2.f, -2.f},
//        { 3.f, -3.f},
//        { 4.f, -4.f},
//        { 5.f, -5.f},
//        { 6.f, -6.f},
//        { 7.f, -7.f},
//        { 8.f, -8.f}
//    };
//    complex_f dst[8];
//    int i;
//
//    for (i = 0; i < 8; i++)
//    {
//        dst[i] = src[i];
//    }
//
//    ffti_copy_shuffle_f(src, dst, 3);  /* 3 = log2(8) */
//
//    ck_assert_flt_eq(dst[0].re, 1.f);
//    ck_assert_flt_eq(dst[0].im, -1.f);
//
//    ck_assert_flt_eq(dst[1].re, 5.f);
//    ck_assert_flt_eq(dst[1].im, -5.f);
//
//    ck_assert_flt_eq(dst[2].re, 3.f);
//    ck_assert_flt_eq(dst[2].im, -3.f);
//
//    ck_assert_flt_eq(dst[3].re, 7.f);
//    ck_assert_flt_eq(dst[3].im, -7.f);
//
//    ck_assert_flt_eq(dst[4].re, 2.f);
//    ck_assert_flt_eq(dst[4].im, -2.f);
//
//    ck_assert_flt_eq(dst[5].re, 6.f);
//    ck_assert_flt_eq(dst[5].im, -6.f);
//
//    ck_assert_flt_eq(dst[6].re, 4.f);
//    ck_assert_flt_eq(dst[6].im, -4.f);
//
//    ck_assert_flt_eq(dst[7].re, 8.f);
//    ck_assert_flt_eq(dst[7].im, -8.f);
//}
//END_TEST
//
//START_TEST (ffti_copy_shuffle_f_copies_sixteen_shuffled_values)
//{
//    /*
//     * 00: 0000 -> 0000 : 00
//     * 01: 0001 -> 1000 : 08
//     * 02: 0010 -> 0100 : 04
//     * 03: 0011 -> 1100 : 12
//     * 04: 0100 -> 0010 : 02
//     * 05: 0101 -> 1010 : 10
//     * 06: 0110 -> 0110 : 06
//     * 07: 0111 -> 1110 : 14
//     * 08: 1000 -> 0001 : 01
//     * 09: 1001 -> 1001 : 09
//     * 10: 1010 -> 0101 : 05
//     * 11: 1011 -> 1101 : 13
//     * 12: 1100 -> 0011 : 03
//     * 13: 1101 -> 1011 : 11
//     * 14: 1110 -> 0111 : 07
//     * 15: 1111 -> 1111 : 15
//     */
//
//    complex_f src[16] = {
//        { 1.f, -1.f},
//        { 2.f, -2.f},
//        { 3.f, -3.f},
//        { 4.f, -4.f},
//        { 5.f, -5.f},
//        { 6.f, -6.f},
//        { 7.f, -7.f},
//        { 8.f, -8.f},
//        { 9.f, -9.f},
//        { 10.f, -10.f},
//        { 11.f, -11.f},
//        { 12.f, -12.f},
//        { 13.f, -13.f},
//        { 14.f, -14.f},
//        { 15.f, -15.f},
//        { 16.f, -16.f}
//    };
//    complex_f dst[16];
//    int i;
//
//    for (i = 0; i < 16; i++)
//    {
//        dst[i] = src[i];
//    }
//
//    ffti_copy_shuffle_f(src, dst, 4);  /* 4 = log2(16) */
//
//    ck_assert_flt_eq(dst[0].re, 1.f);
//    ck_assert_flt_eq(dst[0].im, -1.f);
//
//    ck_assert_flt_eq(dst[1].re, 9.f);
//    ck_assert_flt_eq(dst[1].im, -9.f);
//
//    ck_assert_flt_eq(dst[2].re, 5.f);
//    ck_assert_flt_eq(dst[2].im, -5.f);
//
//    ck_assert_flt_eq(dst[3].re, 13.f);
//    ck_assert_flt_eq(dst[3].im, -13.f);
//
//    ck_assert_flt_eq(dst[4].re, 3.f);
//    ck_assert_flt_eq(dst[4].im, -3.f);
//
//    ck_assert_flt_eq(dst[5].re, 11.f);
//    ck_assert_flt_eq(dst[5].im, -11.f);
//
//    ck_assert_flt_eq(dst[6].re, 7.f);
//    ck_assert_flt_eq(dst[6].im, -7.f);
//
//    ck_assert_flt_eq(dst[7].re, 15.f);
//    ck_assert_flt_eq(dst[7].im, -15.f);
//
//    ck_assert_flt_eq(dst[8].re, 2.f);
//    ck_assert_flt_eq(dst[8].im, -2.f);
//
//    ck_assert_flt_eq(dst[9].re, 10.f);
//    ck_assert_flt_eq(dst[9].im, -10.f);
//
//    ck_assert_flt_eq(dst[10].re, 6.f);
//    ck_assert_flt_eq(dst[10].im, -6.f);
//
//    ck_assert_flt_eq(dst[11].re, 14.f);
//    ck_assert_flt_eq(dst[11].im, -14.f);
//
//    ck_assert_flt_eq(dst[12].re, 4.f);
//    ck_assert_flt_eq(dst[12].im, -4.f);
//
//    ck_assert_flt_eq(dst[13].re, 12.f);
//    ck_assert_flt_eq(dst[13].im, -12.f);
//
//    ck_assert_flt_eq(dst[14].re, 8.f);
//    ck_assert_flt_eq(dst[14].im, -8.f);
//
//    ck_assert_flt_eq(dst[15].re, 16.f);
//    ck_assert_flt_eq(dst[15].im, -16.f);
//}
//END_TEST
//

