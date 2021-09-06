use crate::values::Complex;

pub fn ffti_shuffle(xs: &mut [Complex], log2n: usize) {
    let num = 1 << log2n;
    let nd2 = num >> 1; /* N/2 = number range midpoint */
    let nm1 = num - 1; /* N-1 = digit mask */
    let mut j = 0; // index for next element swap location
    for i in 0..num {
        if j > i {
            let t = xs[i];
            xs[i] = xs[j];
            xs[j] = t;
        }

        // find least significant zero bit
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

        // Toggle bits with bit-reverse mask
        let bits = nm1 & !(if mszb == 0 { usize::MAX } else { mszb - 1 });
        j ^= bits;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shuffle_works_with_four_values() {
        /*
         * 0: 00 -> 00 : 0
         * 1: 01 -> 10 : 2
         * 2: 10 -> 01 : 1
         * 3: 11 -> 11 : 3
         */

        let mut xs = [
            Complex::new(1., -1.),
            Complex::new(2., -2.),
            Complex::new(3., -3.),
            Complex::new(4., -4.),
        ];

        ffti_shuffle(&mut xs, 2); // 2 = log2(4)

        assert_eq!(xs[0].re, 1.);
        assert_eq!(xs[0].im, -1.);

        assert_eq!(xs[1].re, 3.);
        assert_eq!(xs[1].im, -3.);

        assert_eq!(xs[2].re, 2.);
        assert_eq!(xs[2].im, -2.);

        assert_eq!(xs[3].re, 4.);
        assert_eq!(xs[3].im, -4.);
    }

    #[test]
    fn shuffle_works_with_eight_values() {
        /*
         * 0: 000 -> 000 : 0
         * 1: 001 -> 100 : 4
         * 2: 010 -> 010 : 2
         * 3: 011 -> 110 : 6
         * 4: 100 -> 001 : 1
         * 5: 101 -> 101 : 5
         * 6: 110 -> 011 : 3
         * 7: 111 -> 111 : 7
         */

        let mut xs = [
            Complex::new(1., -1.),
            Complex::new(2., -2.),
            Complex::new(3., -3.),
            Complex::new(4., -4.),
            Complex::new(5., -5.),
            Complex::new(6., -6.),
            Complex::new(7., -7.),
            Complex::new(8., -8.),
        ];

        ffti_shuffle(&mut xs, 3); // 3 = log2(8)
        assert_eq!(xs[0].re, 1.);
        assert_eq!(xs[0].im, -1.);

        assert_eq!(xs[1].re, 5.);
        assert_eq!(xs[1].im, -5.);

        assert_eq!(xs[2].re, 3.);
        assert_eq!(xs[2].im, -3.);

        assert_eq!(xs[3].re, 7.);
        assert_eq!(xs[3].im, -7.);

        assert_eq!(xs[4].re, 2.);
        assert_eq!(xs[4].im, -2.);

        assert_eq!(xs[5].re, 6.);
        assert_eq!(xs[5].im, -6.);

        assert_eq!(xs[6].re, 4.);
        assert_eq!(xs[6].im, -4.);

        assert_eq!(xs[7].re, 8.);
        assert_eq!(xs[7].im, -8.);
    }
}

