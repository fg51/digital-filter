mod polyval;
pub use polyval::polyval;

mod linspace;
pub use linspace::linspace;

mod prod;
pub use prod::prod;

mod lp2lp_zpk;
pub use lp2lp_zpk::lp2lp_zpk;

mod zpk2tf;
pub use zpk2tf::zpk2tf;

pub mod sp_fft {
    pub fn next_fast_len(n: usize) -> usize {
        good_size_cmplx(n)
        //good_size_real(n)
    }

    /* returns the smallest composite of 2, 3, 5 which is >= n */
    #[allow(dead_code)]
    pub fn good_size_real(n: usize) -> usize {
        if n <= 6 {
            return n;
        }

        let mut best_fac = 2 * n;
        let mut f5 = 1;
        while f5 < best_fac {
            let mut x = f5; //  size_t x = f5;
            while x < n {
                x <<= 1; //x *= 2;
            }
            loop {
                if x < n {
                    x *= 3;
                } else if x > n {
                    if x < best_fac {
                        best_fac = x;
                    }
                    if (x & 1) == 0b1 {
                        break;
                    };
                    x >>= 1;
                } else {
                    return n;
                }
            }
            f5 *= 5;
        }
        return best_fac;
    }

    // returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n
    fn good_size_cmplx(n: usize) -> usize {
        if n <= 12 {
            return n;
        }

        let mut best_fac: usize = 2 * n;
        let mut f11 = 1;
        while f11 < best_fac {
            let mut f11_7 = f11;
            while f11_7 < best_fac {
                let mut f11_7_5 = f11_7;
                while f11_7_5 < best_fac {
                    let mut x = f11_7_5;
                    while x < n {
                        x *= 2;
                    }
                    loop {
                        if x < n {
                            x *= 3;
                        } else if n < x {
                            if x < best_fac {
                                best_fac = x;
                            }
                            if (x & 0b1) == 0b1 {
                                break;
                            }
                            x >>= 1;
                        } else {
                            return n;
                        }
                    }
                    f11_7_5 *= 5;
                }
                f11_7 *= 7;
            }
            f11 *= 11;
        }
        return best_fac;
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn it_works() {
            let (n, expect) = (512, 512);
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (1944, 1944); // (2.pow(3)) * (3.pow(5))
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (513, 525);
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (514, 525);
            assert_eq!(next_fast_len(n), expect);

            let (n, expect) = (526, 528);
            assert_eq!(next_fast_len(n), expect);
        }
    }
}
