#![allow(dead_code)]

mod create;
pub mod fft1;

use crate::complex::Complex;

pub struct FFT {
    w: Vec<Complex>,  //complex_t*  w;
    t0: Vec<Complex>, //complex_t*  t0;
    t1: Vec<Complex>, //complex_t*  t1;

    /* radix-2 twiddle factors vectors */
    w32: [Complex; 32],
    w64: [Complex; 64],
    w128: [Complex; 128],
    w256: [Complex; 256],
    w512: [Complex; 512],
    w1024: [Complex; 1024],
    w2048: [Complex; 2048],
    w4096: [Complex; 4096],
    n: usize,
}

impl FFT {
    pub const fn zero() -> Self {
        Self {
            w: vec![],
            t0: vec![],
            t1: vec![],
            w32: [Complex::zero(); 32],
            w64: [Complex::zero(); 64],
            w128: [Complex::zero(); 128],
            w256: [Complex::zero(); 256],
            w512: [Complex::zero(); 512],
            w1024: [Complex::zero(); 1024],
            w2048: [Complex::zero(); 2048],
            w4096: [Complex::zero(); 4096],

            n: 0,
        }
    }
}

impl FFT {
    //int DSPL_API fft(double* x, int n, fft_t* pfft, complex_t* y)
    pub fn fft(&self, xs: &[f64], num: usize) -> Vec<Complex> {
        //    int err;
        //
        //    if(!x || !pfft || !y)
        //        return ERROR_PTR;
        //    if(n<1)
        //        return ERROR_SIZE;

        //    err = fft_create(pfft, n);
        //    if(err != RES_OK)
        //        return err;

        //    re2cmplx(x, n, pfft->t1);
        let mut ys = vec![Complex::zero(); num];
        for i in 0..num {
            ys[i] = Complex::new(xs[i], 0.);
        }

        return fft_krn(&self.t1, ys, &self.w, num, 0);
        //    return fft_krn(pfft->t1, y, pfft, n, 0);
        //}
    }
}

use create::cal_n2;

//int fft_krn(complex_t* t0, complex_t* t1, fft_t* p, int n, int addr)
fn fft_krn(
    t0: &[Complex],
    mut t1: Vec<Complex>,
    pw: &[Complex],
    num: usize,
    addr: usize,
) -> Vec<Complex> {
    //    int n1, n2, k, m, i;
    //    complex_t *pw = p->w+addr;
    //    complex_t tmp;

    let n1 = cal_n2(num);

    //label_size:
    if n1 == 1 {
        for k in 0..num {
            // RE(t1[k]) = IM(t1[k]) = 0.0;
            t1[k] = Complex::zero();
            for m in 0..num {
                let i = (k * m) % num;
                //     RE(tmp) = CMRE(t0[m], pw[i]);
                //     IM(tmp) = CMIM(t0[m], pw[i]);
                //     RE(t1[k]) += RE(tmp);
                //     IM(t1[k]) += IM(tmp);
                t1[k] += Complex::new(
                    t0[m].re() * pw[addr + i].re() - t0[m].im() * pw[addr + i].im(),
                    t0[m].re() * pw[addr + i].im() + t0[m].im() * pw[addr + i].re(),
                );
            }
        }
        return t1;
    } else {
        let n2 = num / n1;

        if n2 > 1 {
            //            memcpy(t1, t0, n*sizeof(complex_t));
            //            matrix_transpose_cmplx(t1, n2, n1, t0);
            todo!();
        }

        if n1 == 4096 {
            //            for(k = 0; k < n2; k++)
            //                dft4096(t0+4096*k, t1+4096*k, p->w4096, p->w256);
            todo!();
        }

        if n1 == 2048 {
            //            for(k = 0; k < n2; k++)
            //                dft2048(t0+2048*k, t1+2048*k, p->w2048, p->w32, p->w64);
            todo!();
        }

        if n1 == 1024 {
            //            for(k = 0; k < n2; k++)
            //                dft1024(t0+1024*k, t1+1024*k, p->w1024, p->w32);
            todo!();
        }
        //
        //        if(n1 == 512)
        //            for(k = 0; k < n2; k++)
        //                dft512(t0+512*k, t1+512*k, p->w512, p->w32);
        //
        //        if(n1 == 256)
        //            for(k = 0; k < n2; k++)
        //                dft256(t0+256*k, t1+256*k, p->w256);
        //
        //        if(n1 == 128)
        //            for(k = 0; k < n2; k++)
        //                dft128(t0+128*k, t1+128*k, p->w128);
        //
        //        if(n1 == 64)
        //            for(k = 0; k < n2; k++)
        //                dft64(t0+64*k, t1+64*k, p->w64);
        //
        //        if(n1 == 32)
        //            for(k = 0; k < n2; k++)
        //                dft32(t0+32*k, t1+32*k, p->w32);
        //
        //        if(n1 == 16)
        //            for(k = 0; k < n2; k++)
        //                dft16(t0+16*k, t1+16*k);
        //
        //        if(n1 == 7)
        //            for(k = 0; k < n2; k++)
        //                dft7(t0+7*k, t1+7*k);
        //
        //        if(n1 == 8)
        //            for(k = 0; k < n2; k++)
        //                dft8(t0+8*k, t1+8*k);
        //
        //        if(n1 == 5)
        //            for(k = 0; k < n2; k++)
        //                dft5(t0+5*k, t1+5*k);
        //
        //        if(n1 == 4)
        //            for(k = 0; k < n2; k++)
        //                dft4(t0+4*k, t1+4*k);
        //
        //        if(n1 == 3)
        //            for(k = 0; k < n2; k++)
        //                dft3(t0+3*k, t1+3*k);
        //
        //        if(n1 == 2)
        //            for(k = 0; k < n2; k++)
        //                dft2(t0+2*k, t1+2*k);

        if n2 > 1 {
            //
            //            for(k =0; k < n; k++)
            //            {
            //                RE(t0[k]) = CMRE(t1[k], pw[k]);
            //                IM(t0[k]) = CMIM(t1[k], pw[k]);
            //            }

            //            matrix_transpose_cmplx(t0, n1, n2, t1);

            //            for(k = 0; k < n1; k++)
            //            {
            //                fft_krn(t1+k*n2, t0+k*n2, p, n2, addr+n);
            //            }
            //
            //            matrix_transpose_cmplx(t0, n2, n1, t1);
        }
        return t1;
    }
    //    return RES_OK;
}

#[cfg(test)]
mod tests {

    use super::*;

    #[ignore]
    #[test]
    fn fft_works() {
        let num = 14;

        let xs = stub(num);

        let fft = FFT::new(num).unwrap();
        let ys = fft.fft(&xs, num);

        let expects = expects();
        for i in 0..num {
            assert!((expects[i].0 - ys[i].re()).abs() < 1E-2);
            assert!((expects[i].1 - ys[i].im()).abs() < 1E-2);
        }
    }

    fn stub(num: usize) -> Vec<f64> {
        (0..num).map(|x| x as f64).collect()
    }

    fn expects() -> Vec<(f64, f64)> {
        vec![
            (91.000, 0.000),
            (-7.000, 30.669),
            (-7.000, 14.536),
            (-7.000, 8.778),
            (-7.000, 5.582),
            (-7.000, 3.371),
            (-7.000, 1.598),
            (-7.000, 0.000),
            (-7.000, -1.598),
            (-7.000, -3.371),
            (-7.000, -5.582),
            (-7.000, -8.778),
            (-7.000, -14.536),
            (-7.000, -30.669),
        ]
    }
}
