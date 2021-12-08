use std::f64::consts::PI;

use crate::complex::Complex;

pub type Result<T> = std::result::Result<T, String>;

use super::FFT;

const FFT_COMPOSITE_MAX: usize = 46340; // sqrt(1 << 31)

impl FFT {
    //int DSPL_API fft_create(fft_t* pfft, int n)
    pub fn new(n: usize) -> Result<Self> {
        //    int n1, n2, addr, s, k, m, nw, err;
        let mut fft = FFT::zero();
        let mut nw = 0;
        let mut addr = 0;

        let mut s = n;
        while s > 1 {
            let n2 = cal_n2(s);

            if n2 == 1 {
                if s > FFT_COMPOSITE_MAX {
                    return Err("ERROR FFT SIZE".to_string());
                }

                nw += s;
                //pfft->w = if pfft->w {
                //          (complex_t*) realloc(pfft->w,  nw*sizeof(complex_t)) } else {
                //          (complex_t*) malloc(           nw*sizeof(complex_t))};
                fft.w = vec![Complex::zero(); nw];
                for k in 0..s {
                    let phi = -2. * PI * k as f64 / s as f64;
                    fft.w[addr] = Complex::new(phi.cos(), phi.sin());
                    addr += 1;
                }
                s = 1;
            } else {
                let n1 = s / n2;
                nw += s;
                //pfft->w = pfft->w ?
                //          (complex_t*) realloc(pfft->w,    nw*sizeof(complex_t)):
                //          (complex_t*) malloc(             nw*sizeof(complex_t));
                fft.w = vec![Complex::zero(); nw];

                for k in 0..n1 {
                    for m in 0..n2 {
                        let phi = -2. * PI * (k * m) as f64 / s as f64;
                        fft.w[addr] = Complex::new(phi.cos(), phi.sin());
                        addr += 1;
                    }
                }
            }
            s /= n2;
        }

        //    pfft->t0 = pfft->t0 ? (complex_t*) realloc(pfft->t0, n*sizeof(complex_t)):
        //                          (complex_t*) malloc(           n*sizeof(complex_t));
        fft.t0 = vec![Complex::zero(); n];
        //
        //    pfft->t1 = pfft->t1 ? (complex_t*) realloc(pfft->t1, n*sizeof(complex_t)):
        //                          (complex_t*) malloc(           n*sizeof(complex_t));
        fft.t1 = vec![Complex::zero(); n];
        fft.n = n;

        fill_w32(&mut fft);
        fill_w64(&mut fft);
        fill_w128(&mut fft);
        fill_w256(&mut fft);
        fill_w512(&mut fft);
        fill_w1024(&mut fft);
        fill_w2048(&mut fft);
        fill_w4096(&mut fft);

        //error_proc:
        //    if(pfft->t0) free(pfft->t0);
        //    if(pfft->t1) free(pfft->t1);
        //    if(pfft->w)    free(pfft->w);
        //    pfft->n = 0;
        //    return err;
        return Ok(fft);
    }
}

pub fn cal_n2(s: usize) -> usize {
    for i in [
        4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 7, 8, 5, 4, 3, 2,
    ] {
        if s % i == 0 {
            return i;
        }
    }
    1
}

fn fill_w32(x: &mut FFT) {
    let (m, n) = (4, 8);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w32[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w64(x: &mut FFT) {
    let (m, n) = (8, 8);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w64[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w128(x: &mut FFT) {
    let (m, n) = (8, 16);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w128[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w256(x: &mut FFT) {
    let (m, n) = (16, 16);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w256[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w512(x: &mut FFT) {
    let (m, n) = (16, 32);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w512[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w1024(x: &mut FFT) {
    let (m, n) = (32, 32);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w1024[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w2048(x: &mut FFT) {
    let (m, n) = (32, 64);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w2048[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}

fn fill_w4096(x: &mut FFT) {
    let (m, n) = (16, 256);
    let mn = (m * n) as f64;
    let mut addr = 0;
    for i in 0..m {
        for j in 0..n {
            let phi = -2. * PI * (i * j) as f64 / mn;
            x.w4096[addr] = Complex::new(phi.cos(), phi.sin());
            addr += 1;
        }
    }
}
