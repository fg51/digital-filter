use crate::errors::Result;
use crate::values::Complex;
use crate::values::PI;

use super::freqs;
use super::freqz;
use super::phase_unwrap;

const DSPL_FLAG_DIGITAL: u32 = 0x00000000;
const DSPL_FLAG_ANALOG: u32 = 0x00000001;
const DSPL_FLAG_LOGMAG: u32 = 0x00000002;
const DSPL_FLAG_UNWRAP: u32 = 0x00000004;
const DSPL_FLAG_FFT_SHIFT: u32 = 0x00000008;
const DSPL_FLAG_PSD_TWOSIDED: u32 = DSPL_FLAG_FFT_SHIFT;

//int DSPL_API filter_freq_resp(double* b, double* a, int ord, double* w, int n, int flag, double* mag, double* phi, double* tau)
//pub fn freq_response(
//    b: &[f64],
//    a: Option<&[f64]>,
//    order: usize,
//    w: &[f64],
//    n: usize,
//    flag: u32,
//    mag: &mut [f64],
//    phi: &mut [f64],
//    _tau: &mut [f64],
//) -> Result<()> {
//    //    int res, k, flag_analog;
//
//    //    complex_t *hc  = NULL;
//    //    double *phi0   = NULL;
//    //    double *phi1   = NULL;
//    //    double *w0     = NULL;
//    //    double *w1     = NULL;
//
//    //    if(!b || !w)
//    //        return ERROR_PTR;
//    //    if(ord < 1)
//    //        return ERROR_FILTER_ORD;
//    //    if(n < 1)
//    //        return ERROR_SIZE;
//
//    let is_analog = (flag & DSPL_FLAG_ANALOG) > 0;
//
//    //    hc = (complex_t*) malloc (n*sizeof(complex_t));
//
//    let mut hc = vec![Complex::new(0., 0.); n];
//    if is_analog {
//        freqs(b, a.unwrap(), order, w, n, &mut hc)?;
//    } else {
//        freqz(b, a, order, w, n, &mut hc)?;
//    }
//
//    //    if(res != RES_OK)
//    //        goto exit_label;
//
//    //    if(mag)
//    //    {
//    //        if(flag & DSPL_FLAG_LOGMAG)
//    //        {
//    for i in 0..n {
//        mag[i] = 10.0 * (hc[i].norm_sqr().log10());
//    }
//
//    //        } else {
//    //            for(k = 0; k < n; k++)
//    //                mag[k] = sqrt(ABSSQR(hc[k]));
//    //        }
//    //    }
//
//    //    if(phi)
//    //    {
//    //        for(k = 0; k < n; k++)
//    //            phi[k] = atan2(IM(hc[k]), RE(hc[k]));
//    //
//    //        if(flag & DSPL_FLAG_UNWRAP)
//    //        {
//    //            res = unwrap(phi, n, M_2PI, 0.8);
//    //            if(res != RES_OK)
//    //                goto exit_label;
//    //        }
//    //    }
//
//    //    if(tau)
//    //        res = group_delay(b, a, ord, flag, w, n, tau);
//
//    //exit_label:
//    //    if(hc)
//    //        free(hc);
//    //    if(phi0)
//    //        free(phi0);
//    //    if(phi1)
//    //        free(phi1);
//    //    if(w0)
//    //        free(w0);
//    //    if(w1)
//    //        free(w1);
//    //    return res;
//    return Ok(());
//}

pub struct FreqResponse {
    hc: Vec<Complex>,
    num: usize,
}

impl FreqResponse {
    pub fn analog(b: &[f64], a: &[f64], order: usize, w: &[f64], num: usize) -> Result<Self> {
        let mut hc = vec![Complex::new(0., 0.); num];
        freqs(b, a, order, w, num, &mut hc)?;
        Ok(Self { hc, num })
    }

    pub fn fir(b: &[f64], order: usize, w: &[f64], num: usize) -> Result<Self> {
        let mut hc = vec![Complex::new(0., 0.); num];
        freqz(b, None, order, w, num, &mut hc)?;
        Ok(Self { hc, num })
    }

    pub fn iir(b: &[f64], a: &[f64], order: usize, w: &[f64], num: usize) -> Result<Self> {
        let mut hc = vec![Complex::new(0., 0.); num];
        freqz(b, Some(a), order, w, num, &mut hc)?;
        Ok(Self { hc, num })
    }

    pub fn mag(&self, is_log: bool) -> Vec<f64> {
        if is_log == true {
            (0..self.num)
                .map(|i| 10.0 * self.hc[i].norm_sqr().log10())
                .collect()
        } else {
            // mag[k] = sqrt(ABSSQR(hc[k]));
            (0..self.num)
                .map(|i| self.hc[i].norm_sqr().sqrt())
                .collect()
        }
    }

    pub fn phi(&self, is_unwrap: bool) -> Vec<f64> {
        let mut phi: Vec<f64> = vec![];
        for i in 0..self.num {
            //phi[i] = atan2(self.hc[i].im, self.hc[i].re);
            phi[i] = self.hc[i].im.atan2(self.hc[i].re);
        }

        if is_unwrap == true {
            phase_unwrap(&mut phi, self.num, 2. * PI, 0.8);
            //  if(res != RES_OK)
            //  goto exit_label;
        }
        return phi;
    }

    pub fn tau(&self) -> Vec<f64> {
        //        res = group_delay(b, a, ord, flag, w, n, tau);
        todo!();
    }
}
