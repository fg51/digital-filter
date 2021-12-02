use crate::errors::Result;
use crate::values::Complex;
use crate::values::PI;

use super::freqs;
use super::freqz;
use super::group_delay;
use super::phase_unwrap;

//use super:: DSPL_FLAG_DIGITAL;
use super::DSPL_FLAG_ANALOG;
use super::DSPL_FLAG_LOGMAG;
use super::DSPL_FLAG_UNWRAP;
//use super:: DSPL_FLAG_FFT_SHIFT;
//use super:: DSPL_FLAG_PSD_TWOSIDED;

//int DSPL_API filter_freq_resp(double* b, double* a, int ord, double* w, int n, int flag, double* mag, double* phi, double* tau)
pub fn freq_response(
    b: &[f64],
    a: Option<&[f64]>,
    order: usize,
    w: &[f64],
    n: usize,
    flag: u32,
    mag: Option<&mut [f64]>,
    phi: Option<&mut [f64]>,
    tau: Option<&mut [f64]>,
) -> Result<()> {
    //    //    int res, k, flag_analog;
    //
    //    //    complex_t *hc  = NULL;
    //    //    double *phi0   = NULL;
    //    //    double *phi1   = NULL;
    //    //    double *w0     = NULL;
    //    //    double *w1     = NULL;

    //    //    if(!b || !w)
    //    //        return ERROR_PTR;
    //    //    if(ord < 1)
    //    //        return ERROR_FILTER_ORD;
    //    //    if(n < 1)
    //    //        return ERROR_SIZE;

    let is_analog = (flag & DSPL_FLAG_ANALOG) > 0;

    //    //    hc = (complex_t*) malloc (n*sizeof(complex_t));
    //
    //    let mut hc = vec![Complex::new(0., 0.); n];
    let fr = if is_analog {
        //        freqs(b, a.unwrap(), order, w, n, &mut hc)?;
        FreqResponse::analog(b, a.unwrap(), order, w, n)?
    } else {
        match a {
            Some(a) => FreqResponse::iir(b, a, order, w, n)?,
            None => FreqResponse::fir(b, order, w, n)?,
        }
    };

    if let Some(mag) = mag {
        if (flag & DSPL_FLAG_LOGMAG) > 0 {
            for i in 0..n {
                mag[i] = 10.0 * (fr.hc[i].norm_sqr().log10());
            }
        } else {
            for _k in 0..n {
                //mag[k] = sqrt(ABSSQR(hc[k]));
                todo!();
            }
        }
    }

    if let Some(phi) = phi {
        for (i, p) in fr
            .phi((flag & DSPL_FLAG_UNWRAP) > 0)?
            .into_iter()
            .enumerate()
        {
            phi[i] = p;
        }
    }

    if let Some(tau) = tau {
        //        res = group_delay(b, a, ord, flag, w, n, tau);
        for (i, t) in fr.tau(b, a, order, flag, w, n)?.into_iter().enumerate() {
            tau[i] = t;
        }
    }

    //exit_label:
    //    if(hc)
    //        free(hc);
    //    if(phi0)
    //        free(phi0);
    //    if(phi1)
    //        free(phi1);
    //    if(w0)
    //        free(w0);
    //    if(w1)
    //        free(w1);
    //    return res;
    //  return Ok(());
    todo!();
}

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

    pub fn phi(&self, is_unwrap: bool) -> Result<Vec<f64>> {
        let mut phi: Vec<f64> = vec![];
        for i in 0..self.num {
            //phi[i] = atan2(self.hc[i].im, self.hc[i].re);
            phi[i] = self.hc[i].im.atan2(self.hc[i].re);
        }

        if is_unwrap == true {
            phase_unwrap(&mut phi, self.num, 2. * PI, 0.8)?;
            //  if(res != RES_OK)
            //  goto exit_label;
        }
        return Ok(phi);
    }

    pub fn tau(
        &self,
        b: &[f64],
        a: Option<&[f64]>,
        order: usize,
        flag: u32,
        w: &[f64],
        n: usize,
    ) -> Result<Vec<f64>> {
        //res = group_delay(b, a, ord, flag, w, n, tau);
        let mut tau = vec![];
        let _res = group_delay(b, a, order, flag, w, n, &mut tau)?;
        return Ok(tau);
    }
}
