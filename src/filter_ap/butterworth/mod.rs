use crate::values::Complex;

type Result<T> = std::result::Result<T, String>;

mod butter_ap_zp;
use butter_ap_zp::butter_ap_zp;

use crate::filter::filter_zp2ab;

//int DSPL_API butter_ap(double rp, int ord, double* b, double* a)
pub fn butter_ap(ripple: f64, ord: usize, b: &mut [f64], a: &mut [f64]) -> Result<()> {
    //    int res;
    //    complex_t *z = NULL;
    //    complex_t *p = NULL;

    //    if(rp < 0.0)
    //        return ERROR_FILTER_RP;
    //    if(ord < 1)
    //        return ERROR_FILTER_ORD;
    //    if(!a || !b)
    //        return ERROR_PTR;

    let mut zeros = vec![Complex::new(0., 0.); ord];
    let mut poles = vec![Complex::new(0., 0.); ord];

    let (num_of_zeros, num_of_poles) = butter_ap_zp(ord, ripple, &mut zeros, &mut poles)?;
    //    if(res != RES_OK)
    //        goto exit_label;

    filter_zp2ab(&zeros, num_of_zeros, &poles, num_of_poles, ord, b, a)?;
    //    if(res != RES_OK)
    //        goto exit_label;

    b[0] = a[0];

    //exit_label:
    //    if(z)
    //        free(z);
    //    if(p)
    //        free(p);
    //    return res;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn butter_ap_works() {
        //#include <stdio.h>
        //#include <stdlib.h>
        //#include <string.h>
        //#include "dspl.h"

        // Filter order
        let ord = 3;

        // Frequency response vector size
        //let num = 1000;

        //int main(int argc, char* argv[])
        //{
        //    void* hdspl;    /* DSPL handle */
        //    void* hplot;    /* GNUPLOT handle */
        //
        //    /* Load DSPL functions */
        //    hdspl = dspl_load();

        let mut a = vec![0.0; ord + 1]; // H(s) numerator     coefficients vector
        let mut b = vec![0.0; ord + 1]; // H(s) denominator coefficients vector
        let ripple = 1.0; // Magnitude ripple from 0 to 1 rad/s
                          //    double w[N];     /* Angular frequency (rad/s)               */
                          //    double mag[N];   /* Filter Magnitude (dB)                   */
                          //    double phi[N];   /* Phase response                          */
                          //    double tau[N];   /* Group delay                             */
                          //    int k;

        // H(s) coefficients calculation
        //    int res = butter_ap(Rp, ORD, b, a);
        butter_ap(ripple, ord, &mut b, &mut a).unwrap();

        //b[ 0] =     1.965     a[ 0] =     1.965
        //b[ 1] =     0.000     a[ 1] =     3.138
        //b[ 2] =     0.000     a[ 2] =     2.505
        //b[ 3] =     0.000     a[ 3] =     1.000
        for (i, expect) in [1.965, 0.0, 0.0, 0.0].iter().enumerate() {
            assert!(
                (b[i] - expect).abs() < 1E-3,
                "b[{}]: except: {}, actual: {}",
                i,
                expect,
                b[i]
            );
        }
        for (i, expect) in [1.965, 3.138, 2.505, 1.000].iter().enumerate() {
            assert!(
                (a[i] - expect).abs() < 1E-3,
                "a[{}]: except: {}, actual: {}",
                i,
                expect,
                a[i]
            );
        }

        //    /* Frequency in logarithmic scale from 0.01 to 100 rad/s */
        //    logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);
        //
        //    /* Filter frequency parameter calculation */
        //    filter_freq_resp(b, a, ORD, w, N,
        //                     DSPL_FLAG_LOGMAG|DSPL_FLAG_UNWRAP|DSPL_FLAG_ANALOG,
        //                     mag, phi, tau);
        //
        //    /* Write Magnitude, phase response and group delay to the files */
        //    writetxt(w, mag, N, "dat/butter_ap_test_mag.txt");
        //    writetxt(w, phi, N, "dat/butter_ap_test_phi.txt");
        //    writetxt(w, tau, N, "dat/butter_ap_test_tau.txt");
    }
}
