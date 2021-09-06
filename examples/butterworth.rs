use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::Result;

use digital_filter as lib;

pub fn main() -> Result<()> {
    let order = 3;

    let mut a = vec![0.0; order + 1]; // H(s) numerator     coefficients vector
    let mut b = vec![0.0; order + 1]; // H(s) denominator coefficients vector
    let ripple = 1.0; // Magnitude ripple from 0 to 1 rad/s
                      //    double w[N];     /* Angular frequency (rad/s)               */
                      //    double mag[N];   /* Filter Magnitude (dB)                   */
                      //    double phi[N];   /* Phase response                          */
                      //    double tau[N];   /* Group delay                             */
                      //    int k;

    // H(s) coefficients calculation
    //    int res = butter_ap(Rp, ORD, b, a);
    lib::filter_ap::butterworth::butter_ap(ripple, order, &mut b, &mut a).unwrap();

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

    // Frequency response vector size
    let num = 1000;

    //    /* Frequency in logarithmic scale from 0.01 to 100 rad/s */
    let w = logspace(-2.0, 2.0, num, DSPL_SYMMETRIC);

    //    /* Filter frequency parameter calculation */
    //lib::filter_an::freq_response::FreqResponse
    //filter_freq_resp(
    //    b,
    //    a,
    //    ORD,
    //    w,
    //    N,
    //    DSPL_FLAG_LOGMAG | DSPL_FLAG_UNWRAP | DSPL_FLAG_ANALOG,
    //    mag,
    //    phi,
    //    tau,
    //);
    let res = lib::filter_an::freq_response::FreqResponse::analog(&b, &a, order, w, num).unwrap();
    let mag = res.mag(true);

    //    /* Write Magnitude, phase response and group delay to the files */
    let writer = BufWriter::new(File::create("dat/butter_ap_test_mag.txt")?);
    for i in 0..num {
        writeln!(writer, "{},{}", w[i], mag[i],)?;
    }
    //    writetxt(w, phi, N, "dat/butter_ap_test_phi.txt");
    //    writetxt(w, tau, N, "dat/butter_ap_test_tau.txt");
    Ok(())
}
