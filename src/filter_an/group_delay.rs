use crate::errors::Result;

use super::DSPL_FLAG_ANALOG;

//int DSPL_API group_delay(double* pb, double* pa, int ord, int flag,
//                         double* w, int n, double* tau)
pub fn group_delay(
    pb: &[f64],
    pa: Option<&[f64]>,
    order: usize,
    flag: u32,
    w: &[f64],
    num: usize,
    tau: &mut [f64],
) -> Result<()> {
    //    double a, b, c, d, da, db, dc, dd, f, e;
    //    int t, m;
    //
    //    double *qa = NULL;
    //
    //    if(!pb || !w || !tau || (!pa && (flag & DSPL_FLAG_ANALOG)))
    //        return ERROR_PTR;
    //    if(ord < 1)
    //        return ERROR_FILTER_ORD;
    //    if(n < 1)
    //        return ERROR_SIZE;

    let mut qa1 = vec![0.; order + 1];
    qa1[0] = 1.0;

    let qa = match pa {
        Some(pa) => pa,
        None => &qa1,
    };

    for t in 0..num {
        let (mut a, mut b, mut c, mut d) = (0., 0., 0., 0.);
        let (mut da, mut db, mut dc, mut dd) = (0., 0., 0., 0.);
        if (flag & DSPL_FLAG_ANALOG) > 0 {
            // for(m = 0; m < ord+1; m+=4) {
            let mut m = 0;
            while m < order + 1 {
                a += pb[m] * (w[t].powf(m as f64));
                c += qa[m] * (w[t].powf(m as f64));
                da += pb[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                dc += qa[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                m += 4;
            }
            // for(m = 2; m < ord+1; m+=4) {
            let mut m = 2;
            while m < order + 1 {
                a -= pb[m] * (w[t].powf(m as f64));
                c -= qa[m] * (w[t].powf(m as f64));
                da -= pb[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                dc -= qa[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                m += 4;
            }

            // for(m = 1; m < ord+1; m+=4) {
            let mut m = 1;
            while m < order + 1 {
                b += pb[m] * (w[t].powf(m as f64));
                d += qa[m] * (w[t].powf(m as f64));
                db += pb[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                dd += qa[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                m += 4;
            }

            // for(m = 3; m < ord+1; m+=4) {
            let mut m = 3;
            while m < order + 1 {
                b -= pb[m] * (w[t].powf(m as f64));
                d -= qa[m] * (w[t].powf(m as f64));
                db -= pb[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                dd -= qa[m] * (m as f64) * (w[t].powf((m - 1) as f64));
                m += 4;
            }
        } else {
            for m in 0..order + 1 {
                let (sinw, cosw) = ((w[t] * (m as f64)).sin(), (w[t] * (m as f64)).cos());
                a += pb[m] * cosw;
                b -= pb[m] * sinw;
                c += qa[m] * cosw;
                d -= qa[m] * sinw;

                da -= pb[m] * (m as f64) * sinw;
                db -= pb[m] * (m as f64) * cosw;
                dc -= qa[m] * (m as f64) * sinw;
                dd -= qa[m] * (m as f64) * cosw;
            }
        }

        let f = da * c + a * dc + db * d + b * dd;
        let e = db * c + b * dc - da * d - a * dd;
        tau[t] = (f * (b * c - a * d) - e * (a * c + b * d)) / ((a * a + b * b) * (c * c + d * d));
    }
    //
    //    if(qa != pa)
    //      free(qa);
    //
    //    return RES_OK;
    todo!();
}
