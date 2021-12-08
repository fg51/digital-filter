pub type Result<T> = std::result::Result<T, String>;

pub fn phase_unwrap(phi: &mut [f64], num: usize, lev: f64, mar: f64) -> Result<()> {
    //    double a[2] = {0.0, 0.0};
    //    double d;
    //    double th;
    //    int k;
    //    int flag = 1;

    //    if(!phi)
    //        return ERROR_PTR;
    //
    //    if(n<1)
    //        return ERROR_SIZE;
    //
    //    if(lev<=0 || mar <=0)
    //        return ERROR_UNWRAP;

    let threshold = mar * lev;
    let mut flag = false;
    while flag == false {
        flag = true;
        let mut a = [0.0, 0.0];
        for i in 0..num - 1 {
            let d = phi[i + 1] - phi[i];
            if d > threshold {
                a[0] -= lev;
                flag = false;
            }
            if d < -threshold {
                a[0] += lev;
                flag = false;
            }
            phi[i] += a[1];
            a[1] = a[0];
        }
        phi[num - 1] += a[1];
    }

    return Ok(());
}
