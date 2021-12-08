use std::f64::consts::PI;

use crate::complex::ComplexF64 as Complex;

pub type Result<T> = std::result::Result<T, String>;

pub fn dft(xs: &[f64], num: usize) -> Result<Vec<Complex>> {
    if num < 1 {
        return Err("ERROR SIZE".to_string());
    }

    let step = 1.0 / num as f64;

    let mut ys = vec![Complex::zero(); num];
    for i in 0..num {
        for j in 0..num {
            let phi = -2. * PI * step * i as f64 * j as f64;
            ys[i] += Complex::new(xs[j] * phi.cos(), xs[j] * phi.sin());
        }
    }
    return Ok(ys);
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn dft_works() {
        let num = 16;
        let xs = stub(num);
        let ys = dft(&xs, num).unwrap();

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
            (120., 0.),
            (-8.000, 40.219),
            (-8.000, 19.314),
            (-8.000, 11.973),
            (-8.000, 8.000),
            (-8.000, 5.345),
            (-8.000, 3.314),
            (-8.000, 1.591),
            (-8.000, 0.000),
            (-8.000, -1.591),
            (-8.000, -3.314),
            (-8.000, -5.345),
            (-8.000, -8.000),
            (-8.000, -11.973),
            (-8.000, -19.314),
            (-8.000, -40.219),
        ]
    }
}
