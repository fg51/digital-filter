use std::f64::consts::{E, PI};

use num_complex::Complex64 as Complex;

pub fn fft(xs: &[f64]) -> Vec<Complex> {
    //    #f:サイズNの入力データ
    let num = xs.len();

    if num == 1 {
        return vec![Complex::new(xs[0], 0.)];
    }

    //let f_even = f[0:N:2]    #fの偶数番目の要素
    let xs_even = {
        let mut ts = vec![];
        let mut i = 0;
        while i < num {
            ts.push(xs[i]);
            i += 2;
        }
        ts
    };
    //    f_odd = f[1:N:2]    #fの奇数番目の要素
    let xs_odd = {
        let mut ts = vec![];
        let mut i = 1;
        while i < num {
            ts.push(xs[i]);
            i += 2;
        }
        ts
    };
    let ys_even = fft(&xs_even); //  #(3)偶数番目の要素でFFT
    let ys_odd = fft(&xs_odd); //    #(4)偶数番目の要素でFFT

    //tが0~N/2-1番目までのWを計算した配列
    //let w_n = np.exp(-1j * (2 * PI * np.arange(0, num / 2)) / N)
    //let w_n = (-1j * (2. * PI * np.arange(0, num / 2)) / num).exp();
    let w_n: Vec<Complex> = (0..(num / 2))
        .map(|x| Complex::new(0., -1. * (2. * PI * (x as f64)) / (num as f64)).expf(E))
        .collect();

    //    F = np.zeros(N, dtype ='complex')    #FFTの出力
    //    F[0:N//2] = F_even + W_N * F_odd    #(9)を計算(t:0~N/2-1)
    let mut zs: Vec<Complex> = (0..(num / 2))
        .map(|i| ys_even[i] + w_n[i] * ys_odd[i])
        .collect();
    //F[N//2:N] = F_even - W_N * F_odd    #(10)を計算(t:N/2~N-1)
    let zs1: Vec<Complex> = ((num / 2)..num)
        .map(|i| ys_even[i] - w_n[i] * ys_odd[i])
        .collect();
    zs.extend(zs1);
    return zs;
}
