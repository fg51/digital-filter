pub fn linspace(start: f64, stop: f64, num: usize, endpoint: Option<bool>) -> Vec<f64> {
    let div = match endpoint {
        Some(true) | None => (num - 1) as f64,
        Some(false) => num as f64,
    };

    let delta = stop - start;
    //let y = _nx.arange(0, num, dtype=dt).reshape((-1,) + (1,) * ndim(delta))
    //let y: Vec<f64> = (0..num).map(|i| (i as f64) * delta).collect();
    //# In-place multiplication y *= delta/div is faster, but prevents the multiplicant
    //# from overriding what class is produced, and thus prevents, e.g. use of Quantities,
    //# see gh-7142. Hence, we multiply in place only for standard scalar types.
    //_mult_inplace = _nx.isscalar(delta)
    //let y = if div > 0. {
    //    let step = delta / div;
    //    //if _nx.any(step == 0):
    //    //    // Special handling for denormal numbers, gh-5437
    //    //    y /= div
    //    //    if _mult_inplace:
    //    //        y *= delta
    //    //    else:
    //    //        y = y * delta
    //    //else:
    //    //    if _mult_inplace:
    //    //        y *= step
    //    //    else:
    //    //        y = y * step
    //    todo!();
    //} else {
    //    //# sequences with 0 items or 1 item with endpoint=True (i.e. div <= 0)
    //    //# have an undefined step
    //    //step = NaN
    //    //# Multiply with delta to allow possible override of output class.
    //    //y = y * delta
    //    todo!();
    //};

    let step = delta / div;
    let mut y: Vec<f64> = (0..num).map(|i| start + step * (i as f64)).collect();

    if let Some(endpoint) = endpoint {
        if endpoint == true {
            if num > 1 {
                y[num - 1] = stop;
            }
        }
    }

    //if axis != 0:
    //    y = _nx.moveaxis(y, 0, axis)

    //if _nx.issubdtype(dtype, _nx.integer):
    //    _nx.floor(y, out=y)

    //if retstep:
    //    return y.astype(dtype, copy=False), step
    //else:
    //    return y.astype(dtype, copy=False)
    return y;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linspace_works1() {
        let expect = [0., 5., 10.];
        let actual = linspace(0., 10., 3, None);

        for i in 0..3 {
            assert!(
                (actual[i] - expect[i]).abs() < 1E-4,
                "index: {}, actual: {},expect: {}",
                i,
                actual[i],
                expect[i],
            );
        }
    }

    #[test]
    fn linspace_works2() {
        let expect = [0., 3.33333, 6.6666, 10.];
        let actual = linspace(0., 10., 4, None);

        for i in 0..4 {
            assert!(
                (actual[i] - expect[i]).abs() < 1E-4,
                "index: {}, actual: {},expect: {}",
                i,
                actual[i],
                expect[i],
            );
        }
    }
}
