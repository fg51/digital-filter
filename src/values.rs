pub use std::f64::consts::PI;

pub use num_complex::Complex64 as Complex;

pub enum Direction {
    Forward, // kernel uses "-1" sign
    Inverse, // kernel uses "+1" sign
}
