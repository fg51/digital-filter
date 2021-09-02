pub type ComplexF64 = Complex;

#[derive(Clone, Copy)]
pub struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub const fn zero() -> Self {
        Self { re: 0., im: 0. }
    }

    pub const fn re(&self) -> f64 {
        self.re
    }

    pub const fn im(&self) -> f64 {
        self.im
    }

    pub fn sub(&mut self, rh: Complex) {
        self.re -= rh.re;
        self.im -= rh.im;
    }
}

impl std::ops::Add for Complex {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl std::ops::AddAssign for Complex {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            re: self.re + other.re,
            im: self.im + other.im,
        };
    }
}

impl std::ops::Sub for Complex {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

impl std::ops::SubAssign for Complex {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            re: self.re - other.re,
            im: self.im - other.im,
        };
    }
}
