pub mod errors;
mod service;
pub mod values;

pub mod complex;
pub mod conv;
pub mod polyval;

pub mod dft;
pub mod fft;
pub mod ffti;

pub mod butter;
mod filter;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
