pub mod errors;
mod service;
pub mod values;

pub mod complex;
pub mod conv;
pub mod polyval;

pub mod dft;
pub mod fft;
pub mod ffti;

mod filter;
pub mod filter_an;
pub mod filter_ap;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
