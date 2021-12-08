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

pub mod filter_design;
pub mod fir_filter_design;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
