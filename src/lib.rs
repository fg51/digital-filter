pub mod errors;
mod service;
pub mod values;

pub mod complex;
pub mod dft;
pub mod fft;
pub mod ffti;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
