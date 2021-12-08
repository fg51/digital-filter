use digital_filter as lib;

use lib::fir_filter_design::{firwin, values::WindowKind};

pub fn main() {
    let b = firwin(
        80,
        &[0.5],
        None,
        Some(WindowKind::Kaiser(8.)),
        None,
        None,
        None,
        None,
    )
    .unwrap();

    for i in b {
        println!("{}", i);
    }
}
