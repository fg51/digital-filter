use anyhow::Result;

use digital_filter as lib;

use lib::filter_design::{
    iirfilter::iirfilter,
    values::{FilterForm, FilterKind, IIRFilterKind},
};

pub fn main() -> Result<()> {
    let fs = 48000.;
    let bpfc1 = 100.;
    let order = 2;
    let rs = 40.;
    let rp = 1.0;
    let analog = false;
    let ftype = IIRFilterKind::Chebyshev1;
    let output = FilterForm::Ba;
    let (b, a) = iirfilter(
        order,
        &[bpfc1],
        Some(rp),
        Some(rs),
        Some(FilterKind::LowPass),
        Some(analog),
        Some(ftype),
        Some(output),
        Some(fs),
    )?;

    println!("b");
    for i in b {
        println!("{:e}", i);
    }
    println!("a");
    for i in a {
        println!("{:e}", i);
    }
    Ok(())
}
