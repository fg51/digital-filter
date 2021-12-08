#![allow(dead_code)]
pub mod freq_response;

mod freqs;
pub use freqs::freqs;

pub mod freqz;
pub use freqz::freqz;

pub mod phase_unwrap;
pub use phase_unwrap::phase_unwrap;

pub mod group_delay;
pub use group_delay::group_delay;

const DSPL_FLAG_DIGITAL: u32 = 0x00000000;
const DSPL_FLAG_ANALOG: u32 = 0x00000001;
const DSPL_FLAG_LOGMAG: u32 = 0b1 << 1; //0x00000002;
const DSPL_FLAG_UNWRAP: u32 = 0b1 << 2; // 0x00000004;
const DSPL_FLAG_FFT_SHIFT: u32 = 0b1 << 3; // 0x00000008;
const DSPL_FLAG_PSD_TWOSIDED: u32 = DSPL_FLAG_FFT_SHIFT;
