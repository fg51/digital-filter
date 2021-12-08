//def _get_fs(fs, nyq):
// Utility for replacing the argument 'nyq' (with default 1) with 'fs'.
pub fn get_fs(fs: Option<f64>, nyq: Option<f64>) -> f64 {
    match nyq {
        Some(nyq) => 2. * nyq,
        None => match fs {
            Some(fs) => fs,
            None => 2.,
        },
    }
}

pub mod signaltools;

mod linspace;
pub use linspace::linspace;
