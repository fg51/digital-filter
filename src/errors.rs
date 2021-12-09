use thiserror::Error;

//pub type Result<T> = std::result::Result<T, String>;
pub type Result<T> = std::result::Result<T, ErrorKind>;

#[derive(Debug, Error)]
pub enum ErrorKind {
    #[error("Value Error: {}", 0)]
    ValueError(String),
}
