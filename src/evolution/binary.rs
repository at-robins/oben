//! The `binary` module contains aliases and helper functions for working with binary
//! substrate data.
extern crate bitvec;

use bitvec::{boxed::BitBox, order::Msb0};

/// A type alias for the underlying binary representation of [`Substrate`]s.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
pub type BinarySubstrate = BitBox<Msb0, u8>;

/// Converts the [`Substrate`] into a 64 bit array.
///
/// # Parameters
///
/// * `substrate` - the [`Substrate`] to convert
///
/// # Panics
///
/// If the [`Substrate`] is not of length `64`.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
pub fn as_64(substrate: &BinarySubstrate) -> [u8; 8] {
    if substrate.len() != 64 {
        panic!("64 bits are required for conversion, but {} were passed.", substrate.len());
    }
    let mut array = [0; 8];
    let bytes = &substrate.as_slice()[..array.len()];
    array.copy_from_slice(bytes);
    array
}

/// Converts the [`Substrate`] into a 64 bit float.
///
/// # Parameters
///
/// * `substrate` - the [`Substrate`] to convert
///
/// # Panics
///
/// If the [`Substrate`] is not of length `64`.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
pub fn as_f64(substrate: &BinarySubstrate) -> f64 {
    f64::from_be_bytes(as_64(substrate))
}
