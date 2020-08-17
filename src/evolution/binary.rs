//! The `binary` module contains aliases and helper functions for working with binary
//! substrate data.
extern crate bitvec;

use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use super::gene::CrossOver;
use super::helper::{a_or_b, do_a_or_b};
use std::cell::RefCell;

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

impl CrossOver for BinarySubstrate {
    fn is_similar(&self, other: &Self) -> bool {
        self.len() == other.len()
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            // If the two substates are similar, randomly select bits of each parent.
            let recombined = RefCell::new(BitVec::new());
            for i in 0..self.len() {
                // A RefCell is used since only one of the functions will ever be executed at the
                // time.
                do_a_or_b(|| recombined.borrow_mut().push(*self.get(i).unwrap()),
                    || recombined.borrow_mut().push(*other.get(i).unwrap()));
            }
            recombined.into_inner().into()
        } else {
            // If the two substrates are not similar return a random one.
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

// This implementation improves the usability of the `CrossOver` trait for genomic elements.
impl CrossOver for usize {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        a_or_b(*self, *other)
    }
}

// This implementation improves the usability of the `CrossOver` trait for genomic elements.
impl<T> CrossOver for Option<T> where T: CrossOver + Sized + Clone {
    fn is_similar(&self, other: &Self) -> bool {
        match (self, other) {
            (None, None) => true,
            (Some(a), Some(b)) => a.is_similar(b),
            _ => false,
        }
    }

    fn cross_over(&self, other: &Self) -> Self {
        match (self, other) {
            (None, None) => None,
            (Some(a), Some(b)) => Some(a.cross_over(b)),
            _ => do_a_or_b(|| self.clone(), || other.clone()),
        }

    }
}

#[cfg(test)]
mod tests;
