//! The `binary` module contains aliases and helper functions for working with binary
//! substrate data.
extern crate bitvec;

pub use binary_chemistry::{BinaryReaction, BinaryState};
use rand::{thread_rng, Rng};
//pub use binary_mutation::BinaryMutation;

use super::chemistry::Information;
use super::gene::CrossOver;
use super::helper::do_a_or_b;
use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use std::cell::RefCell;

/// A type alias for the underlying binary representation of [`Substrate`]s.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
pub type BinarySubstrate = BitBox<u8, Msb0>;
/*
/// A type alias for the representation of a binary [`Genome`].
///
/// [`Genome`]: ../gene/struct.Genome.html
pub type BinaryGenome = Genome<BinaryReaction, BinaryState, BinarySubstrate>;
*/

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
    let bytes = &substrate.as_raw_slice()[..array.len()];
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

/// Converts a 64 bit float into its binary reprsentation
/// as a [`BinarySubstrate`] that allows cross-over.
///
/// # Parameters
///
/// * `value` - the number to convert to its binary representation.
pub fn f64_to_binary(value: f64) -> BinarySubstrate {
    BitBox::from_boxed_slice(Box::new(value.to_be_bytes()))
}

/// Converts the [`Substrate`] into a 64 bit unsigned integer.
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
pub fn as_u64(substrate: &BinarySubstrate) -> u64 {
    u64::from_be_bytes(as_64(substrate))
}

/// Converts a 64 bit unsigned integer into its binary reprsentation
/// as a [`BinarySubstrate`] that allows cross-over.
///
/// # Parameters
///
/// * `value` - the number to convert to its binary representation.
pub fn u64_to_binary(value: u64) -> BinarySubstrate {
    BitBox::from_boxed_slice(Box::new(value.to_be_bytes()))
}

/// Returns a copy of the [`BinarySubstrate`] with a random bit flipped.
///
/// # Parameters
///
/// * `base` - the substrate to mutate
pub fn flip_random_bit(base: &BinarySubstrate) -> BinarySubstrate {
    let mut binary_base = base.clone();
    if binary_base.len() > 0 {
        let random_bit_index = thread_rng().gen_range(0..binary_base.len());
        // The unwrap should always work, since the index being in range was checked for.
        // The clone should be cheap as a bool primitive is cloned.
        let random_bit = *binary_base.get(random_bit_index).unwrap();
        // Flip the bit.
        binary_base.set(random_bit_index, random_bit ^ true);
    }
    binary_base
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
                do_a_or_b(
                    || recombined.borrow_mut().push(*self.get(i).unwrap()),
                    || recombined.borrow_mut().push(*other.get(i).unwrap()),
                );
            }
            recombined.into_inner().into()
        } else {
            // If the two substrates are not similar return a random one.
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

impl Information for BinarySubstrate {
    fn update_value(&mut self, _time_passed: i32) {
        // Does nothing as no update of the internal value is required.
    }
}

mod binary_chemistry;
//mod binary_mutation;
#[cfg(test)]
mod tests;
