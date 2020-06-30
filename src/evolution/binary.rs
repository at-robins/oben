//! The `binary` module contains aliases and helper functions for working with binary
//! substrate data.
extern crate bitvec;

use bitvec::{boxed::BitBox, order::Msb0};

pub type BinarySubstrate = BitBox<Msb0, u8>;
