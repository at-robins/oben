//! The `helper` module contains helper constructs for general workflow.
extern crate bitvec;
extern crate rand;
extern crate serde;

use rand::{Rng, thread_rng};
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;

/// Randomly returns one of the specified values.
///
/// # Parameters
///
/// * `a` - the first value
/// * `b` - the second value
pub fn a_or_b<T>(a: T, b: T) -> T {
    if thread_rng().gen() {
        a
    } else {
        b
    }
}

/// Randomly calls one of the specified functions.
///
/// # Parameters
///
/// * `a` - the first function
/// * `b` - the second function
pub fn do_a_or_b<F,G,T>(a: F, b: G) -> T where
    F: FnOnce() -> T,
    G: FnOnce() -> T {
    if thread_rng().gen() {
        a()
    } else {
        b()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
/// An `Iteration` is a sequential datatype for determining the absolute difference between
/// to objects of an iterative process, where the maximum distance is capped.
pub struct Iteration {
    current_iteration: i128,
}

impl Iteration {
    /// Creates a new `Iteration`.
    pub fn new() -> Self {
        Self{current_iteration: 0}
    }

    /// Increments the `Iteration` to its sequentially next step.
    pub fn increment(&self) -> Self {
        Self{current_iteration: self.current_iteration + 1}
    }

    /// Returns the difference in steps between both `Iteration`s. If the difference is out
    /// of range of the returned datatype, the corresponding datatype's maximum or minimum
    /// value is returned.
    ///
    /// # Parameters
    ///
    /// * `other` - the `Iteration` to calculate the difference from
    pub fn difference<T: std::borrow::Borrow<Iteration>> (&self, other: T) -> i32 {
        let real_difference = self.current_iteration - other.borrow().current_iteration;
        i32::try_from(real_difference).unwrap_or_else(|_| -> i32 {
            if real_difference > i32::MAX as i128 {
                i32::MAX
            } else {
                i32::MIN
            }
        })
    }
}

impl std::ops::Sub for Iteration {
    type Output = i32;

    fn sub(self, other: Self) -> Self::Output {
        self.difference(other)
    }
}

impl std::ops::Sub for &Iteration {
    type Output = i32;

    fn sub(self, other: Self) -> Self::Output {
        self.difference(other)
    }
}

impl Iterator for Iteration {
    type Item = Iteration;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.increment())
    }
}

#[cfg(test)]
mod tests;
