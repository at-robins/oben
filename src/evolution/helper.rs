//! The `helper` module contains helper constructs for general workflow.
extern crate rand;

use rand::{Rng, thread_rng};

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


#[cfg(test)]
mod tests;
