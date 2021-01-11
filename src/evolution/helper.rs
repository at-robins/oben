//! The `helper` module contains helper constructs for general workflow.
extern crate bitvec;
extern crate rand;

use rand::{Rng, thread_rng};
use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use std::cmp::Ordering;

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

/// A type alias for the underlying binary representation of
/// [`Iteration`](oben::evolution::helper::Iteration)s.
pub type BinaryIteration = BitBox<Msb0, usize>;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Iteration {
    current_iteration: BinaryIteration,
}

impl Iteration {
    pub fn new() -> Self {
        Self{current_iteration: BitBox::from_element(0)}
    }

    pub fn increment(&self) -> Self {
        let mut next_value = BitVec::from_boxed_bitslice(self.current_iteration.clone());
        let mut carry = true;
        let mut index = 0;
        while carry {
            match next_value.get(index) {
                None => {
                    next_value.insert(index, carry);
                    carry = false;
                },
                Some(val) => {
                    carry = *val && carry;
                    if carry {
                        next_value.set(index, false);
                    } else {
                        next_value.set(index, true);
                    }
                },
            }
            index += 1;
        }
        Self{current_iteration: next_value.into()}
    }

    // pub fn difference<T: std::borrow::Borrow<Iteration>> (&self, other: T) -> i32 {
    //
    // }
}

impl PartialOrd for Iteration {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Iteration {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.current_iteration.len() > other.current_iteration.len() {
            Ordering::Greater
        } else if self.current_iteration.len() < other.current_iteration.len() {
            Ordering::Less
        } else {
            for (index, value) in self.current_iteration.iter().enumerate().rev() {
                if *value && !other.current_iteration[index] {
                    return Ordering::Greater;
                } else if !value && other.current_iteration[index] {
                    return Ordering::Less;
                }
            }
            Ordering::Equal
        }
    }
}

// impl std::ops::Add for Iteration {
//     type Output = Self;
//
//     fn add(self, other: Self) -> Self::Output {
//         Self {
//             current_iteration: self.current_iteration.overflowing_add(other.current_iteration).0
//         }
//     }
// }

// impl std::ops::Add for &Iteration {
//     type Output = Iteration;
//
//     fn add(self, other: Self) -> Self::Output {
//         Iteration {
//             current_iteration: self.current_iteration.overflowing_add(other.current_iteration).0
//         }
//     }
// }

// impl std::ops::Sub for Iteration {
//     type Output = Self;
//
//     fn sub(self, other: Self) -> Self::Output {
//         let cycle = self.current_iteration.overflowing_div(other.current_iteration);
//         let difference = cycle.0;
//         if cycle.1 {
//
//         }
//         Self {
//             current_iteration: cycle
//         }
//     }
// }

impl Iterator for Iteration {
    type Item = Iteration;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.increment())
    }
}

#[cfg(test)]
mod tests;
