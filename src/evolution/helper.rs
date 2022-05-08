//! The `helper` module contains helper constructs for general workflow.

use rand::{thread_rng, Rng};
use rand_distr::{Distribution, Standard};
use serde::{Deserialize, Serialize};
use std::{
    borrow::Borrow,
    convert::TryFrom,
    ops::{Add, AddAssign, Sub, SubAssign},
};

use super::{
    binary::{as_u64, flip_random_bit, u64_to_binary},
    gene::CrossOver,
};

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
pub fn do_a_or_b<F, G, T>(a: F, b: G) -> T
where
    F: FnOnce() -> T,
    G: FnOnce() -> T,
{
    if thread_rng().gen() {
        a()
    } else {
        b()
    }
}

/// Returns the supplied value if the specified number is positve and normal;
/// returns `0.0` if not to introduce nonlinearity to the function.
///
/// # Parameters
///
/// * `value` - the number to transform
pub fn nonlinear_normal_positve(value: f64) -> f64 {
    if value.is_normal() && value.is_sign_positive() {
        value
    } else {
        // Introduces non-linearity.
        0.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
/// An `Iteration` is a sequential datatype for determining the absolute difference between
/// two objects of an iterative process, where the maximum distance is capped.
pub struct Iteration {
    current_iteration: i128,
}

impl Iteration {
    /// Creates a new `Iteration`.
    pub fn new() -> Self {
        Self {
            current_iteration: 0,
        }
    }

    /// Increments the `Iteration` to its sequentially next step.
    pub fn increment(&self) -> Self {
        Self {
            current_iteration: self.current_iteration + 1,
        }
    }

    /// Returns the difference in steps between both `Iteration`s. If the difference is out
    /// of range of the returned datatype, the corresponding datatype's maximum or minimum
    /// value is returned.
    ///
    /// # Parameters
    ///
    /// * `other` - the `Iteration` to calculate the difference from
    pub fn difference<T: std::borrow::Borrow<Iteration>>(&self, other: T) -> i32 {
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

impl Sub for Iteration {
    type Output = i32;

    fn sub(self, other: Self) -> Self::Output {
        self.difference(other)
    }
}

impl Sub for &Iteration {
    type Output = i32;

    fn sub(self, other: Self) -> Self::Output {
        self.difference(other)
    }
}

impl Default for Iteration {
    fn default() -> Self {
        Iteration::new()
    }
}

/// An `ActionChain` contains actions that are executed at the same
/// [`Iteration`](crate::evolution::helper::Iteration).
pub struct ActionChain<T> {
    mean_size: usize,
    iteration: Iteration,
    actions: Vec<T>,
}

impl<T: PartialEq> ActionChain<T> {
    /// Creates a new `ActionChain` starting a new
    /// [`Iteration`](crate::evolution::helper::Iteration) cycle.
    pub fn new() -> Self {
        Self {
            mean_size: 0,
            iteration: Iteration::new(),
            actions: Vec::new(),
        }
    }

    /// Returns the current [`Iteration`](crate::evolution::helper::Iteration) the actions
    /// are performed at.
    pub fn current_iteration(&self) -> Iteration {
        self.iteration
    }

    /// Returns all actions and starts a new [`Iteration`](crate::evolution::helper::Iteration)
    /// cycle.
    pub fn pop_actions(&mut self) -> Vec<T> {
        self.mean_size = (self.mean_size + self.actions.len()) / 2;
        self.iteration = self.iteration.increment();
        let mut new_actions = Vec::with_capacity(self.mean_size);
        std::mem::swap(&mut self.actions, &mut new_actions);
        new_actions
    }

    /// Adds an action to the `ActionChain` if the action is not already present and returns
    /// `true` if the action was successfully added.
    ///
    /// # Parameters
    ///
    /// * `action` - the action to add
    pub fn push_action(&mut self, action: T) -> bool {
        if self.actions.contains(&action) {
            false
        } else {
            self.actions.push(action);
            true
        }
    }

    /// Returns `true` if no actions are in the `ActionChain`.
    pub fn is_empty(&self) -> bool {
        self.actions.is_empty()
    }
}

impl<T, I: Borrow<Iteration>> From<I> for ActionChain<T> {
    fn from(starting_iteration: I) -> Self {
        Self {
            mean_size: 0,
            iteration: *starting_iteration.borrow(),
            actions: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
/// A scaling factor to scale the fitness function based on the current mean fitness of the population.
pub struct ScalingFactor {
    factor: i32,
    base: f64,
}

impl ScalingFactor {
    /// Creates a new `ScalingFactor` with the specified base value.
    ///
    /// # Parameters
    ///
    /// * `base` - the base to scale by
    pub fn new(base: f64) -> Self {
        Self { factor: 0, base }
    }

    /// Creates a new `ScalingFactor` with the specified base value and exponent.
    ///
    /// # Parameters
    ///
    /// * `base` - the base to scale by
    /// * `exponent` - the exponent to initialise scaling with
    pub fn new_with_exponent(base: f64, exponent: i32) -> Self {
        Self {
            factor: exponent,
            base,
        }
    }

    /// Returns the scaling exponent.
    pub fn exponent(&self) -> i32 {
        self.factor
    }

    /// Increases the scaling factor.
    pub fn increment(&mut self) {
        self.factor += 1;
    }

    /// Decreases the scaling factor.
    pub fn decrement(&mut self) {
        self.factor -= 1;
    }

    /// Returns the current scaling factor.
    pub fn value(&self) -> f64 {
        self.base.powi(self.factor)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
/// A non linear, equally spaced floating point representation bound between 0 and 1.
pub struct Nlbf64 {
    value: u64,
}

impl Nlbf64 {
    /// The minimum value of 0.
    pub const MIN: Nlbf64 = Nlbf64 { value: u64::MIN };
    /// The maximum value of 1.
    pub const MAX: Nlbf64 = Nlbf64 { value: u64::MAX };

    /// Retruns the value as floating point number.
    pub fn value(&self) -> f64 {
        if self.value == 0 {
            0.0
        } else if self.value == u64::MAX {
            1.0
        } else {
            (self.value as f64) / (u64::MAX as f64)
        }
    }

    /// Retruns the value as internal u64 representation.
    pub fn value_as_u64(&self) -> u64 {
        self.value
    }

    /// Flips a random bit of the [`Nlbf64`] and returns the result.
    ///
    /// # Parameters
    ///
    /// * `base` - the number to mutate
    pub fn flip_random_bit(base: Self) -> Self {
        let binary_base = u64_to_binary(base.value);
        Self {
            value: as_u64(&flip_random_bit(&binary_base)),
        }
    }
}

impl From<f64> for Nlbf64 {
    fn from(value: f64) -> Self {
        if value <= 0.0 {
            Self { value: 0 }
        } else if value >= 1.0 {
            Self { value: u64::MAX }
        } else {
            Self {
                value: (value * u64::MAX as f64) as u64,
            }
        }
    }
}

impl From<&f64> for Nlbf64 {
    fn from(value: &f64) -> Self {
        if *value <= 0.0 {
            Self { value: 0 }
        } else if *value >= 1.0 {
            Self { value: u64::MAX }
        } else {
            Self {
                value: (value * u64::MAX as f64) as u64,
            }
        }
    }
}

impl Add for Nlbf64 {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            value: self.value.checked_add(other.value).unwrap_or(u64::MAX),
        }
    }
}

impl Add for &Nlbf64 {
    type Output = Nlbf64;

    fn add(self, other: Self) -> Self::Output {
        *self + *other
    }
}

impl<T: Borrow<Nlbf64>> AddAssign<T> for Nlbf64 {
    fn add_assign(&mut self, other: T) {
        *self = *self + *other.borrow();
    }
}

impl Sub for Nlbf64 {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            value: self.value.checked_sub(other.value).unwrap_or(0),
        }
    }
}

impl Sub for &Nlbf64 {
    type Output = Nlbf64;

    fn sub(self, other: Self) -> Self::Output {
        *self - *other
    }
}

impl<T: Borrow<Nlbf64>> SubAssign<T> for Nlbf64 {
    fn sub_assign(&mut self, other: T) {
        *self = *self - *other.borrow();
    }
}

impl Default for Nlbf64 {
    fn default() -> Self {
        Self { value: 0 }
    }
}

impl Distribution<Nlbf64> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Nlbf64 {
        Nlbf64 { value: rng.gen() }
    }
}

impl CrossOver for Nlbf64 {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        Nlbf64 {
            value: as_u64(&u64_to_binary(self.value).cross_over(&u64_to_binary(other.value))),
        }
    }
}

pub mod noop;
pub mod testing;
#[cfg(test)]
mod tests;
