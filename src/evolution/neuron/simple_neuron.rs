//! The `neuron` module contains consturcts for working with neuronal networks.
extern crate bitvec;

use crate::evolution::helper::Nlbf64;
use rand::Rng;
use rand_distr::{Distribution, Standard};
use serde::{Deserialize, Serialize};

use super::super::chemistry::Information;
use super::super::gene::CrossOver;
use std::ops::Add;

/// The halflife time of potential changes (from the base potential) in iteration steps.
const POTENTIAL_HALFLIFE_TIME: f64 = 16.0;

/// A `SimpleNeuron` is the basic representation of a biological neuron with an excitation and a base potential.
#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleNeuron {
    /// The current potential of the neuron.
    current_potential: Nlbf64,
    /// The base potential represents the initial value as well as the value that the neuron will settle on after sufficient time.
    base_potential: Nlbf64,
}

impl SimpleNeuron {
    /// Creates a new `SimpleNeuron` with the same current and base potential.
    ///
    /// # Parmeters
    ///
    /// `potential` - the base and starting potential of the neuron
    pub fn new<N: Into<Nlbf64>>(potential: N) -> Self {
        let converted_potential: Nlbf64 = potential.into();
        SimpleNeuron {
            current_potential: converted_potential,
            base_potential: converted_potential,
        }
    }

    /// Creates a new `SimpleNeuron` with the specified current potential
    /// and the base potential of this neuron.
    ///
    /// # Parmeters
    ///
    /// `value` - the current potential of the created neuron
    pub fn with_new_current_potential<N: Into<Nlbf64>>(&self, value: N) -> Self {
        SimpleNeuron {
            current_potential: value.into(),
            base_potential: self.base_potential(),
        }
    }

    /// Returns the current potential of this neuron.
    pub fn current_potential(&self) -> Nlbf64 {
        self.current_potential
    }

    /// Returns the base potential of this neuron.
    /// This defines the initial state of the neuron as well as the value the neuron
    /// will settle on after sufficient time passed.
    pub fn base_potential(&self) -> Nlbf64 {
        self.base_potential
    }

    /// Returns a random neuron in its initial state.
    pub fn random() -> Self {
        SimpleNeuron::new(f64::from_be_bytes(rand::random()))
    }
}

impl Add for SimpleNeuron {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.with_new_current_potential(self.current_potential() + other.current_potential())
    }
}

impl Distribution<SimpleNeuron> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> SimpleNeuron {
        SimpleNeuron::new(rng.gen::<Nlbf64>())
    }
}

impl CrossOver for SimpleNeuron {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        let recombined_value = self.base_potential().cross_over(&other.base_potential());
        SimpleNeuron::new(recombined_value)
    }
}

impl Information for SimpleNeuron {
    fn update_value(&mut self, time_passed: i32) {
        if self.current_potential != self.base_potential && time_passed != 0 {
            let diff_potential: f64 = self.current_potential.value() - self.base_potential.value();
            let change: f64 = diff_potential * 0.5f64.powf(time_passed as f64 / POTENTIAL_HALFLIFE_TIME);
            self.current_potential = (self.base_potential.value() + change).into();
        }
    }
}

#[cfg(test)]
mod tests;
