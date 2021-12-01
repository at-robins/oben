//! The `neuron` module contains consturcts for working with neuronal networks.
extern crate bitvec;

use serde::{Deserialize, Serialize};
use super::super::binary::{as_f64, f64_to_binary};
use super::super::chemistry::Information;
use super::super::gene::CrossOver;
use super::super::helper::nonlinear_normal_positve;
use std::ops::Add;

/// A `SimpleNeuron` is the basic representation of a biological neuron with an excitation and a base potential.
#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleNeuron {
    /// The current potential of the neuron.
    current_potential: f64,
    /// The base potential represents the initial value as well as the value that the neuron will settle on after sufficient time.
    base_potential: f64,
}

impl SimpleNeuron {
    /// Creates a new `SimpleNeuron` with the same current and base potential.
    pub fn new(potential: f64) -> Self {
        let positive_normal_or_zero: f64 = nonlinear_normal_positve(potential);
        SimpleNeuron{current_potential: positive_normal_or_zero, base_potential: positive_normal_or_zero}
    }

    /// Creates a new `SimpleNeuron` with the specified current potential
    /// and the base potential of this neuron.
    pub fn with_new_current_potential(&self, value: f64) -> Self {
        SimpleNeuron{current_potential: nonlinear_normal_positve(value), base_potential: self.base_potential()}
    }

    /// Returns the current potential of this neuron.
    pub fn current_potential(&self) -> f64 {
        self.current_potential
    }

    /// Returns the base potential of this neuron.
    /// This defines the initial state of the neuron as well as the value the neuron 
    /// will settle on after sufficient time passed.
    pub fn base_potential(&self) -> f64 {
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

impl CrossOver for SimpleNeuron {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        let recombined_binary_value = f64_to_binary(self.base_potential()).cross_over(&f64_to_binary(other.base_potential()));
        let recombined_value = as_f64(&recombined_binary_value);
        SimpleNeuron::new(recombined_value)
    }
}

impl Information for SimpleNeuron {
    fn update_value(&mut self, time_passed: i32) {
        // TODO: Implement
    }
}

#[cfg(test)]
mod tests;