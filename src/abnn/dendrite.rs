use std::sync::Arc;

use parking_lot::Mutex;

use crate::evolution::helper::Iteration;

use super::neuron::Neuron;

/// The value that the connection weight is reinforced if an activation coincides with
/// an action potential in the [`Neuron`].
pub const ACTIVATION_POTENTIAL_REINFORCEMENT: f64 = 0.01;

pub struct Dendrite {
    weight: Mutex<f64>,
    target: Arc<Neuron>,
    times_activated: Mutex<u64>,
}

impl Dendrite {
    /// Creates a new `Dendrite`.
    ///
    /// # Parameters
    ///
    /// * `target` - the neuron targeted by this `Dendrite`
    /// * `weight` - the strenght of the dendrite influence on the target neuron  
    pub fn new(target: Arc<Neuron>, weight: f64) -> Self {
        let normalised_weight = if weight.is_nan() {
            0.0
        } else if weight <= -1.0 {
            -1.0
        } else if weight >= 1.0 {
            1.0
        } else {
            weight
        };
        Self {
            weight: Mutex::new(normalised_weight),
            target,
            times_activated: Mutex::new(0),
        }
    }

    /// Returns the weight of the dendrite connection.
    pub fn weight(&self) -> f64 {
        *self.weight.lock()
    }

    /// Sets the weight of the `Dendrite` to the specified value.
    ///
    /// # Parameters
    ///
    /// * `weight` - the new weight value to set
    pub fn set_weight(&self, weight: f64) {
        let normalised_weight = if weight.is_nan() {
            0.0
        } else if weight <= -1.0 {
            -1.0
        } else if weight >= 1.0 {
            1.0
        } else {
            weight
        };
        *self.weight.lock() = normalised_weight;
    }

    /// Returns the target [`Neuron`].
    pub fn target(&self) -> Arc<Neuron> {
        Arc::clone(&self.target)
    }

    /// Returns how often the dendrite was activated.
    pub fn times_activated(&self) -> u64 {
        *self.times_activated.lock()
    }

    /// Sets the amount of timesthe `Dendrite` was activated.
    ///
    /// # Parameters
    ///
    /// * `value` - how many times the `Dendrite` was activated
    pub fn set_times_activated(&self, value: u64) {
        *self.times_activated.lock() = value;
    }

    /// Triggers the dendrite passing the weighted action potential on to its target [`Neuron`].
    ///
    /// # Parameters
    ///
    /// * `time` - the current timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    pub fn trigger(&self, time: Iteration) {
        self.target.add_value(self.weight(), time);
        self.set_times_activated(self.times_activated() + 1);
    }
}

unsafe impl Send for Dendrite {}

unsafe impl Sync for Dendrite {}
