use std::sync::Arc;

use parking_lot::Mutex;

use crate::evolution::helper::Iteration;

use super::{neuron::Neuron, parameters::ConfigurableParameters};

pub struct Dendrite {
    weight: Mutex<f64>,
    target: Arc<Neuron>,
    times_activated: Mutex<u64>,
    weight_modifier: Mutex<f64>,
    inhibitory: bool,
    configuration: Arc<ConfigurableParameters>,
}

impl Dendrite {
    /// Creates a new `Dendrite`.
    ///
    /// # Parameters
    ///
    /// * `target` - the neuron targeted by this `Dendrite`
    /// * `weight` - the strenght of the dendrite influence on the target neuron  
    /// * `inhibitory` - if the `Dendrite` is inhibitory
    /// * `configuration` - a set of configuration parameters
    pub fn new(
        target: Arc<Neuron>,
        weight: f64,
        inhibitory: bool,
        configuration: Arc<ConfigurableParameters>,
    ) -> Self {
        let normalised_weight = Self::normalise_weight(weight);
        Self {
            weight: Mutex::new(normalised_weight),
            target,
            times_activated: Mutex::new(0),
            weight_modifier: Mutex::new(0.0),
            inhibitory,
            configuration,
        }
    }

    /// Checks if the `Dendrite` is inhibitory.
    pub fn is_inhibitory(&self) -> bool {
        self.inhibitory
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
        *self.weight.lock() = Self::normalise_weight(weight);
    }

    /// Returns the current weight modifier of the dendrite connection.
    pub fn weight_modifier(&self) -> f64 {
        *self.weight_modifier.lock()
    }

    /// Sets the weight modifier of the `Dendrite` to the specified value.
    ///
    /// # Parameters
    ///
    /// * `weight_modifier` - the new weight modifier value to set
    pub fn set_weight_modifier(&self, weight_modifier: f64) {
        *self.weight_modifier.lock() = Self::normalise_weight(weight_modifier);
    }

    /// Sets the weight of the `Dendrite` to the specified value
    /// and includes some normalisation based on the previous corrections.
    ///
    /// # Parameters
    ///
    /// * `weight` - the value to add to the weight
    pub fn add_weight_after_evaluation(&self, weight: f64) {
        self.set_weight_modifier(self.weight_modifier() * 0.9 + weight * 0.1);
        self.set_weight(self.weight() + self.weight_modifier());
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
    pub fn trigger(&self, source_value: f64, network_activity: f64, time: Iteration) {
        let current_weight = self.weight();
        let scaling_factor = (self
            .configuration
            .dendrite_global_activity_regulation_midpoint()
            .exp()
            / network_activity.exp())
        .powf(
            self.configuration
                .dendrite_global_activity_regulation_exponent(),
        );
        // println!("activity: {} ; scaling factor: {}", network_activity, scaling_factor);
        let scaled_weight = if current_weight == 0.0 {
            0.0
        } else if !self.is_inhibitory() {
            current_weight
                * self
                    .configuration
                    .dendrite_global_activity_regulation_scaling_reinforcment()
                * scaling_factor
        } else {
            -current_weight
                * self
                    .configuration
                    .dendrite_global_activity_regulation_scaling_depression()
                / scaling_factor
        };
        if !self.target.add_value(scaled_weight, time) {
            self.set_weight(
                self.weight()
                    - self
                        .configuration
                        .dendrite_activation_potential_depression(),
            );
        }
        self.set_times_activated(self.times_activated() + 1);
    }

    fn normalise_weight(weight: f64) -> f64 {
        if weight.is_nan() {
            0.0
        } else if weight <= 0.0 {
            0.0
        } else if weight >= 1.0 {
            1.0
        } else {
            weight
        }
    }
}

unsafe impl Send for Dendrite {}

unsafe impl Sync for Dendrite {}
