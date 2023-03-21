use std::{collections::HashMap, sync::Arc};

use parking_lot::Mutex;

use crate::evolution::helper::Iteration;

use super::{network::ErrorPropagationImpulse, neuron::Neuron, parameters::ConfigurableParameters};

pub struct Dendrite {
    weight: Mutex<f64>,
    target: Arc<Neuron>,
    source: Arc<Neuron>,
    times_activated: Mutex<u64>,
    propagation_errors: Mutex<HashMap<usize, f64>>,
    inhibitory: bool,
    configuration: Arc<ConfigurableParameters>,
}

impl Dendrite {
    /// Creates a new `Dendrite`.
    ///
    /// # Parameters
    ///
    /// * `target` - the neuron targeted by this `Dendrite`
    /// * `target` - the neuron that is the source of this `Dendrite`
    /// * `weight` - the strenght of the dendrite influence on the target neuron  
    /// * `inhibitory` - if the `Dendrite` is inhibitory
    /// * `configuration` - a set of configuration parameters
    pub fn new(
        target: Arc<Neuron>,
        source: Arc<Neuron>,
        weight: f64,
        inhibitory: bool,
        configuration: Arc<ConfigurableParameters>,
    ) -> Self {
        let normalised_weight = Self::normalise_weight(weight);
        Self {
            weight: Mutex::new(normalised_weight),
            target,
            source,
            times_activated: Mutex::new(0),
            propagation_errors: Mutex::new(HashMap::new()),
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

    /// Sets the weight of the `Dendrite` to the specified value
    /// and includes some normalisation based on the previous corrections.
    ///
    /// # Parameters
    ///
    /// * `weight` - the value to add to the weight
    pub fn add_weight_after_evaluation(&self, weight: f64) {
        self.set_weight(self.weight() + weight);
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

    /// Returns the target [`Neuron`] of the dendrite.
    pub fn target(&self) -> Arc<Neuron> {
        // TODO: think about where to use weak references to prevent memory leaks
        Arc::clone(&self.target)
    }

    /// Returns the source [`Neuron`] of the dendrite.
    pub fn source(&self) -> Arc<Neuron> {
        // TODO: think about where to use weak references to prevent memory leaks
        Arc::clone(&self.source)
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
    pub fn trigger(&self, network_activity: f64, time: Iteration) {
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
            // * scaling_factor
        } else {
            -current_weight
                * self
                    .configuration
                    .dendrite_global_activity_regulation_scaling_depression()
            // / scaling_factor
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

    pub fn propagate_error(&self, error: ErrorPropagationImpulse, total_dendrite_activations: u64) {
        let mut error_map = self.propagation_errors.lock();
        if !error_map.contains_key(&error.output_id()) {
            error_map.insert(error.output_id(), error.impact_factor() * error.error());
            drop(error_map);
            self.source()
                .propagate_error(error, total_dendrite_activations);
        }
    }

    pub fn evaluate_error_propagation(&self) {
        let mut error_map = self.propagation_errors.lock();
        let error_sum: f64 = error_map.values().sum();
        // println!("{}", error_sum);
        if self.is_inhibitory() {
            self.set_weight(self.weight() + error_sum);
        } else {
            self.set_weight(self.weight() - error_sum);
        }
        error_map.clear();
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
