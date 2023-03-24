use std::{collections::HashMap, sync::{Arc, Weak}};

use parking_lot::Mutex;

use super::{network::ErrorPropagationImpulse, neuron::Neuron, parameters::ConfigurableParameters};

const NEURON_WEAK_UPGRADE_ERROR: &str = "Neurons are not dropped, so it must be present.";

pub struct Dendrite {
    weight: Mutex<f64>,
    target: Weak<Neuron>,
    source: Weak<Neuron>,
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
            target: Arc::<Neuron>::downgrade(&target),
            source: Arc::<Neuron>::downgrade(&source),
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

    /// Returns the target [`Neuron`] of the dendrite.
    pub fn target(&self) -> Arc<Neuron> {
        Arc::clone(&self.target.upgrade().expect(NEURON_WEAK_UPGRADE_ERROR))
    }

    /// Returns the source [`Neuron`] of the dendrite.
    pub fn source(&self) -> Arc<Neuron> {
        Arc::clone(&self.source.upgrade().expect(NEURON_WEAK_UPGRADE_ERROR))
    }

    /// Returns the triggering frequency of this dendrite based on its weight and source [`Neuron`].
    pub fn triggering_frequency(&self) -> f64{
        let activation_potential_frequency = self.source().value() * self.weight();
        if self.is_inhibitory() {
            -activation_potential_frequency
        } else { 
            activation_potential_frequency
        }
    }

    pub fn propagate_error(&self, error: ErrorPropagationImpulse, b: f64) {
        let mut error_map = self.propagation_errors.lock();
        if !error_map.contains_key(&error.output_id()) {
            println!("Weighted error: {} ; error: {} ; distance: {}; b: {}", error.error() * (1.0 / b) * 0.99999f64.powi(error.distance() as i32), error.error(), error.distance(), b);
            error_map.insert(error.output_id(), error.error() * (1.0 / b) * 0.99f64.powi(error.distance() as i32));
            drop(error_map);
            self.source()
                .propagate_error(error);
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
