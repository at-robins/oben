use std::sync::{Arc, Weak};

use parking_lot::Mutex;

use super::{
    dendrite::Dendrite, network::ErrorPropagationImpulse, parameters::ConfigurableParameters,
};

const DENDRITE_WEAK_UPGRADE_ERROR: &str = "Dendrites are not dropped, so it must be present.";

/// A neuron.
pub struct Neuron {
    /// The relative frequency of action potential triggering.
    value: Mutex<f64>,
    outgoing_dendrites: Mutex<Vec<Arc<Dendrite>>>,
    ingoing_dendrites: Mutex<Vec<Weak<Dendrite>>>,
    configuration: Arc<ConfigurableParameters>,
}

impl Neuron {
    /// Creates a new `Neuron`.
    ///
    /// # Parameters
    ///
    /// * `configuration` - a set of configuration parameters
    pub fn new(configuration: Arc<ConfigurableParameters>) -> Self {
        Self {
            value: Mutex::new(configuration.neuron_base_value()),
            outgoing_dendrites: Mutex::new(Vec::new()),
            ingoing_dendrites: Mutex::new(Vec::new()),
            configuration,
        }
    }

    /// Returns the current relative action potential triggering frequency of the `Neuron`.
    pub fn value(&self) -> f64 {
        *self.value.lock()
    }

    /// Updates the `Neuron`'s value based on its ingoing connections.
    pub fn update_value(&self) -> bool {
        let current_value = self.value();
        self.set_value(
            self.ingoing_dendrites
                .lock()
                .iter()
                .map(|dendrite| {
                    dendrite
                        .upgrade()
                        .expect(DENDRITE_WEAK_UPGRADE_ERROR)
                        .triggering_frequency()
                })
                .sum(),
        );
        current_value != self.value()
    }

    /// Sets the [`Neuron`]'s value as specified.
    ///
    /// # Parameters
    ///
    /// * `value` - the new value
    pub fn set_value(&self, value: f64) {
        let new_value = if value <= 0.0 || value.is_nan() {
            0.0
        } else if value >= 1.0 {
            1.0
        } else {
            value
        };
        *self.value.lock() = new_value;
    }

    /// Adds an outgoing dendrite to the [`Neuron`].
    ///
    /// # Parameters
    ///
    /// * `dendrite` - the dendrite to add
    pub fn add_outgoing_dendrite(&self, dendrite: Arc<Dendrite>) {
        self.outgoing_dendrites.lock().push(dendrite);
    }

    /// Adds an ingoing dendrite to the [`Neuron`].
    ///
    /// # Parameters
    ///
    /// * `dendrite` - the dendrite to add
    pub fn add_ingoing_dendrite(&self, dendrite: Arc<Dendrite>) {
        self.ingoing_dendrites.lock().push(Arc::<Dendrite>::downgrade(&dendrite));
    }

    /// Returns the outgoing [`Dendrite`]s.
    pub fn outgoing_dendrites(&self) -> Vec<Arc<Dendrite>> {
        self.outgoing_dendrites.lock().clone()
    }

    /// Returns the ingoing [`Dendrite`]s.
    pub fn ingoing_dendrites(&self) -> Vec<Arc<Dendrite>> {
        self.ingoing_dendrites
            .lock()
            .iter()
            .map(|dendrite| dendrite.upgrade().expect(DENDRITE_WEAK_UPGRADE_ERROR))
            .collect()
    }

    /// Returns all [`Neuron`]s targeted by this [`Neuron`] by [`Dendrite`]s with a weight greater than `0.0`.
    pub fn targets(&self) -> Vec<Arc<Neuron>> {
        self.outgoing_dendrites
            .lock()
            .iter()
            //.filter(|dendrite| dendrite.weight().is_normal())
            .map(|dendrite| dendrite.target())
            .collect()
    }

    pub fn propagate_error(&self, error: ErrorPropagationImpulse) {
        // TODO: rework error propagation
        let ingoing_dendrites = self.ingoing_dendrites();
        let number_of_ingoing_dendrites = ingoing_dendrites.len() as f64;
        for dendrite in ingoing_dendrites.into_iter() {
            dendrite.propagate_error(error.with_updated_distance(), number_of_ingoing_dendrites);
        }
    }
}

unsafe impl Send for Neuron {}

unsafe impl Sync for Neuron {}

#[cfg(test)]
mod tests;
