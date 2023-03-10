use std::sync::Arc;

use parking_lot::Mutex;

use crate::evolution::helper::Iteration;

use super::{dendrite::Dendrite, parameters::ConfigurableParameters};

/// A neuron.
pub struct Neuron {
    value: Mutex<f64>,
    last_value_update: Mutex<Iteration>,
    dendrites: Mutex<Vec<Arc<Dendrite>>>,
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
            last_value_update: Mutex::new(Iteration::new()),
            dendrites: Mutex::new(Vec::new()),
            configuration,
        }
    }

    /// Returns the current value of the `Neuron`.
    fn value(&self) -> f64 {
        *self.value.lock()
    }

    /// Updates the `Neuron`'s value based on the specified time.
    ///
    /// # Parameters
    ///
    /// * `time` - the new timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    fn update_value(&self, time: Iteration) {
        let time_difference = time - self.last_value_update();
        let current_value = self.value();
        if time_difference < 0 {
            panic!(
                "The specified timepoint {:?} is before the last update timepoint {:?}.",
                time,
                self.last_value_update()
            );
        } else if time_difference > 0 {
            if current_value != self.configuration.neuron_base_value() {
                let halflife = if current_value < self.configuration.neuron_base_value() {
                    self.configuration.neuron_value_halflife_refractory()
                } else {
                    self.configuration.neuron_value_halflife()
                };
                let value_difference: f64 = current_value - self.configuration.neuron_base_value();
                let change: f64 =
                    value_difference * 0.5f64.powf(time_difference as f64 / halflife);
                self.set_value(self.configuration.neuron_base_value() + change);
            }
            self.set_time(time);
        }
    }

    /// Returns the value of the `Neuron` at the specified timepoint.
    /// This will update the neurons value.
    ///
    /// # Parameters
    ///
    /// * `time` - the new timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    pub fn value_at_timepoint(&self, time: Iteration) -> f64 {
        self.update_value(time);
        self.value()
    }

    /// Returns the timepoint when the value of the `Neuron` was last updated.
    pub fn last_value_update(&self) -> Iteration {
        *self.last_value_update.lock()
    }

    /// Sets the [`Neuron`]'s value as specified without updating the time.
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

    /// Adds a dendrite to the [`Neuron`].
    ///
    /// # Parameters
    ///
    /// * `dendrite` - the dendrite to add
    pub fn add_dendrite(&self, dendrite: Arc<Dendrite>) {
        self.dendrites.lock().push(dendrite);
    }

    pub fn dendrites(&self) -> Vec<Arc<Dendrite>> {
        self.dendrites.lock().clone()
    }

    /// Checks if the `Neuron`'s potential is high enough to trigger an activation potential.
    ///
    /// # Parameters
    ///
    /// * `time` - the current timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    pub fn should_trigger_activation_potential(&self, time: Iteration) -> bool {
        self.value_at_timepoint(time) >= self.configuration.neuron_activation_threshold()
    }

    /// Tries to trigger an action potential and returns the influenced dendrites if successful.
    ///
    /// # Parameters
    ///
    /// * `time` - the current timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    pub fn try_trigger_action_potential(&self, time: Iteration) -> Vec<Arc<Dendrite>> {
        if self.should_trigger_activation_potential(time) {
            // Hyperpolarisation.
            self.set_value(0.0);
            self.dendrites()
        } else {
            Vec::new()
        }
    }

    /// Sets the timepoint the [`Neuron`] was last updated without checking for the timepoint it was last updated.
    ///
    /// # Parameters
    ///
    /// * `time` - the timepoint to set
    pub fn set_time(&self, time: Iteration) {
        *self.last_value_update.lock() = time;
    }

    /// Adds the specified value at the specified timepoint and returns 
    /// `true` if the addition was successfull or `false` if the `Neuron` was hyperpolarised.
    ///
    /// # Parameters
    ///
    /// * `value` - the value to add
    /// * `time` - the new timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    pub fn add_value(&self, value: f64, time: Iteration) -> bool {
        let current_value = self.value_at_timepoint(time);
        if current_value > self.configuration.neuron_refractory_limit() {
            self.set_value(current_value + value);
            false
        } else {
            true
        }
    }
}

unsafe impl Send for Neuron {}

unsafe impl Sync for Neuron {}

#[cfg(test)]
mod tests;
