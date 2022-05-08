use std::sync::Arc;

use parking_lot::Mutex;

use crate::evolution::helper::Iteration;

use super::dendrite::Dendrite;

/// The base value / potential of a [`Neuron`].
pub const NEURON_BASE_VALUE: f64 = 0.2;
/// The threshold for triggering and activation potential.
pub const NEURON_ACTIVATION_POTENTIAL_THRESHOLD: f64 = 0.95;
/// The halflife time of the value / potential change of a [`Neuron`]
/// from the base level in [`Iteration`]s.
pub const NEURON_VALUE_HALFLIFE: f64 = 16.0;

/// A neuron.
pub struct Neuron {
    value: Mutex<f64>,
    last_value_update: Mutex<Iteration>,
    dendrites: Mutex<Vec<Arc<Dendrite>>>,
}

impl Neuron {
    /// Creates a new `Neuron`.
    pub fn new() -> Self {
        Self {
            value: Mutex::new(NEURON_BASE_VALUE),
            last_value_update: Mutex::new(Iteration::new()),
            dendrites: Mutex::new(Vec::new()),
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
            if current_value != NEURON_BASE_VALUE {
                let value_difference: f64 = current_value - NEURON_BASE_VALUE;
                let change: f64 =
                    value_difference * 0.5f64.powf(time_difference as f64 / NEURON_VALUE_HALFLIFE);
                self.set_value(NEURON_BASE_VALUE + change);
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
    fn set_value(&self, value: f64) {
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
        self.value_at_timepoint(time) >= NEURON_ACTIVATION_POTENTIAL_THRESHOLD
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
    fn set_time(&self, time: Iteration) {
        *self.last_value_update.lock() = time;
    }

    /// Adds the specified value at the specified timepoint.
    ///
    /// # Parameters
    ///
    /// * `value` - the value to add
    /// * `time` - the new timepoint
    ///
    /// # Panics
    ///
    /// If the specified timepoint is earlier than the last update.
    pub fn add_value(&self, value: f64, time: Iteration) {
        let current_value = self.value_at_timepoint(time);
        self.set_value(current_value + value);
    }
}

/*impl Clone for Neuron {
    fn clone(&self) -> Self {
        let cloned_value = self.value();
        let cloned_last_value_update = self.last_value_update();
        Self {
            value: Mutex::new(cloned_value),
            last_value_update: Mutex::new(cloned_last_value_update),
            dendrites: Mutex::new(self.dendrites()),
            last_action_potential: Mutex::new(self.last_action_potential()),
        }
    }
}

impl std::fmt::Debug for Neuron {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Neuron")
            .field("value", &self.value())
            .field("last_value_update", &self.last_value_update())
            .field("dendrites", &self.dendrites())
            .field("last_action_potential", &self.last_action_potential())
            .finish()
    }
}

impl PartialEq for Neuron {
    fn eq(&self, other: &Self) -> bool {
        let dendrite_equality: bool = self
            .dendrites()
            .into_iter()
            .zip(other.dendrites().into_iter())
            .all(|(a, b)| Arc::ptr_eq(&a, &b));
        self.value() == other.value()
            && self.last_value_update() == other.last_value_update()
            && dendrite_equality
            && self.last_action_potential() == other.last_action_potential()
    }
} */

unsafe impl Send for Neuron {}

unsafe impl Sync for Neuron {}

#[cfg(test)]
mod tests;
