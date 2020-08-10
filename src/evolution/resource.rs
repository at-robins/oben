//! The `resource` module contains the growth limiting framework of the evolutionary network.
extern crate serde;

use serde::{Serialize, Deserialize};


#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
/// A `Resource` is needed for [`ClonalPopulation`]s to produce offspring.
///
/// [`ClonalPopulation`]: ../population/struct.Population.html
pub struct Resource {
    available: f64,
    recycling: f64,
    half_life: f64,
}

impl Resource {
    /// Creates a new `Resource` with all resources being completely available initially.
    /// Inavailable `Resource`s are recycled with the specified half life time.
    ///
    /// # Parameters
    ///
    /// * `total` - the total amount of `Resource`s available
    /// * `half_life` - the half life time for recycling `Resource`s in generations
    pub fn new(total: f64, half_life: f64) -> Self {
        Resource {
            available: total,
            recycling: 0.0,
            half_life,
        }
    }

    /// Returns the total amount of `Resource`s, combining available `Resource`s and
    /// `Resources`, which are currently recycled.
    pub fn total(&self) -> f64 {
        self.recycling + self.available
    }

    /// Returns the amount of available `Resource`s.
    pub fn available(&self) -> f64 {
        self.available
    }

    /// Returns the amount of inavailable `Resource`s that are currently in the process of
    /// being recycled.
    pub fn recycling(&self) -> f64 {
        self.recycling
    }

    /// Returns the time in generations it takes for half of the inavailable `Resource`s
    /// to be converted back into available `Resource`s.
    pub fn half_life(&self) -> f64 {
        self.half_life
    }

    /// Recycle inavailable `Resource`s by forwarding the time a single generation.
    pub fn recycle(&mut self) {
        if self.recycling() > 0.0 {
            let recycled = self.recycling() * (1.0 - 0.5f64.powf(1.0 / self.half_life()));
            self.recycling -= recycled;
            self.available += recycled;
        }
    }

    /// Retrieve the specified amount of `Resource`s if available, otherwise return the maximum
    /// available amount.
    ///
    /// # Parameters
    ///
    /// * `amount` - the amount of `Resource`s to be claimed
    ///
    /// # Panics
    ///
    /// If the specified `amount` is not a valid positive number.
    pub fn claim_resources(&mut self, amount: f64) -> f64 {
        if amount >= 0.0 {
            if amount >= self.available() {
                let all_resources = self.available();
                self.available = 0.0;
                all_resources
            } else {
                self.available -= amount;
                amount
            }
        } else {
            panic!("{} resources cannot be claimed.", amount);
        }
    }

    /// Add the specified amount of `Resource`s.
    ///
    /// # Parameters
    ///
    /// * `amount` - the amount of `Resource`s to add
    ///
    /// # Panics
    ///
    /// If the specified `amount` is not a valid positive number.
    pub fn repatriate_resources(&mut self, amount: f64) {
        if amount >= 0.0 {
            self.recycling += amount;
        } else {
            panic!("{} resources cannot be repatriated.", amount);
        }
    }
}

impl Default for Resource {
    fn default() -> Self {
        Resource::new(100_000.0, 3.0)
    }
}
