use derive_builder::Builder;
use rand::{thread_rng, Rng};
use rand_distr::{Distribution, Standard};
use serde::{Deserialize, Serialize};

use crate::evolution::{chemistry::Information, gene::CrossOver, helper::Nlbf64};

#[derive(Builder, Copy, Clone, Serialize, Deserialize, PartialEq)]
#[builder(default)]
/// This struct provides configurable parameters for neural networks.
pub struct ConfigurableParameters {
    /// The base value / potential of a [`Neuron`].
    neuron_base_value: Nlbf64,
    /// The threshold for triggering and activation potential.
    neuron_activation_threshold: Nlbf64,
    /// The halflife time of the value / potential change of a [`Neuron`]
    /// from the base level in [`Iteration`]s when above the base level.
    neuron_value_halflife: Nlbf64,
    /// The halflife time of the value / potential change of a [`Neuron`]
    /// from the base level in [`Iteration`]s when below the base level.
    neuron_value_halflife_refractory: Nlbf64,
    /// The limit below which no value alterations can occur.
    neuron_refractory_limit: Nlbf64,
    /// The value that the connection weight is reinforced if an activation coincides with
    /// an action potential in the [`Neuron`].
    dendrite_activation_potential_reinforcement: Nlbf64,
    /// The value that the connection weight is depressed if an activation coincides with
    /// a [`Neuron`] in the refractory state.
    dendrite_activation_potential_depression: Nlbf64,
    /// The relative activity of the network where no weight scaling is performed.
    dendrite_global_activity_regulation_midpoint: Nlbf64,
    /// The exponent of relative network activity related weight scaling.
    dendrite_global_activity_regulation_exponent: Nlbf64,
    /// The scaing factor applied to positive [`Dendrite`] weights.
    dendrite_global_activity_regulation_scaling_reinforcment: Nlbf64,
    /// The scaling factor applied to negative [`Dendrite`] weights.
    dendrite_global_activity_regulation_scaling_depression: Nlbf64,
}

impl ConfigurableParameters {
    pub fn neuron_base_value(&self) -> f64 {
        self.neuron_base_value.value()
    }

    pub fn neuron_activation_threshold(&self) -> f64 {
        self.neuron_activation_threshold.value()
    }

    pub fn neuron_value_halflife(&self) -> f64 {
        Self::nlbf64_to_range(self.neuron_value_halflife, 0.0000000001, 10000.0)
    }

    pub fn neuron_value_halflife_refractory(&self) -> f64 {
        Self::nlbf64_to_range(self.neuron_value_halflife_refractory, 0.0000000001, 10000.0)
    }

    pub fn neuron_refractory_limit(&self) -> f64 {
        self.neuron_refractory_limit.value()
    }

    pub fn dendrite_activation_potential_reinforcement(&self) -> f64 {
        Self::nlbf64_to_range(self.dendrite_activation_potential_reinforcement, -1.0, 1.0)
    }

    pub fn dendrite_activation_potential_depression(&self) -> f64 {
        Self::nlbf64_to_range(self.dendrite_activation_potential_depression, -1.0, 1.0)
    }

    pub fn dendrite_global_activity_regulation_midpoint(&self) -> f64 {
        self.dendrite_global_activity_regulation_midpoint.value()
    }

    pub fn dendrite_global_activity_regulation_exponent(&self) -> f64 {
        Self::nlbf64_to_range(self.dendrite_global_activity_regulation_exponent, -20.0, 20.0)
    }

    pub fn dendrite_global_activity_regulation_scaling_reinforcment(&self) -> f64 {
        Self::nlbf64_to_range(
            self.dendrite_global_activity_regulation_scaling_reinforcment,
            0.0,
            10000.0,
        )
    }

    pub fn dendrite_global_activity_regulation_scaling_depression(&self) -> f64 {
        Self::nlbf64_to_range(
            self.dendrite_global_activity_regulation_scaling_depression,
            0.0,
            10000.0,
        )
    }

    pub fn mutate(&self) -> Self {
        let mut mutated = self.clone();
        match thread_rng().gen_range(0..=10u8) {
            0 => mutated.neuron_base_value = Nlbf64::flip_random_bit(mutated.neuron_base_value),
            1 => {
                mutated.neuron_activation_threshold =
                    Nlbf64::flip_random_bit(mutated.neuron_activation_threshold)
            },
            2 => {
                mutated.neuron_value_halflife =
                    Nlbf64::flip_random_bit(mutated.neuron_value_halflife)
            },
            3 => {
                mutated.neuron_value_halflife_refractory =
                    Nlbf64::flip_random_bit(mutated.neuron_value_halflife_refractory)
            },
            4 => {
                mutated.neuron_refractory_limit =
                    Nlbf64::flip_random_bit(mutated.neuron_refractory_limit)
            },
            5 => {
                mutated.dendrite_activation_potential_reinforcement =
                    Nlbf64::flip_random_bit(mutated.dendrite_activation_potential_reinforcement)
            },
            6 => {
                mutated.dendrite_activation_potential_depression =
                    Nlbf64::flip_random_bit(mutated.dendrite_activation_potential_depression)
            },
            7 => {
                mutated.dendrite_global_activity_regulation_midpoint =
                    Nlbf64::flip_random_bit(mutated.dendrite_global_activity_regulation_midpoint)
            },
            8 => {
                mutated.dendrite_global_activity_regulation_exponent =
                    Nlbf64::flip_random_bit(mutated.dendrite_global_activity_regulation_exponent)
            },
            9 => {
                mutated.dendrite_global_activity_regulation_scaling_reinforcment =
                    Nlbf64::flip_random_bit(
                        mutated.dendrite_global_activity_regulation_scaling_reinforcment,
                    )
            },
            10 => {
                mutated.dendrite_global_activity_regulation_scaling_depression =
                    Nlbf64::flip_random_bit(
                        mutated.dendrite_global_activity_regulation_scaling_depression,
                    )
            },
            a => panic!("Out of range: {}", a),
        }
        mutated
    }

    /// Scales the [`Nlbf64`] to the specified boundaries.
    /// 
    /// # Parameters
    /// 
    /// * `lower_bound` - the lowest possible value (inclusive)
    /// * `upper_bound` - the highest possible value (inclusive)
    /// 
    /// # Panics
    /// 
    /// If the lower bound is greater than the upper bound.
    fn nlbf64_to_range(value: Nlbf64, lower_bound: f64, upper_bound: f64) -> f64 {
        if upper_bound <= lower_bound {
            panic!("The lower bound is higher than the upper bound.");
        }
        let difference = upper_bound - lower_bound;
        lower_bound + difference * value.value()
    }
}

impl Information for ConfigurableParameters {
    fn update_value(&mut self, _time_passed: i32) {}
}

impl CrossOver for ConfigurableParameters {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        Self {
            neuron_base_value: self.neuron_base_value.cross_over(&other.neuron_base_value),
            neuron_activation_threshold: self
                .neuron_activation_threshold
                .cross_over(&other.neuron_activation_threshold),
            neuron_value_halflife: self
                .neuron_value_halflife
                .cross_over(&other.neuron_value_halflife),
            neuron_value_halflife_refractory: self
                .neuron_value_halflife_refractory
                .cross_over(&other.neuron_value_halflife_refractory),
            neuron_refractory_limit: self
                .neuron_refractory_limit
                .cross_over(&other.neuron_refractory_limit),
            dendrite_activation_potential_reinforcement: self
                .dendrite_activation_potential_reinforcement
                .cross_over(&other.dendrite_activation_potential_reinforcement),
            dendrite_activation_potential_depression: self
                .dendrite_activation_potential_depression
                .cross_over(&other.dendrite_activation_potential_depression),
            dendrite_global_activity_regulation_midpoint: self
                .dendrite_global_activity_regulation_midpoint
                .cross_over(&other.dendrite_global_activity_regulation_midpoint),
            dendrite_global_activity_regulation_exponent: self
                .dendrite_global_activity_regulation_exponent
                .cross_over(&other.dendrite_global_activity_regulation_exponent),
            dendrite_global_activity_regulation_scaling_reinforcment: self
                .dendrite_global_activity_regulation_scaling_reinforcment
                .cross_over(&other.dendrite_global_activity_regulation_scaling_reinforcment),
            dendrite_global_activity_regulation_scaling_depression: self
                .dendrite_global_activity_regulation_scaling_depression
                .cross_over(&other.dendrite_global_activity_regulation_scaling_depression),
        }
    }
}

impl Default for ConfigurableParameters {
    fn default() -> Self {
        Self {
            neuron_base_value: 0.2.into(),
            neuron_activation_threshold: 0.95.into(),
            neuron_value_halflife: (16000.0 / 10000.0).into(),
            neuron_value_halflife_refractory: (4.0 / 10000.0).into(),
            neuron_refractory_limit: 0.15.into(),
            dendrite_activation_potential_reinforcement: 0.5005.into(),
            dendrite_activation_potential_depression: 0.5005.into(),
            dendrite_global_activity_regulation_midpoint: 0.5.into(),
            dendrite_global_activity_regulation_exponent: 0.58.into(),
            dendrite_global_activity_regulation_scaling_reinforcment: (2.0 / 10000.0).into(),
            dendrite_global_activity_regulation_scaling_depression: (2.0 / 10000.0).into(),
        }
    }
}

impl std::fmt::Debug for ConfigurableParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ConfigurableParameters")
            .field("neuron_base_value", &self.neuron_base_value())
            .field("neuron_activation_threshold", &self.neuron_activation_threshold())
            .field("neuron_value_halflife", &self.neuron_value_halflife())
            .field("neuron_value_halflife_refractory", &self.neuron_value_halflife_refractory())
            .field("neuron_refractory_limit", &self.neuron_refractory_limit())
            .field(
                "dendrite_activation_potential_reinforcement",
                &self.dendrite_activation_potential_reinforcement(),
            )
            .field(
                "dendrite_activation_potential_depression",
                &self.dendrite_activation_potential_depression(),
            )
            .field(
                "dendrite_global_activity_regulation_midpoint",
                &self.dendrite_global_activity_regulation_midpoint(),
            )
            .field(
                "dendrite_global_activity_regulation_exponent",
                &self.dendrite_global_activity_regulation_exponent(),
            )
            .field(
                "dendrite_global_activity_regulation_scaling_reinforcment",
                &self.dendrite_global_activity_regulation_scaling_reinforcment(),
            )
            .field(
                "dendrite_global_activity_regulation_scaling_depression",
                &self.dendrite_global_activity_regulation_scaling_depression(),
            )
            .finish()
    }
}

impl Distribution<ConfigurableParameters> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ConfigurableParameters {
        ConfigurableParameters {
            neuron_base_value: rng.gen(),
            neuron_activation_threshold: rng.gen(),
            neuron_value_halflife: rng.gen(),
            neuron_value_halflife_refractory: rng.gen(),
            neuron_refractory_limit: rng.gen(),
            dendrite_activation_potential_reinforcement: rng.gen(),
            dendrite_activation_potential_depression: rng.gen(),
            dendrite_global_activity_regulation_midpoint: rng.gen(),
            dendrite_global_activity_regulation_exponent: rng.gen(),
            dendrite_global_activity_regulation_scaling_reinforcment: rng.gen(),
            dendrite_global_activity_regulation_scaling_depression: rng.gen(),
        }
    }
}

unsafe impl Send for ConfigurableParameters {}

unsafe impl Sync for ConfigurableParameters {}
