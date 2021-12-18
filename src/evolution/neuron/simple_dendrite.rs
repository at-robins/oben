//! The `simple_dendrite` module contains consturcts for emulating a simplified biological dendrite.
extern crate rand;

use crate::evolution::helper::{a_or_b, Nlbf64};
use rand::{thread_rng, Rng};
use rand_distr::{Distribution, Standard};
use serde::{Deserialize, Serialize};

use super::super::chemistry::{Reaction, State};
use super::super::gene::CrossOver;
use super::super::helper::Iteration;
use super::SimpleNeuron;

/// A threshold triggering an activation potential of a simplified biological neuron.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleDendriteThreshold {
    /// The threshold of the dendrite. Values higher or equal to the threshold will cause an activation potential.
    threshold: Nlbf64,
}

impl SimpleDendriteThreshold {
    /// Creates a new dendrite threshold with the specified threshold.
    ///
    /// # Parameters
    ///
    /// * `threshold` - the threshold value
    pub fn new<N: Into<Nlbf64>>(threshold: N) -> Self {
        SimpleDendriteThreshold {
            threshold: threshold.into(),
        }
    }

    /// Returns the threshold of this dendrite.
    pub fn threshold(&self) -> Nlbf64 {
        self.threshold
    }
}

impl CrossOver for SimpleDendriteThreshold {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        SimpleDendriteThreshold {
            threshold: self.threshold.cross_over(&other.threshold),
        }
    }
}

impl State<SimpleNeuron> for SimpleDendriteThreshold {
    fn get_substrate_number(&self) -> usize {
        1
    }

    fn random() -> Self {
        SimpleDendriteThreshold {
            threshold: thread_rng().gen(),
        }
    }

    fn detect(&self, substrates: &[&SimpleNeuron], _detection_time: Iteration) -> bool {
        substrates[0].current_potential() >= self.threshold()
    }
}

impl Distribution<SimpleDendriteThreshold> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> SimpleDendriteThreshold {
        SimpleDendriteThreshold::new(rng.gen::<Nlbf64>())
    }
}

/// An activiation potential of a simplified biological neuron.
#[derive(Debug, PartialEq, Eq, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleDendriteActivationPotential {
    /// The weight of the activation.
    weight: Nlbf64,
    /// `true` if the signal is inhibitory.
    is_inhibitory: bool,
}

impl SimpleDendriteActivationPotential {
    /// Creates a new activation potential traveling  along a single dendrite with the specified weight.
    ///
    /// # Parameters
    ///
    /// * `weight` - the weight value
    /// * `is_inhibitory` - if the activation potential is inhibitory for the target neuron
    pub fn new<B: Into<bool>, N: Into<Nlbf64>>(weight: N, is_inhibitory: B) -> Self {
        SimpleDendriteActivationPotential {
            weight: weight.into(),
            is_inhibitory: is_inhibitory.into(),
        }
    }

    /// The weight of the activation potential.
    /// This represents the activative force an activation potential possese.
    pub fn weight(&self) -> Nlbf64 {
        self.weight
    }

    /// Returns if the activation potential of this dendrite is inhibitory
    /// for the connected neuron.
    pub fn is_inhibitory(&self) -> bool {
        self.is_inhibitory
    }
}

impl CrossOver for SimpleDendriteActivationPotential {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        let recombined_weight = self.weight.cross_over(&other.weight);
        let recombined_inhibitory_state = a_or_b(self.is_inhibitory, other.is_inhibitory);
        SimpleDendriteActivationPotential {
            weight: recombined_weight,
            is_inhibitory: recombined_inhibitory_state,
        }
    }
}

impl Reaction<SimpleNeuron> for SimpleDendriteActivationPotential {
    fn get_educt_number(&self) -> usize {
        1
    }

    fn get_product_number(&self) -> usize {
        1
    }

    fn random() -> Self {
        SimpleDendriteActivationPotential {
            weight: thread_rng().gen(),
            is_inhibitory: thread_rng().gen(),
        }
    }

    fn react(&self, educts: &[&SimpleNeuron], _reaction_time: Iteration) -> Vec<SimpleNeuron> {
        educts
            .iter()
            .map(|educt| {
                if self.is_inhibitory() {
                    educt.with_new_current_potential(educt.current_potential() - self.weight())
                } else {
                    educt.with_new_current_potential(educt.current_potential() + self.weight())
                }
            })
            .collect()
    }
}

impl Distribution<SimpleDendriteActivationPotential> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> SimpleDendriteActivationPotential {
        SimpleDendriteActivationPotential::new(rng.gen::<Nlbf64>(), rng.gen::<bool>())
    }
}

#[cfg(test)]
mod test_activation;
#[cfg(test)]
mod test_threshold;
