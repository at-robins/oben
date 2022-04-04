//! The `neuron` module contains consturcts for working with neuronal networks.

pub use simple_dendrite::{SimpleDendriteActivationPotential, SimpleDendriteThreshold};
pub use simple_neuron::SimpleNeuron;
pub use simple_neuronal_mutation::{
    add_dendrite, add_neuron, mutate_dendrite_inhibitory_state, mutate_dendrite_target,
    mutate_dendrite_threshold, mutate_dendrite_weight, mutate_neuron_value, remove_dendrite,
    remove_neuron,
};

mod simple_dendrite;
mod simple_neuron;
mod simple_neuronal_mutation;
