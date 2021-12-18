//! The `neuron` module contains consturcts for working with neuronal networks.

pub use simple_neuron::SimpleNeuron;
pub use simple_dendrite::{SimpleDendriteActivationPotential, SimpleDendriteThreshold};

mod simple_dendrite;
mod simple_neuron;
//mod simple_neuronal_mutation;
