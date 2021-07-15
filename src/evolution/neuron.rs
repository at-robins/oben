//! The `neuron` module contains consturcts for working with neuronal networks.
//use super::gene::Genome;

pub use simple_neuron::SimpleNeuron;
//pub use simple_dendrite::{SimpleDendriteActivation, SimpleDendriteThreshold};

//pub type SimpleNeuronalGenome = Genome<SimpleDendriteActivation, SimpleDendriteThreshold, SimpleNeuron>;

//mod simple_dendrite;
mod simple_neuron;
//mod simple_neuronal_mutation;
