//! The `neuron` module contains consturcts for working with neuronal networks.

pub use simple_dendrite::{SimpleDendriteActivationPotential, SimpleDendriteThreshold};
pub use simple_neuron::SimpleNeuron;
pub use simple_neuron_parameter_sensor::parameter_input::SimpleNeuronParameterInputSensor;
pub use simple_neuron_parameter_sensor::parameter_output::SimpleNeuronParameterOutputSensor;
pub use simple_neuron_parameter_sensor::random_genome;
pub use simple_neuronal_mutation::{
    add_dendrite, add_neuron, mutate_associate_input, mutate_associate_output,
    mutate_dendrite_inhibitory_state, mutate_dendrite_source, mutate_dendrite_target,
    mutate_dendrite_threshold, mutate_dendrite_weight, mutate_disociate_input,
    mutate_neuron_add_input_feedback, mutate_neuron_potential_halflife_time,
    mutate_neuron_remove_input_feedback, mutate_neuron_value, mutation_associate_finish_substrate,
    mutation_disociate_output, remove_dendrite, remove_neuron,
};

mod simple_dendrite;
mod simple_neuron;
mod simple_neuron_parameter_sensor;
mod simple_neuron_text_sensor;
mod simple_neuronal_mutation;
