use std::sync::Arc;

use crate::evolution::helper::Iteration;

use super::{neuron::Neuron, network::ErrorPropagationImpulse};

pub trait Interactor {
    fn id(&self) -> u32;

    fn number_of_input_neurons(&self) -> usize;

    fn number_of_output_neurons(&self) -> usize;

    fn set_input_neurons(&mut self, neurons: Vec<Arc<Neuron>>);

    fn set_output_neurons(&mut self, neurons: Vec<Arc<Neuron>>);

    fn input_neurons(&self) -> Vec<Arc<Neuron>>;

    fn output_neurons(&self) -> Vec<Arc<Neuron>>;

    fn update_iteration(&mut self, time: Iteration);

    fn request_evaluation(&mut self) -> bool;

    fn evalute_results(&mut self) -> Option<Vec<(Arc<Neuron>, ErrorPropagationImpulse)>>;

    fn initialise_new_evaluation(&mut self);
}
