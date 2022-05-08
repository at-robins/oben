use std::sync::Arc;

use crate::evolution::helper::Iteration;

use super::neuron::Neuron;

pub trait Interactor {
    fn number_of_input_neurons(&self) -> usize;

    fn number_of_output_neurons(&self) -> usize;

    fn set_input_neurons(&mut self, neurons: Vec<Arc<Neuron>>);

    fn set_output_neurons(&mut self, neurons: Vec<Arc<Neuron>>);

    fn update_iteration(&mut self, time: Iteration);

    fn request_evaluation(&self) -> bool;
}
