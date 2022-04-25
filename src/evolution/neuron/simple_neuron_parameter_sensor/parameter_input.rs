use serde::{Deserialize, Serialize};

use crate::evolution::{
    chemistry::Input,
    gene::CrossOver,
    helper::Nlbf64,
    neuron::{simple_neuron::POTENTIAL_HALFLIFE_TIME_MAXIMUM, SimpleNeuron},
};

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleNeuronParameterInputSensor {}

impl Input<Vec<Nlbf64>, SimpleNeuron> for SimpleNeuronParameterInputSensor {
    fn set_input(&mut self, input: Vec<Nlbf64>) -> Vec<SimpleNeuron> {
        input
            .iter()
            .map(|param_value| SimpleNeuron::new(*param_value, POTENTIAL_HALFLIFE_TIME_MAXIMUM))
            .collect()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        _changes: std::collections::HashMap<usize, SimpleNeuron>,
    ) -> Option<Vec<SimpleNeuron>> {
        None
    }

    fn random() -> Self {
        SimpleNeuronParameterInputSensor {}
    }
}

impl CrossOver for SimpleNeuronParameterInputSensor {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        SimpleNeuronParameterInputSensor {}
    }
}
