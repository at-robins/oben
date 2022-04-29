use serde::{Deserialize, Serialize};

use crate::evolution::{chemistry::Output, gene::CrossOver, helper::Nlbf64, neuron::SimpleNeuron};

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleNeuronParameterOutputSensor {}

impl Output<Vec<Nlbf64>, SimpleNeuron> for SimpleNeuronParameterOutputSensor {
    fn random() -> Self {
        SimpleNeuronParameterOutputSensor {}
    }

    fn get_output(&mut self, information: Vec<Option<SimpleNeuron>>) -> Vec<Nlbf64> {
        information
            .iter()
            .map(|neuron_option| {
                neuron_option.map_or(Nlbf64::default(), |neuron| neuron.current_potential())
            })
            .collect()
    }

    fn is_finished(&self, information: SimpleNeuron) -> bool {
        information.current_potential() >= 0.9.into()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        _changes: std::collections::HashMap<usize, SimpleNeuron>,
        _current_output_information: Vec<Option<SimpleNeuron>>,
    ) -> Option<Vec<SimpleNeuron>> {
        None
    }
}

impl CrossOver for SimpleNeuronParameterOutputSensor {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        SimpleNeuronParameterOutputSensor {}
    }
}
