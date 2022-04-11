use rand::{thread_rng, Rng};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::evolution::{
    chemistry::{Input, Output},
    gene::{CrossOver, Genome},
    helper::Nlbf64,
    neuron::{
        simple_neuron::POTENTIAL_HALFLIFE_TIME_MAXIMUM, SimpleDendriteActivationPotential,
        SimpleDendriteThreshold, SimpleNeuron,
    },
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

pub fn mutation_associate_input<
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        Vec<Nlbf64>,
        SimpleNeuronParameterInputSensor,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        Vec<Nlbf64>,
        SimpleNeuronParameterInputSensor,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let new_input_substrate = Some(mutated_genome.random_gene_substrate());
    let number_of_input_substrates = mutated_genome.input().number_of_input_substrates();
    if number_of_input_substrates > 0 {
        mutated_genome.input_mut().set_input_substrate(
            thread_rng().gen_range(0..number_of_input_substrates),
            new_input_substrate,
        );
        Some(mutated_genome)
    } else {
        None
    }
}

pub fn mutation_disociate_input<
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        Vec<Nlbf64>,
        SimpleNeuronParameterInputSensor,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        Vec<Nlbf64>,
        SimpleNeuronParameterInputSensor,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let number_of_input_substrates = mutated_genome.input().number_of_input_substrates();
    if number_of_input_substrates > 0 {
        mutated_genome.input_mut().set_input_substrate(
            thread_rng().gen_range(0..number_of_input_substrates),
            None,
        );
        Some(mutated_genome)
    } else {
        None
    }
}
