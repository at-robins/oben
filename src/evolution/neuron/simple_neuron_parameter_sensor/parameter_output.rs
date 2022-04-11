use rand::{thread_rng, Rng};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::evolution::{
    chemistry::{Input, Output},
    gene::{CrossOver, Genome},
    helper::Nlbf64,
    neuron::{SimpleDendriteActivationPotential, SimpleDendriteThreshold, SimpleNeuron},
};

#[derive(Debug, PartialEq, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleNeuronParameterOutputSensor {}

impl Output<Vec<Nlbf64>, SimpleNeuron> for SimpleNeuronParameterOutputSensor {
    fn random() -> Self {
        SimpleNeuronParameterOutputSensor {}
    }

    fn get_output(&self, information: Vec<Option<SimpleNeuron>>) -> Vec<Nlbf64> {
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
}

impl CrossOver for SimpleNeuronParameterOutputSensor {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        SimpleNeuronParameterOutputSensor {}
    }
}

pub fn mutation_associate_output<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        Vec<Nlbf64>,
        SimpleNeuronParameterOutputSensor,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        Vec<Nlbf64>,
        SimpleNeuronParameterOutputSensor,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let new_output_substrate = Some(mutated_genome.random_gene_substrate());
    let number_of_output_substrates = mutated_genome.output().number_of_output_substrates();
    if number_of_output_substrates > 0 {
        mutated_genome.output_mut().set_output_substrate(
            thread_rng().gen_range(0..number_of_output_substrates),
            new_output_substrate,
        );
        Some(mutated_genome)
    } else {
        None
    }
}

pub fn mutation_disociate_output<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        Vec<Nlbf64>,
        SimpleNeuronParameterOutputSensor,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        Vec<Nlbf64>,
        SimpleNeuronParameterOutputSensor,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let number_of_output_substrates = mutated_genome.output().number_of_output_substrates();
    if number_of_output_substrates > 0 {
        mutated_genome
            .output_mut()
            .set_output_substrate(thread_rng().gen_range(0..number_of_output_substrates), None);
        Some(mutated_genome)
    } else {
        None
    }
}
