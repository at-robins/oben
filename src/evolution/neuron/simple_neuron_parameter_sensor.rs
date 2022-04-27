use std::{collections::HashMap, num::NonZeroU32};

use rand::{thread_rng, Rng};

use crate::evolution::{
    gene::{
        Gene, GeneSubstrate, Genome, GenomicCatalyticCentre, GenomicInputSensor,
        GenomicOutputSensor, GenomicReceptor,
    },
    helper::Nlbf64,
};

use super::{
    SimpleDendriteActivationPotential, SimpleDendriteThreshold, SimpleNeuron,
    SimpleNeuronParameterInputSensor, SimpleNeuronParameterOutputSensor,
};

pub(crate) mod parameter_input;
pub(crate) mod parameter_output;

/// Creates a random genome containing a simple neural network and parameter IO.
/// 
/// # Parameters
/// 
/// * `number_of_neurons` - the number of neurons the initial network should contain
/// * `number_of_dendrites` - the number of dendrites the initial network should contain
/// * `number_of_input_parameters` - the length of the input vector
/// * `number_of_output_parameters` - the length of the output vector
pub fn random_genome(
    number_of_neurons: NonZeroU32,
    number_of_dendrites: NonZeroU32,
    number_of_input_parameters: u32,
    number_of_output_parameters: u32,
) -> Genome<
    SimpleDendriteActivationPotential,
    SimpleDendriteThreshold,
    SimpleNeuron,
    Vec<Nlbf64>,
    SimpleNeuronParameterInputSensor,
    Vec<Nlbf64>,
    SimpleNeuronParameterOutputSensor,
> {
    let substrates = (0..number_of_neurons.get())
        .map(|_| thread_rng().gen::<SimpleNeuron>())
        .collect();
    let mut genes = vec![Gene::new(substrates)];
    let number_of_substrates = genes[0].number_of_substrates().get();
    for _ in 0..number_of_dendrites.get() {
        let educts = vec![thread_rng().gen_range(0..number_of_substrates)];
        let products = vec![thread_rng().gen_range(0..number_of_substrates)];
        let activation = GenomicCatalyticCentre::new(
            educts.clone(),
            products,
            thread_rng().gen::<SimpleDendriteActivationPotential>(),
        );
        let dendrite = GenomicReceptor::new(
            educts.clone(),
            educts,
            thread_rng().gen::<SimpleDendriteThreshold>(),
            activation,
        );
        genes[0].add_receptor(dendrite);
    }
    let input_substrates = (0..number_of_input_parameters)
        .map(|_| Some(GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates))))
        .collect();
    let input = GenomicInputSensor::new(
        input_substrates,
        HashMap::new(),
        SimpleNeuronParameterInputSensor {},
    );
    let output_substrates = (0..number_of_output_parameters)
        .map(|_| Some(GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates))))
        .collect();
    let finish_substrate =
        Some(GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates)));
    let output = GenomicOutputSensor::new(
        output_substrates,
        HashMap::new(),
        finish_substrate,
        SimpleNeuronParameterOutputSensor {},
    );
    Genome::new(input, output, genes)
}
