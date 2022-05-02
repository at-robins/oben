use std::{collections::HashMap, num::NonZeroU32};

use rand::{thread_rng, Rng};

use crate::evolution::{
    chemistry::Output,
    gene::{
        Gene, GeneSubstrate, Genome, GenomicCatalyticCentre, GenomicInputSensor,
        GenomicOutputSensor, GenomicReceptor,
    },
};

use self::audio_output_sixteen_bit::SimpleNeuronAudioSixteenOutputSensor;

use super::{
    simple_neuron_text_sensor::text_input::{SimpleNeuronTextInputSensor, READ_LENGTH},
    SimpleDendriteActivationPotential, SimpleDendriteThreshold, SimpleNeuron,
};

pub(crate) mod audio_output_sixteen_bit;

/// Creates a random genome containing a simple neural network and parameter IO.
///
/// # Parameters
///
/// * `number_of_neurons` - the number of neurons the initial network should contain
/// * `number_of_dendrites` - the number of dendrites the initial network should contain
/// * `number_of_input_parameters` - the length of the input vector
/// * `number_of_output_parameters` - the length of the output vector
pub fn random_genome_tts_sixteen(
    number_of_neurons: NonZeroU32,
    number_of_dendrites: NonZeroU32,
) -> Genome<
    SimpleDendriteActivationPotential,
    SimpleDendriteThreshold,
    SimpleNeuron,
    String,
    SimpleNeuronTextInputSensor,
    Vec<i16>,
    SimpleNeuronAudioSixteenOutputSensor,
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
    let input_substrates = (0..READ_LENGTH)
        .map(|_| Some(GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates))))
        .collect();
    let mut input_feedback = HashMap::new();
    input_feedback
        .insert(0, GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates)));
    let input = GenomicInputSensor::new(
        input_substrates,
        input_feedback,
        SimpleNeuronTextInputSensor::new(),
    );
    let output_substrates = (0..2)
        .map(|_| Some(GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates))))
        .collect();
    let mut output_feedback = HashMap::new();
    output_feedback
        .insert(0, GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates)));
    let finish_substrate =
        Some(GeneSubstrate::new(0, thread_rng().gen_range(0..number_of_substrates)));
    let output = GenomicOutputSensor::new(
        output_substrates,
        output_feedback,
        finish_substrate,
        SimpleNeuronAudioSixteenOutputSensor::random(),
    );
    Genome::new(input, output, genes)
}
