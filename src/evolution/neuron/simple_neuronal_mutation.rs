//! The `simple_neuronal_mutation` module contains functions used to mutate simple neuronal networks.

use std::fmt::Debug;

use crate::evolution::chemistry::{Input, Output};
use crate::evolution::gene::{CrossOver, Genome};
use crate::evolution::helper::Nlbf64;
use serde::de::DeserializeOwned;
use serde::Serialize;

use super::super::chemistry::{Reaction, State};
use super::super::gene::{GenomicCatalyticCentre, GenomicReceptor};
use super::{SimpleDendriteActivationPotential, SimpleDendriteThreshold, SimpleNeuron};

pub fn add_dendrite<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let gene = mutated_genome.get_random_gene();
    let source_neuron = mutated_genome.get_gene(gene).get_random_substrate();
    let target_neuron = mutated_genome.get_gene(gene).get_random_substrate();
    let dendrite_activation = GenomicCatalyticCentre::new(
        vec![source_neuron, target_neuron],
        vec![target_neuron],
        SimpleDendriteActivationPotential::random(),
    );
    let dendrite_threshold = GenomicReceptor::new(
        vec![source_neuron],
        vec![source_neuron],
        SimpleDendriteThreshold::random(),
        dendrite_activation,
    );
    mutated_genome
        .get_gene_mut(gene)
        .add_receptor(dendrite_threshold);
    Some(mutated_genome)
}

pub fn remove_dendrite<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let random_gene_index = genome.get_random_gene();
    let mut mutated_genome = genome.duplicate();
    mutated_genome
        .get_gene(random_gene_index)
        .get_random_receptor()
        .and_then(|random_receptor| {
            mutated_genome
                .get_gene_mut(random_gene_index)
                .remove_receptor(random_receptor);
            Some(mutated_genome)
        })
}

pub fn mutate_neuron_value<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let neuron_position = mutated_genome.random_gene_substrate();
    let neuron = mutated_genome.get_substrate_mut(neuron_position).unwrap();
    *neuron = SimpleNeuron::new(Nlbf64::flip_random_bit(neuron.base_potential()));
    Some(mutated_genome)
}

pub fn mutate_dendrite_threshold<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let neuron_position = mutated_genome.random_gene_substrate();
    mutated_genome
        .get_gene(neuron_position.gene())
        .get_random_receptor()
        .and_then(|dendrit_index| {
            let dendrit = mutated_genome
                .get_gene_mut(neuron_position.gene())
                .receptor_mut(dendrit_index)
                .unwrap();
            let threshold =
                SimpleDendriteThreshold::new(Nlbf64::flip_random_bit(dendrit.state().threshold()));
            let substrates = dendrit.substrates().clone();
            dendrit.replace_state(threshold, substrates);
            Some(mutated_genome)
        })
}

pub fn mutate_dendrite_weight<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let neuron_position = mutated_genome.random_gene_substrate();
    mutated_genome
        .get_gene(neuron_position.gene())
        .get_random_receptor()
        .and_then(|dendrit_index| {
            let dendrit = mutated_genome
                .get_gene_mut(neuron_position.gene())
                .receptor_mut(dendrit_index)
                .unwrap()
                .enzyme_mut();
            let activation = SimpleDendriteActivationPotential::new(
                Nlbf64::flip_random_bit(dendrit.reaction().weight()),
                dendrit.reaction().is_inhibitory(),
            );
            let educts = dendrit.educts().clone();
            let products = dendrit.products().clone();
            *dendrit = GenomicCatalyticCentre::new(educts, products, activation);
            Some(mutated_genome)
        })
}

pub fn mutate_dendrite_inhibitory_state<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let neuron_position = mutated_genome.random_gene_substrate();
    mutated_genome
        .get_gene(neuron_position.gene())
        .get_random_receptor()
        .and_then(|dendrit_index| {
            let dendrit = mutated_genome
                .get_gene_mut(neuron_position.gene())
                .receptor_mut(dendrit_index)
                .unwrap()
                .enzyme_mut();
            let activation = SimpleDendriteActivationPotential::new(
                dendrit.reaction().weight(),
                !dendrit.reaction().is_inhibitory(),
            );
            let educts = dendrit.educts().clone();
            let products = dendrit.products().clone();
            *dendrit = GenomicCatalyticCentre::new(educts, products, activation);
            Some(mutated_genome)
        })
}

pub fn mutate_dendrite_target<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let neuron_position = mutated_genome.random_gene_substrate();
    let target_neuron = mutated_genome
        .get_gene(neuron_position.gene())
        .get_random_substrate();
    mutated_genome
        .get_gene(neuron_position.gene())
        .get_random_receptor()
        .and_then(|dendrit_index| {
            let dendrit = mutated_genome
                .get_gene_mut(neuron_position.gene())
                .receptor_mut(dendrit_index)
                .unwrap()
                .enzyme_mut();

            let educts = vec![
                *dendrit
                    .educts()
                    .get(0)
                    .expect("A dendrite must have a source and a target."),
                target_neuron,
            ];
            let products = vec![target_neuron];
            let activation = dendrit.reaction().clone();
            *dendrit = GenomicCatalyticCentre::new(educts, products, activation);
            Some(mutated_genome)
        })
}

// #[cfg(test)]
// mod tests;
