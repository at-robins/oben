//! The `simple_neuronal_mutation` module contains functions used to mutate simple neuronal networks.

use std::fmt::Debug;

use crate::evolution::binary::{as_f64, f64_to_binary, flip_random_bit};
use crate::evolution::chemistry::{Input, Output};
use crate::evolution::gene::{GeneSubstrate, Genome};
use crate::evolution::helper::Nlbf64;
use rand::{thread_rng, Rng};
use serde::de::DeserializeOwned;
use serde::Serialize;

use super::super::chemistry::{Reaction, State};
use super::super::gene::{GenomicCatalyticCentre, GenomicReceptor};
use super::{SimpleDendriteActivationPotential, SimpleDendriteThreshold, SimpleNeuron};

pub fn add_dendrite<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
        vec![source_neuron],
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
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    *neuron = SimpleNeuron::new(
        Nlbf64::flip_random_bit(neuron.base_potential()),
        neuron.potential_halflife_time(),
    );
    Some(mutated_genome)
}

pub fn mutate_neuron_potential_halflife_time<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    let mutated_halflife_time =
        as_f64(&flip_random_bit(&f64_to_binary(neuron.potential_halflife_time())));
    *neuron = SimpleNeuron::new(neuron.base_potential(), mutated_halflife_time);
    Some(mutated_genome)
}

pub fn mutate_dendrite_threshold<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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

            let educts = vec![*dendrit
                .educts()
                .get(0)
                .expect("A dendrite must have a source.")];
            let products = vec![target_neuron];
            let activation = dendrit.reaction().clone();
            *dendrit = GenomicCatalyticCentre::new(educts, products, activation);
            Some(mutated_genome)
        })
}

pub fn mutate_dendrite_source<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    let source_neuron = mutated_genome
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

            let products = vec![*dendrit
                .products()
                .get(0)
                .expect("A dendrite must have a target.")];
            let educts = vec![source_neuron];
            let activation = dendrit.reaction().clone();
            *dendrit = GenomicCatalyticCentre::new(educts, products, activation);
            Some(mutated_genome)
        })
}

pub fn add_neuron<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    if mutated_genome
        .add_substrate_to_gene(gene, SimpleNeuron::random())
        .is_some()
    {
        Some(mutated_genome)
    } else {
        None
    }
}

pub fn remove_neuron<
    InputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    let substrates_in_gene = mutated_genome.get_gene(gene).number_of_substrates().get();
    let substrate = mutated_genome.get_gene(gene).get_random_substrate();
    if substrates_in_gene > 1 {
        mutated_genome.remove_substrate(GeneSubstrate::new(gene, substrate));
        Some(mutated_genome)
    } else {
        None
    }
}

pub fn mutate_associate_input<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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

pub fn mutate_disociate_input<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    let number_of_input_substrates = mutated_genome.input().number_of_input_substrates();
    if number_of_input_substrates > 0 {
        mutated_genome
            .input_mut()
            .set_input_substrate(thread_rng().gen_range(0..number_of_input_substrates), None);
        Some(mutated_genome)
    } else {
        None
    }
}

pub fn mutate_associate_output<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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

pub fn mutation_associate_finish_substrate<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
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
    let new_finish_substrate = Some(mutated_genome.random_gene_substrate());
    mutated_genome
        .output_mut()
        .set_finish_substrate(new_finish_substrate);
    Some(mutated_genome)
}

// #[cfg(test)]
// mod tests;
