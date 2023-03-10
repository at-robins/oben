use crate::evolution::binary::{as_f64, f64_to_binary, flip_random_bit};
use crate::evolution::chemistry::{Input, Output};
use crate::evolution::gene::{GeneSubstrate, Genome};
use crate::evolution::helper::Nlbf64;
use rand::{thread_rng, Rng};
use serde::de::DeserializeOwned;
use serde::Serialize;

use super::super::chemistry::{Reaction, State};
use super::MathematicalFunction;

pub fn mutate_function<
    ReactionType: Reaction<MathematicalFunction>,
    StateType: State<MathematicalFunction>,
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, MathematicalFunction>,
    OutputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    OutputSensorType: Output<OutputElementType, MathematicalFunction>,
>(
    genome: &Genome<
        ReactionType,
        StateType,
        MathematicalFunction,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
) -> Option<
    Genome<
        ReactionType,
        StateType,
        MathematicalFunction,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let function_position = mutated_genome.random_gene_substrate();
    let neuron = mutated_genome.get_substrate_mut(function_position).unwrap();
    let external_parameters: Vec<usize> = (0..3usize).collect();
    *neuron = neuron.mutate(&mut thread_rng(), external_parameters, 0.01);
    Some(mutated_genome)
}
