//! The `mutation` module contains structures relevant for mutating genomes.

use std::iter;

use rand::{prelude::SliceRandom, thread_rng, Rng};
use serde::{de::DeserializeOwned, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State},
    gene::{CrossOver, Genome},
    helper::Nlbf64,
};

/// The maximum number of mutation events that might occur
/// for a single [`Mutation`] type upon [`Genome`] duplication.
const MAX_MUTATION_EVENTS: usize = 128;

#[derive(Debug, PartialEq)]
/// A list of [`Mutation`]s.
pub struct MutationCompendium<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
> {
    mutations:
        Vec<Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
    >
    MutationCompendium<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates an empty `MutationCompendium`.
    pub fn new() -> MutationCompendium<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    > {
        MutationCompendium {
            mutations: Vec::new(),
        }
    }

    /// Returns the number of different mutations contained within the compendium.
    pub fn size(&self) -> usize {
        self.mutations.len()
    }

    /// Adds a [`Mutation`] to the list.
    ///
    /// # Parameters
    ///
    /// * `mutation` - the mutation to add to the compendium
    pub fn add(
        &mut self,
        mutation: Mutation<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
    ) {
        self.mutations.push(mutation);
    }

    /// Mutates a [`Genome`] based on the [`Mutation`]s and
    /// their respective mutation frequencies and returns
    /// the mutated [`Genome`] if any [`Mutation`] occured
    /// and it is still in a valid state.
    ///
    /// # Parameters
    ///
    /// * `genome` - the [`Genome`] to mutate
    pub fn mutate(
        &self,
        genome: &Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
    ) -> Option<Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>
    {
        let mut applied_mutations: Vec<
            &Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
        > = self
            .mutations
            .iter()
            .flat_map(|mutation| iter::repeat(mutation).take(mutation.number_of_mutation_events()))
            .collect();
        // Randomise the order in which mutations are applied.
        applied_mutations.shuffle(&mut thread_rng());
        let mut mutated_genome = None;
        for (i, mutation) in applied_mutations.into_iter().enumerate() {
            if i == 0 {
                // If it is the first mutation use the original genome.
                mutated_genome = mutation.mutate(genome);
            } else if let Some(g) = mutated_genome {
                // On any subsequent mutation use the mutated genome.
                mutated_genome = mutation.mutate(&g);
            } else {
                // If any of the mutations was not successful the rest of the loop can be skipped.
                break;
            }
        }
        mutated_genome
    }
}

impl<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    From<Vec<Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>>
    for MutationCompendium<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    >
{
    fn from(
        mutations: Vec<
            Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
        >,
    ) -> Self {
        MutationCompendium { mutations }
    }
}

impl<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    From<
        MutationCompendium<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
    >
    for Vec<Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>
{
    fn from(
        compendium: MutationCompendium<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
    ) -> Self {
        compendium.mutations
    }
}

/// A `Mutation` that occurs with a specified frequency druing
/// [`Genome`] duplication.
pub struct Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType> {
    chance: Nlbf64,
    mutation: Box<
        dyn Fn(
                &Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
            ) -> Option<
                Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
            > + Send
            + Sync
            + 'static,
    >,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
    > Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates a new `Mutation` which occurs with the specified frequency.
    ///
    /// # Parameters
    ///
    /// * `mutation_chance` - the chance of the mutation happening a single time during genome duplication
    /// * `mutation` - the underlying mutation function
    pub fn new<
        N: Into<Nlbf64>,
        F: Fn(
                &Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
            ) -> Option<
                Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
            > + Send
            + Sync
            + 'static,
    >(
        mutation_chance: N,
        mutation: F,
    ) -> Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType> {
        Mutation {
            chance: mutation_chance.into(),
            mutation: Box::new(mutation),
        }
    }

    /// Returns the chance of this mutation happening during [`Genome`] duplication a single time.
    pub fn mutation_chance(&self) -> Nlbf64 {
        self.chance
    }

    /// Apply this `Mutation` to a [`Genome`]
    /// and return the mutated [`Genome`]
    /// if still in a valid state.
    ///
    /// # Parameters
    ///
    /// * `genome` - the [`Genome`] to mutate
    pub fn mutate(
        &self,
        genome: &Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
    ) -> Option<Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>
    {
        (self.mutation)(genome)
    }

    /// Simulates how often a mutation event occurs on a statistical basis.
    /// This is an internal function with a random element.
    fn number_of_mutation_events(&self) -> usize {
        if self.chance == Nlbf64::MIN {
            0
        } else if self.chance == Nlbf64::MAX {
            MAX_MUTATION_EVENTS
        } else {
            let random_chance: f64 = thread_rng().gen_range(0.0..1.0);
            // Calculate the number of mutation events per mutation
            // corresponding to the generated uniform random percentage.
            //    P("n mutations in a single genome") = "mutation rate" ^ n
            // => n = log(base: "mutation rate", value: P)
            let number_of_mutations =
                random_chance.log(self.mutation_chance().value()).floor() as usize;
            // Limits the number of mutation events per mutation to
            // prevent overflow in rare statistical cases.
            if number_of_mutations > MAX_MUTATION_EVENTS {
                MAX_MUTATION_EVENTS
            } else {
                number_of_mutations
            }
        }
    }
}

impl<ReactionType, StateType, InformationType, InputElementType, InputSensorType> std::fmt::Debug
    for Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mutation")
            .field("chance", &self.chance)
            .field("mutation", &format!("{:p}", &self.mutation))
            .finish()
    }
}

impl<ReactionType, StateType, InformationType, InputElementType, InputSensorType> PartialEq
    for Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    fn eq(&self, other: &Self) -> bool {
        self.chance == other.chance && std::ptr::eq(&self.mutation, &other.mutation)
    }
}

impl<ReactionType, StateType, InformationType, InputElementType, InputSensorType> Eq
    for Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
}

#[cfg(test)]
mod tests;
