//! The `mutation` module contains structures relevant for mutating genomes.

use serde::{de::DeserializeOwned, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State},
    gene::{CrossOver, Genome},
    helper::Nlbf64,
};

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
    /// * `mutation_chance` - the chance of the mutation happening during genome duplication
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

    /// Returns the chance of this mutation happening during [`Genome`] duplication.
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
