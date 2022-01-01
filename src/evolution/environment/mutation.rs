//! The `mutation` module contains structures relevant for mutating genomes.

use serde::{de::DeserializeOwned, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State},
    gene::{CrossOver, Genome},
    helper::Nlbf64,
};

///
pub struct MutationList<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    mutations:
        Vec<Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>,
}

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
    pub fn new<N: Into<Nlbf64>>(
        mutation_chance: N,
        mutation: Box<
            dyn Fn(
                    &Genome<
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                    >,
                ) -> Option<
                    Genome<
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                    >,
                > + Send
                + Sync
                + 'static,
        >,
    ) -> Mutation<ReactionType, StateType, InformationType, InputElementType, InputSensorType> {
        Mutation {
            chance: mutation_chance.into(),
            mutation,
        }
    }
}
