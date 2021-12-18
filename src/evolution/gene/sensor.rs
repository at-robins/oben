//! The `sensor` module contains input-output related genetic processing structures.
use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    marker::PhantomData,
    rc::{Rc, Weak},
};

use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State},
    helper::do_a_or_b,
    protein::{InputSensor, Substrate},
};

use super::{adjust_index, CrossOver, Gene, GeneSubstrate, Genome};

/// A `GenomicInputSensor` represents the information of an actual
/// [`InputSensor`](crate::evolution::protein::InputSensor)
/// that provides input
/// [`Substrate`](crate::evolution::protein::Substrate)s.
/// It is contained within a
/// [`Gene`](crate::evolution::gene::Gene).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenomicInputSensor<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
> {
    phantom_reaction: PhantomData<ReactionType>,
    phantom_state: PhantomData<StateType>,
    phantom_information: PhantomData<InformationType>,
    phantom_input_element: PhantomData<InputElementType>,
    input_substrates: Vec<Option<GeneSubstrate>>,
    feedback_substrates: Vec<Option<GeneSubstrate>>,
    input: InputSensorType,
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
    GenomicInputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates a new `GenomicInputSensor`.
    ///
    /// # Parameters
    ///
    /// * `input_substrates` - the input substrates if associated
    /// * `feedback_substrates` - the feedback substrates if associated
    /// * `input` - the underlying input implementation
    pub fn new(
        input_substrates: Vec<Option<GeneSubstrate>>,
        feedback_substrates: Vec<Option<GeneSubstrate>>,
        input: InputSensorType,
    ) -> Self {
        Self {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_input_element: PhantomData,
            input_substrates,
            feedback_substrates,
            input,
        }
    }

    /// Returns the input [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor.
    pub fn input_substrates(&self) -> &Vec<Option<GeneSubstrate>> {
        &self.input_substrates
    }

    /// Returns the feedback [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor.
    pub fn feedback_substrates(&self) -> &Vec<Option<GeneSubstrate>> {
        &self.feedback_substrates
    }

    /// Get the number of associated input [`Substrate`](crate::evolution::protein::Substrate)s.
    pub fn number_of_associated_inputs(&self) -> usize {
        self.input_substrates.iter().filter(|i| i.is_some()).count()
    }

    /// Adjust the input [`GeneSubstrate`] references after the [`Substrate`] of a contained
    /// [`Gene`] was removed.
    ///
    /// # Parameters
    ///
    /// * `removed_substrate` - an index based pointer to the removed [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn adjust_after_gene_substrate_removal(&mut self, removed_substrate: GeneSubstrate) {
        adjust_substrates(&mut self.input_substrates, removed_substrate);
        adjust_substrates(&mut self.feedback_substrates, removed_substrate);
    }

    /// Validates the internal
    /// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)
    /// of the references after a recombination
    /// event and removes invalid ones.
    ///
    /// # Parameters
    ///
    /// * `genes` - the [`Gene`](crate::evolution::gene::Gene)s
    /// of the [`Genome`](crate::evolution::gene::Genome)
    pub fn validate(&mut self, genes: &Vec<Gene<ReactionType, StateType, InformationType>>) {
        validate_substrates::<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >(&mut self.input_substrates, genes);
        validate_substrates::<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >(&mut self.feedback_substrates, genes);
    }

    /// Translates the `GenomicInputSensor` into an [`InputSensor`](crate::evolution::protein::InputSensor).
    ///
    /// # Parameters
    /// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`](crate::evolution::protein::Substrate)s
    /// of the containing [`Genome`](crate::evolution::gene::Genome)
    ///
    /// # Panics
    ///
    /// If the `substrate_lookup` map does not contain one of the requested [`Substrate`](crate::evolution::protein::Substrate)s.
    pub fn translate(
        &self,
        substrate_lookup: &HashMap<
            GeneSubstrate,
            Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
        >,
    ) -> InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    {
        let input_substrates = translate_substrates(&self.input_substrates, substrate_lookup);
        let feedback_substrates = translate_substrates(&self.feedback_substrates, substrate_lookup);
        let input = self.input.clone();
        InputSensor::new(input, input_substrates, feedback_substrates)
    }

    /// Removes all internal
    /// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)s
    /// referencing the [`Gene`](crate::evolution::gene::Gene)
    /// with the specified index.
    ///
    /// # Parameters
    ///
    /// * `gene_index` - the index of the [`Gene`](crate::evolution::gene::Gene)
    /// to remove associations to
    pub fn remove_associations_with_gene(&mut self, gene_index: usize) {
        remove_associations_with_gene_from_substrates(&mut self.input_substrates, gene_index);
        remove_associations_with_gene_from_substrates(&mut self.feedback_substrates, gene_index);
    }
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
        InputSensorType: Input<InputElementType, InformationType> + Default,
    > Default
    for GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    >
{
    fn default() -> Self {
        GenomicInputSensor {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_input_element: PhantomData,
            input_substrates: Vec::default(),
            feedback_substrates: Vec::default(),
            input: InputSensorType::default(),
        }
    }
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
    > CrossOver
    for GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    >
{
    fn is_similar(&self, other: &Self) -> bool {
        self.input.is_similar(&other.input)
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let input_substrates = self.input_substrates.cross_over(&other.input_substrates);
            let feedback_substrates = self
                .feedback_substrates
                .cross_over(&other.feedback_substrates);
            let input = self.input.cross_over(&other.input);
            GenomicInputSensor {
                phantom_reaction: PhantomData,
                phantom_state: PhantomData,
                phantom_information: PhantomData,
                phantom_input_element: PhantomData,
                input_substrates,
                feedback_substrates,
                input,
            }
        } else {
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

/// Removes all specified
/// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)s
/// referencing the [`Gene`](crate::evolution::gene::Gene)
/// with the specified index.
///
/// # Parameters
///
/// * `substrates` - the list of substrates to remove from
/// * `gene_index` - the index of the [`Gene`](crate::evolution::gene::Gene)
/// to remove associations to
fn remove_associations_with_gene_from_substrates(
    substrates: &mut Vec<Option<GeneSubstrate>>,
    gene_index: usize,
) {
    for substrate_option in substrates {
        if let Some(substrate) = substrate_option {
            if substrate.is_gene(gene_index) {
                *substrate_option = None;
            }
        }
    }
}

/// Translates the [`Substrate`](crate::evolution::protein::Substrate)s of a sensor
///
/// # Parameters
/// * `substrates` - the [`Substrate`](crate::evolution::protein::Substrate)s to transcribe
/// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`](crate::evolution::protein::Substrate)s
/// of the containing [`Genome`](crate::evolution::gene::Genome)
///
/// # Panics
///
/// If the `substrate_lookup` map does not contain one of the requested [`Substrate`](crate::evolution::protein::Substrate)s.
fn translate_substrates<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrates: &Vec<Option<GeneSubstrate>>,
    substrate_lookup: &HashMap<
        GeneSubstrate,
        Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
    >,
) -> Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
    substrates
        .iter()
        .map(Option::as_ref)
        .map(|gene_substrate_option| {
            gene_substrate_option.map(|gene_substrate| {
                substrate_lookup.get(gene_substrate).expect(&format!(
                    "The substrate lookup map did not contain {:?}.",
                    gene_substrate
                ))
            })
        })
        .map(|strong_option| strong_option.map(|strong| Rc::downgrade(strong)))
        .collect()
}

/// Validates the specified
/// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)
/// of the references after a recombination
/// event and removes invalid ones.
///
/// # Parameters
///
/// * `substrates` - the substrates to validate
/// * `genes` - the [`Gene`](crate::evolution::gene::Gene)s
/// of the [`Genome`](crate::evolution::gene::Genome)
fn validate_substrates<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, InformationType>,
>(
    substrates: &mut Vec<Option<GeneSubstrate>>,
    genes: &Vec<Gene<ReactionType, StateType, InformationType>>,
) {
    let mut duplicate_checker = HashSet::new();
    for substrate_option in substrates {
        if let Some(substrate) = substrate_option {
            // Remove the association if it points to a duplicate or an invalid substrate.
            if !duplicate_checker.insert(*substrate)
                || !Genome::<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                >::has_substrate(genes, substrate)
            {
                *substrate_option = None;
            }
            // Otherwise, leave it untouched.
        }
    }
}

/// Adjust the pointers and indices of referenced
/// [`Substrate`](crate::evolution::protein::Substrate)s
/// after removal.
///
/// # Parameters
///
/// * `substrates` - the
///   [`Substrate`](crate::evolution::protein::Substrate)s
///   to adjust
/// * `removed_substrate` - the removed [`Substrate`](crate::evolution::protein::Substrate)
fn adjust_substrates(
    substrates: &mut Vec<Option<GeneSubstrate>>,
    removed_substrate: GeneSubstrate,
) {
    for substrate in substrates {
        *substrate = adjust_substrate(*substrate, removed_substrate);
    }
}

/// Adjust a single potential
/// [`Substrate`](crate::evolution::protein::Substrate)
/// after removal.
///
/// # Parameters
///
/// * `substrate_option` - the `Optional` containing a potential
///   [`Substrate`](crate::evolution::protein::Substrate)
/// * `removed_substrate` - the removed [`Substrate`](crate::evolution::protein::Substrate)
fn adjust_substrate(
    substrate_option: Option<GeneSubstrate>,
    removed_substrate: GeneSubstrate,
) -> Option<GeneSubstrate> {
    substrate_option.and_then(|substrate| {
        if substrate == removed_substrate {
            None
        } else if substrate.is_gene(removed_substrate.gene()) {
            Some(GeneSubstrate::new(
                substrate.gene(),
                adjust_index(substrate.substrate(), removed_substrate.substrate()),
            ))
        } else {
            Some(substrate)
        }
    })
}
