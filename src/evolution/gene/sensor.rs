//! The `sensor` module contains input-output related genetic processing structures.
use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    marker::PhantomData,
    rc::{Rc, Weak},
};

use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Output, Reaction, State},
    helper::{a_or_b, do_a_or_b},
    protein::{InputSensor, OutputSensor, Substrate, SubstrateType},
};

use super::{adjust_index, has_substrate, CrossOver, Gene, GeneSubstrate};

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
    feedback_substrates: HashMap<usize, GeneSubstrate>,
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
        feedback_substrates: HashMap<usize, GeneSubstrate>,
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

    /// Returns the number of input [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor.
    pub fn number_of_input_substrates(&self) -> usize {
        self.input_substrates.len()
    }

    /// Sets the input [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor
    /// and returns the old value if any.
    /// 
    /// # Panics
    /// 
    /// If the specified index is out of bounds.
    pub fn set_input_substrate(&mut self, index: usize, substrate: Option<GeneSubstrate>) -> Option<GeneSubstrate> {
        let old_substrate = self.input_substrates[index];
        self.input_substrates[index] = substrate;
        old_substrate
    }

    /// Returns the feedback [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor.
    pub fn feedback_substrates(&self) -> &HashMap<usize, GeneSubstrate> {
        &self.feedback_substrates
    }

    /// Returns the underlying [`Input`](crate::evolution::chemistry::Input) implementation.
    pub fn input(&self) -> &InputSensorType {
        &self.input
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
        self.feedback_substrates = self
            .feedback_substrates
            .iter()
            .map(|(association, substrate)| {
                (association, adjust_substrate(Some(*substrate), removed_substrate))
            })
            .filter_map(|(association, substrate_option)| {
                if let Some(substrate) = substrate_option {
                    Some((*association, substrate))
                } else {
                    None
                }
            })
            .collect();
    }

    /// Adjusts all internal
    /// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)s
    /// after removal of the [`Gene`](crate::evolution::gene::Gene)
    /// with the specified index.
    ///
    /// # Parameters
    ///
    /// * `gene_index` - the index of the removed [`Gene`](crate::evolution::gene::Gene)
    pub fn adjust_after_gene_removal(&mut self, gene_index: usize) {
        adjust_after_gene_removal_from_substrates(&mut self.input_substrates, gene_index);
        self.feedback_substrates = self
            .feedback_substrates
            .iter()
            .filter_map(|(association, substrate)| {
                if let Some(new_substrate) =
                    adjust_single_substrate_after_gene_removal_from_substrates(
                        *substrate, gene_index,
                    )
                {
                    Some((*association, new_substrate))
                } else {
                    None
                }
            })
            .collect();
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
        >(&mut self.input_substrates, genes);
        self.feedback_substrates
            .retain(|_, substrate| has_substrate(genes, substrate));
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
        // Transforms the input feedback associations into a form that is translatable.
        let transformed_input_associations: HashMap<GeneSubstrate, Vec<usize>> = self
            .feedback_substrates()
            .iter()
            .fold(HashMap::new(), |mut map, (association, substrate)| {
                map.entry(*substrate)
                    .or_insert(Vec::new())
                    .push(*association);
                map
            });
        // Sets the feedback associations.
        set_input_substrate_type(transformed_input_associations, substrate_lookup);
        let input = self.input.clone();
        InputSensor::new(input, input_substrates)
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
            feedback_substrates: HashMap::default(),
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
            let associations: HashSet<usize> = self
                .feedback_substrates()
                .keys()
                .chain(other.feedback_substrates().keys())
                .map(|association| *association)
                .collect();
            let feedback_substrates: HashMap<usize, GeneSubstrate> = associations
                .into_iter()
                .flat_map(|association| {
                    a_or_b(
                        self.feedback_substrates().get(&association),
                        other.feedback_substrates().get(&association),
                    )
                    .map(|substrate| (association, *substrate))
                })
                .collect();
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

/// A `GenomicOutputSensor` represents the information of an actual
/// [`OutputSensor`](crate::evolution::protein::OutputSensor)
/// that provides output
/// [`Substrate`](crate::evolution::protein::Substrate)s.
/// It is contained within a
/// [`Gene`](crate::evolution::gene::Gene).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenomicOutputSensor<
    ReactionType,
    StateType,
    InformationType,
    OutputElementType,
    OutputSensorType,
> {
    phantom_reaction: PhantomData<ReactionType>,
    phantom_state: PhantomData<StateType>,
    phantom_information: PhantomData<InformationType>,
    phantom_output_element: PhantomData<OutputElementType>,
    output_substrates: Vec<Option<GeneSubstrate>>,
    finish_substrate: Option<GeneSubstrate>,
    output: OutputSensorType,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    GenomicOutputSensor<
        ReactionType,
        StateType,
        InformationType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Creates a new `GenomicOutputSensor`.
    ///
    /// # Parameters
    ///
    /// * `output_substrates` - the output substrates if associated
    /// * `output` - the underlying output implementation
    pub fn new(
        output_substrates: Vec<Option<GeneSubstrate>>,
        finish_substrate: Option<GeneSubstrate>,
        output: OutputSensorType,
    ) -> Self {
        Self {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_output_element: PhantomData,
            output_substrates,
            finish_substrate,
            output,
        }
    }

    /// Returns the output [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor.
    pub fn output_substrates(&self) -> &Vec<Option<GeneSubstrate>> {
        &self.output_substrates
    }

    /// Returns the number of output [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor.
    pub fn number_of_output_substrates(&self) -> usize {
        self.output_substrates.len()
    }

    /// Sets the output [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor
    /// and returns the old value if any.
    /// 
    /// # Parameters
    /// * `index` - the index at which the new substrate should be set
    /// * `substrate` - the new substrate to set
    /// 
    /// # Panics
    /// 
    /// If the specified index is out of bounds.
    pub fn set_output_substrate(&mut self, index: usize, substrate: Option<GeneSubstrate>) -> Option<GeneSubstrate> {
        let old_substrate = self.output_substrates[index];
        self.output_substrates[index] = substrate;
        old_substrate
    }

    /// Sets the finish [`Substrate`](crate::evolution::protein::Substrate)s referenced by this sensor
    /// and returns the old value if any.
    /// 
    /// # Parameters
    /// 
    /// * `substrate` - the finish substrate to set
    pub fn set_finish_substrate(&mut self, substrate: Option<GeneSubstrate>) -> Option<GeneSubstrate> {
        let old_substrate = self.finish_substrate;
        self.finish_substrate = substrate;
        old_substrate
    }

    /// Returns the [`Substrate`](crate::evolution::protein::Substrate) signaling the termination of the network
    /// referenced by this sensor.
    pub fn finish_substrate(&self) -> &Option<GeneSubstrate> {
        &self.finish_substrate
    }

    /// Returns the underlying [`Output`](crate::evolution::chemistry::Output) implementation.
    pub fn output(&self) -> &OutputSensorType {
        &self.output
    }

    /// Get the number of associated output [`Substrate`](crate::evolution::protein::Substrate)s.
    pub fn number_of_associated_outputs(&self) -> usize {
        self.output_substrates
            .iter()
            .filter(|i| i.is_some())
            .count()
    }

    /// Adjust the output [`GeneSubstrate`] references after the [`Substrate`] of a contained
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
        adjust_substrates(&mut self.output_substrates, removed_substrate);
        self.finish_substrate = adjust_substrate(self.finish_substrate, removed_substrate);
    }

    /// Adjusts all internal
    /// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)s
    /// after removal of the [`Gene`](crate::evolution::gene::Gene)
    /// with the specified index.
    ///
    /// # Parameters
    ///
    /// * `gene_index` - the index of the removed [`Gene`](crate::evolution::gene::Gene)
    pub fn adjust_after_gene_removal(&mut self, gene_index: usize) {
        adjust_after_gene_removal_from_substrates(&mut self.output_substrates, gene_index);
        if let Some(finish_substrate) = self.finish_substrate {
            self.finish_substrate = adjust_single_substrate_after_gene_removal_from_substrates(
                finish_substrate,
                gene_index,
            );
        }
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
        >(&mut self.output_substrates, genes);
        if let Some(finish_substrate) = self.finish_substrate {
            if !has_substrate(genes, &finish_substrate) {
                self.finish_substrate = None;
            }
        }
    }

    /// Translates the `GenomicOutputSensor` into an [`OutputSensor`](crate::evolution::protein::OutputSensor).
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
    ) -> OutputSensor<ReactionType, StateType, InformationType, OutputElementType, OutputSensorType>
    {
        let output_substrates = translate_substrates(&self.output_substrates, substrate_lookup);
        // Sets the finsih associations.
        if let Some(finish_substrate) = self.finish_substrate {
            substrate_lookup
                .get(&finish_substrate)
                .expect(&format!(
                    "The substrate lookup map did not contain {:?}.",
                    finish_substrate
                ))
                .borrow_mut()
                .set_substrate_type(SubstrateType::OutputFinishSubstrate);
        }
        let output = self.output.clone();
        OutputSensor::new(output, output_substrates)
    }
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType> + Default,
    > Default
    for GenomicOutputSensor<
        ReactionType,
        StateType,
        InformationType,
        OutputElementType,
        OutputSensorType,
    >
{
    fn default() -> Self {
        GenomicOutputSensor {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_output_element: PhantomData,
            output_substrates: Vec::default(),
            finish_substrate: None,
            output: OutputSensorType::default(),
        }
    }
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    > CrossOver
    for GenomicOutputSensor<
        ReactionType,
        StateType,
        InformationType,
        OutputElementType,
        OutputSensorType,
    >
{
    fn is_similar(&self, other: &Self) -> bool {
        self.output.is_similar(&other.output)
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let output_substrates = self.output_substrates.cross_over(&other.output_substrates);
            let finish_substrate = a_or_b(self.finish_substrate, other.finish_substrate);
            let output = self.output.cross_over(&other.output);
            GenomicOutputSensor {
                phantom_reaction: PhantomData,
                phantom_state: PhantomData,
                phantom_information: PhantomData,
                phantom_output_element: PhantomData,
                output_substrates,
                finish_substrate,
                output,
            }
        } else {
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

/// Adjusts all specified
/// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)s
/// after removal of the [`Gene`](crate::evolution::gene::Gene)
/// with the specified index and removes invalid ones.
///
/// # Parameters
///
/// * `substrates` - the list of substrates to adjust
/// * `gene_index` - the index of the removed [`Gene`](crate::evolution::gene::Gene)
fn adjust_after_gene_removal_from_substrates(
    substrates: &mut Vec<Option<GeneSubstrate>>,
    gene_index: usize,
) {
    for substrate_option in substrates {
        if let Some(substrate) = substrate_option {
            *substrate_option =
                adjust_single_substrate_after_gene_removal_from_substrates(*substrate, gene_index);
        }
    }
}

/// Adjusts a single specified
/// [`GeneSubstrate`](crate::evolution::gene::GeneSubstrate)s
/// after removal of the [`Gene`](crate::evolution::gene::Gene)
/// with the specified index if still valid.
///
/// # Parameters
///
/// * `substrate` - the list of substrates to adjust
/// * `gene_index` - the index of the removed [`Gene`](crate::evolution::gene::Gene)
fn adjust_single_substrate_after_gene_removal_from_substrates(
    substrate: GeneSubstrate,
    gene_index: usize,
) -> Option<GeneSubstrate> {
    if substrate.is_gene(gene_index) {
        // Remove the reference if it pointed to the removed gene.
        None
    } else if substrate.gene() > gene_index {
        // Adjust the gene index of the reference pointer if relevant.
        Some(GeneSubstrate::new(substrate.gene() - 1, substrate.substrate()))
    } else {
        // Otherwise return the unmodified substrate reference.
        Some(substrate)
    }
}

/// Translates the [`Substrate`](crate::evolution::protein::Substrate)s of a sensor
///
/// # Parameters
///
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

/// Sets the specified input associations for the substrate lookup map.
///
/// # Parameters
///
/// * `input_feedback_associations` - the [`Substrate`](crate::evolution::protein::Substrate)s associated with specific input feedback functions
/// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`](crate::evolution::protein::Substrate)s
/// of the containing [`Genome`](crate::evolution::gene::Genome)
///
/// # Panics
///
/// If the `substrate_lookup` map does not contain one of the requested [`Substrate`](crate::evolution::protein::Substrate)s.
fn set_input_substrate_type<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    input_feedback_associations: HashMap<GeneSubstrate, Vec<usize>>,
    substrate_lookup: &HashMap<
        GeneSubstrate,
        Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
    >,
) {
    for (gene_substrate, associations) in input_feedback_associations.into_iter() {
        let substrate = substrate_lookup
            .get(&gene_substrate)
            .expect(&format!("The substrate lookup map did not contain {:?}.", gene_substrate));
        substrate
            .borrow_mut()
            .set_substrate_type(SubstrateType::InputFeedbackSubstrate(associations));
    }
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
>(
    substrates: &mut Vec<Option<GeneSubstrate>>,
    genes: &Vec<Gene<ReactionType, StateType, InformationType>>,
) {
    let mut duplicate_checker = HashSet::new();
    for substrate_option in substrates {
        if let Some(substrate) = substrate_option {
            // Remove the association if it points to a duplicate or an invalid substrate.
            if !duplicate_checker.insert(*substrate) || !has_substrate(genes, substrate) {
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

#[cfg(test)]
mod test_input;
#[cfg(test)]
mod test_output;
