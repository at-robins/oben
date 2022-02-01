//! The `sensor` module contains the executive part of the evolutionary network
//! that connects to external information.

use std::{
    cell::RefCell,
    marker::PhantomData,
    rc::{Rc, Weak}, collections::{HashMap},
};

use serde::{de::DeserializeOwned, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State, Output},
    gene::CrossOver, helper::Iteration,
};

use super::{substrate_reference_to_receptors, Receptor, Substrate};

/// An `InputSensor` registers and transforms inputs into internal [`Substrate`] information
/// for further processing. It can also operate on the input via feedback [`Substrate`]s.
#[derive(Debug, Clone)]
pub struct InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    phantom_i: PhantomData<InputElementType>,
    input: InputSensorType,
    input_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
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
    > InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates a new `InputSensor` to process input and transform it into internal
    /// [`Substrate`] information.
    ///
    /// # Parameters
    ///
    /// * `input` - the genetic input sensor definition
    /// * `input_substrates` - the input information processing [`Substrate`]s
    pub fn new(
        input: InputSensorType,
        input_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> Self
    {
        InputSensor {
            phantom_i: PhantomData,
            input,
            input_substrates,
        }
    }

    /// Returns the input substrates.
    pub fn input_substrates(
        &self,
    ) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.input_substrates
    }

    /// Returns the internal input implementation.
    pub fn input(&self) -> &InputSensorType {
        &self.input
    }

    /// Set the input currently registered by the `InputSensor`.
    ///
    /// # Parameters
    ///
    /// * `input` - the current input
    pub fn set_input(&mut self, input: InputElementType) {
        let input_representation = self.input.set_input(input);
        self.set_input_substrates(input_representation);
    }

    /// Updates the input [´Substrate´]s based on the supplied information.
    /// 
    /// # Parameters
    /// 
    /// * ´input_representation´ - the internal information representation to update
    /// 
    /// # Panics
    /// 
    /// If the specified representation is not of the same length as the internal input substrate vector.
    fn set_input_substrates(&mut self, input_representation: Vec<InformationType>) {
        if input_representation.len() != self.input_substrates.len() {
            panic!("{} input substrates were specified, but the input was transformed into {} substrates: {:?}", 
                self.input_substrates.len(), 
                input_representation.len(), 
                input_representation
        );
        }
        for (i, info) in input_representation.into_iter().enumerate() {
            if let Some(substrate) = (&self.input_substrates[i]).as_ref() {
                substrate.upgrade().unwrap().borrow_mut().set_value(info);
            }
        }
    }

    /// Checks all feedback [`Substrate`]s for changens and passes the changes on to the underlying input sensor representation.
    /// Updates the input [`Substrate`]s if necessary and returns if changes were detected.
    /// 
    /// # Parameters
    /// 
    /// * `changes` - the changed input [`Substrate`] values
    pub fn feedback_update(&mut self, changes: HashMap<usize, InformationType>) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        if !changes.is_empty() {
            if let Some(changed_input_substrates) = self.input.handle_feedback_substrate_changes(changes) {
                self.set_input_substrates(changed_input_substrates);
                return self.cascading_receptors();
            }
        }
        Vec::new()
    }

    /// Returns all [`Receptor`]s affected by input changes.
    pub fn cascading_receptors(
        &self,
    ) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        // This unwrap must succeed as the containing structure will always be dropped first
        // and no substrate references are leaked.
        self.input_substrates
            .iter()
            .filter(|i| i.is_some())
            .map(|i| i.as_ref())
            .map(|i| i.unwrap())
            .flat_map(substrate_reference_to_receptors)
            .collect()
    }
}

/// An `OutputSensor` registers and transforms internal [`Substrate`] information into output elements.
#[derive(Debug, Clone)]
pub struct OutputSensor<ReactionType, StateType, InformationType, OutputElementType, OutputSensorType>
{
    phantom_i: PhantomData<OutputElementType>,
    output: OutputSensorType,
    output_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
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
            + CrossOver
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    > OutputSensor<ReactionType, StateType, InformationType, OutputElementType, OutputSensorType>
{
    /// Creates a new `InputSensor` to process input and transform it into internal
    /// [`Substrate`] information.
    ///
    /// # Parameters
    ///
    /// * `input` - the genetic input sensor definition
    /// * `input_substrates` - the input information processing [`Substrate`]s
    pub fn new(
        output: OutputSensorType,
        output_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> Self
    {
        OutputSensor {
            phantom_i: PhantomData,
            output,
            output_substrates,
        }
    }

    /// Returns the output at the specified timepoint.
    /// 
    /// # Parameters
    ///
    /// * `time` - the time of the conversion
    pub fn get_output(&self, time: Iteration) -> OutputElementType {
        self.output.get_output(substrates_as_information(&self.output_substrates, time))
    }
}

/// Converts a vector of references to underlying [`Information`]s.
///
/// # Parameters
///
/// * `substrates` - the [`Substrate`] references
/// * `time` - the time of the conversion
fn substrates_as_information<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrates: &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
    time: Iteration,
) -> Vec<Option<InformationType>> {
    substrates
        .iter()
        .map(|substrate| substrate.as_ref())
        .map(|substrate| substrate.map(|s|substrate_reference_to_information(s, time)))
        .collect()
}

/// Converts a reference to a the underlying [`Information`].
///
/// # Parameters
///
/// * `substrate_reference` - the [`Substrate`] reference
/// * `time` - the time of the conversion
fn substrate_reference_to_information<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrate_reference: &Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
    time: Iteration,
) -> InformationType {
    (*(substrate_reference.upgrade().unwrap())).borrow_mut().value(time).clone()
}

#[cfg(test)]
mod test_input;
