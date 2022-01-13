//! The `sensor` module contains the executive part of the evolutionary network
//! that connects to external information.

use std::{
    cell::RefCell,
    marker::PhantomData,
    rc::{Rc, Weak},
};

use serde::{de::DeserializeOwned, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State},
    gene::CrossOver,
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
    feedback_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
    old_feedback_substrate_values: Vec<Option<Substrate<ReactionType, StateType, InformationType>>>,
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
    /// * `feedback_substrates` - the [`Substrate`]s used to react to input information flow
    pub fn new(
        input: InputSensorType,
        input_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
        feedback_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    {
        let old_feedback_substrate_values: Vec<
            Option<Substrate<ReactionType, StateType, InformationType>>,
        > = substrates_as_owned(&feedback_substrates);
        InputSensor {
            phantom_i: PhantomData,
            input,
            input_substrates,
            feedback_substrates,
            old_feedback_substrate_values,
        }
    }

    /// Returns the input substrates.
    pub fn input_substrates(
        &self,
    ) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.input_substrates
    }

    /// Returns the feedback substrates.
    pub fn feedback_substrates(
        &self,
    ) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.feedback_substrates
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
    pub fn feedback_update(&mut self) -> bool {
        let current_feedback_values = substrates_as_owned(self.feedback_substrates());
        let changes: Vec<Option<InformationType>> = self
            .old_feedback_substrate_values
            .iter()
            .map(Option::as_ref)
            .zip(current_feedback_values.iter().map(Option::as_ref))
            .map(|(old_option, current_option)| match (old_option, current_option) {
                (Some(old), Some(current)) if old.value != current.value => {
                    Some(current.value.clone())
                }
                (None, Some(current)) => Some(current.value.clone()),
                _ => None,
            })
            .collect();
        let mut changes_detected = false;
        if changes.iter().any(|o| o.is_some()) {
            self.old_feedback_substrate_values = current_feedback_values;
            if let Some(changed_input_substrates) = self.input.handle_feedback_substrate_changes(changes) {
                changes_detected = true;
                self.set_input_substrates(changed_input_substrates);
        }
        }
        changes_detected
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

/// Converts a vector of references to [`Substrate`]s.
///
/// # Parameters
///
/// * `substrates` - the [`Substrate`] references
fn substrates_as_owned<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrates: &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
) -> Vec<Option<Substrate<ReactionType, StateType, InformationType>>> {
    substrates
        .iter()
        .map(|substrate| substrate.as_ref())
        .map(|substrate| substrate.map(substrate_reference_to_owned))
        .collect()
}

/// Converts a reference to a [`Substrate`].
///
/// # Parameters
///
/// * `substrate_reference` - the [`Substrate`] reference
fn substrate_reference_to_owned<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrate_reference: &Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
) -> Substrate<ReactionType, StateType, InformationType> {
    (*(substrate_reference.upgrade().unwrap())).borrow().clone()
}

#[cfg(test)]
mod test_input;
