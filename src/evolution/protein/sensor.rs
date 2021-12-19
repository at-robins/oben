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
    information_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
    feedback_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
    old_feedback_substrate_values: Vec<Option<Substrate<ReactionType, StateType, InformationType>>>,
    changed: bool,
}

impl<
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
    > InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates a new `InputSensor` to process input and transform it into internal
    /// [`Substrate`] information.
    ///
    /// # Parameters
    ///
    /// * `input` - the genetic input sensor definition
    /// * `information_substrates` - the information processing [`Substrate`]s
    /// * `feedback_substrates` - the [`Substrate`]s used to react to input information flow
    pub fn new(
        input: InputSensorType,
        information_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
        feedback_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    {
        let old_feedback_substrate_values: Vec<
            Option<Substrate<ReactionType, StateType, InformationType>>,
        > = Self::substrates_as_owned(&feedback_substrates);
        InputSensor {
            phantom_i: PhantomData,
            input,
            information_substrates,
            feedback_substrates,
            old_feedback_substrate_values,
            changed: false,
        }
    }

    /// Returns the input substrates.
    pub fn input_substrates(
        &self,
    ) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.information_substrates
    }

    /// Returns the feedback substrates.
    pub fn feedback_substrates(
        &self,
    ) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.feedback_substrates
    }

    /// Set the input currently registered by the `InputSensor`.
    ///
    /// # Parameters
    ///
    /// * `input` - the current input
    pub fn set_input(&mut self, input: InputElementType) {
        self.changed = true;
        self.input.set_input(input);
    }

    /// Checks if the input was changed since the last query.
    /// The change flag is reset upon calling this method.
    pub fn was_changed(&mut self) -> bool {
        if self.changed {
            self.changed = false;
            true
        } else {
            false
        }
    }

    /// Checks all feedback [`Substrate`]s for changens and passes the changes on to the underlying input sensor representation.
    pub fn feedback_update(&mut self) {
        let current_feedback_values = self.feedback_substrates_as_owned();
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

        if changes.iter().any(|o| o.is_some()) {
            self.changed = self.input.handle_feedback_substrate_changes(changes);
            self.old_feedback_substrate_values = current_feedback_values;
        }
    }

    /// Returns all [`Receptor`]s affected by input changes.
    pub fn cascading_receptors(
        &self,
    ) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        // This unwrap must succeed as the containing structure will always be dropped first
        // and no substrate references are leaked.
        self.information_substrates
            .iter()
            .filter(|i| i.is_some())
            .map(|i| i.as_ref())
            .map(|i| i.unwrap())
            .flat_map(substrate_reference_to_receptors)
            .collect()
    }

    /// Returns the feedback [`Substrate`]s as owned structures.
    fn feedback_substrates_as_owned(
        &self,
    ) -> Vec<Option<Substrate<ReactionType, StateType, InformationType>>> {
        InputSensor::<ReactionType, StateType, InformationType, InputElementType, InputSensorType>::substrates_as_owned(&self.feedback_substrates)
    }

    /// Converts a vector of references to [`Substrate`]s.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the [`Substrate`] references
    fn substrates_as_owned(
        substrates: &Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> Vec<Option<Substrate<ReactionType, StateType, InformationType>>> {
        substrates
            .iter()
            .map(|substrate| substrate.as_ref())
            .map(|substrate| substrate.map(Self::substrate_reference_to_owned))
            .collect()
    }

    /// Converts a reference to a [`Substrate`].
    ///
    /// # Parameters
    ///
    /// * `substrate_reference` - the [`Substrate`] reference
    fn substrate_reference_to_owned(
        substrate_reference: &Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
    ) -> Substrate<ReactionType, StateType, InformationType> {
        (*(substrate_reference.upgrade().unwrap())).borrow().clone()
    }
}
