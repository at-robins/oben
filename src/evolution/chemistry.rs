//! The `chemistry` module contains elementary actions
//! for the evolutionary network.
// extern crate serde;

use super::gene::CrossOver;
use super::helper::Iteration;
use serde::{de::DeserializeOwned, Serialize};
use std::{collections::HashMap, fmt::Debug};

/// An `Information` is a generic recombinable piece of data containing any informative value.
pub trait Information:
    Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned
{
    /// Updates the internal information value based on the time that passed since the last update.
    ///
    /// # Parameters
    ///
    /// * `time_passed` - the time passed since the last update of the value in iteration steps
    fn update_value(&mut self, time_passed: i32);
}

/// A `State` is an elementary operation for comparing substrates.
pub trait State<T: Information>:
    Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned
{
    /// Compares a number of substrates for a logical property.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the substrates to perform the logical operation on
    /// * `detection_time` - the timepoint at which the detection happens
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    ///
    /// # Panics
    ///
    /// If the number of substrates is not exactly equal to the required one.
    fn detect(&self, substrates: &[&T], detection_time: Iteration) -> bool;

    /// Returns the number of substrates required to perform the logical comparison.
    fn get_substrate_number(&self) -> usize;

    /// Creates a random `Reaction`.
    fn random() -> Self;
}

/// A `Reaction` represents an elementary operation for modification of substrates.
pub trait Reaction<T: Information>:
    Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned
{
    /// Performs the specified reaction and returns the products.
    ///
    /// # Parameters
    ///
    /// * `educts` - the educts to convert into products
    /// * `reaction_time` - the timepoint at which the reaction happens
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    ///
    /// # Panics
    ///
    /// If the number of supplied educts and created products is not
    /// exactly equal to the required one.
    fn react(&self, educts: &[&T], reaction_time: Iteration) -> Vec<T>;

    /// Returns the number of educts required to perform the reaction.
    fn get_educt_number(&self) -> usize;

    /// Returns the number of products required to perform the reaction.
    fn get_product_number(&self) -> usize;

    /// Creates a random `Reaction`.
    fn random() -> Self;
}

/// An `Input` represents an element that can convert external stimuli to
/// [`Substrate`](crate::evolution::protein::Substrate)s.
pub trait Input<
    InputElement: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InformationType: Information,
>: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned
{
    /// Sets the input currently registered by the input sensor and returns the
    /// internal corresponding [`Information`] representation.
    ///
    /// # Parameters
    ///
    /// * `input` - the input element that is processed and translated into internal
    /// [`Information`]
    fn set_input(&mut self, input: InputElement) -> Vec<InformationType>;

    /// Handles changes in the feedback substrates and returns the internal corresponding [`Information`] representation
    /// if any change occured.
    ///
    /// # Parameters
    ///
    /// *`changes` - the changed feedback
    /// [`Substrate`](crate::evolution::protein::Substrate)
    /// values if a change occurred
    fn handle_feedback_substrate_changes(
        &mut self,
        changes: HashMap<usize, InformationType>,
    ) -> Option<Vec<InformationType>>;

    /// Creates an random `InputSensor`.
    fn random() -> Self;
}

/// An `Output` represents an element that can convert internal
/// [`Substrate`](crate::evolution::protein::Substrate)s.
/// to an external representation.
pub trait Output<
    OutputElement: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InformationType: Information,
>: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned
{
    /// Returns the output element corresponding to the current
    /// output [`Information`] representation.
    /// 
    /// # Parameters
    /// 
    /// * `information` - the [`Information`] representation of the output
    fn get_output(&self, information: Vec<Option<InformationType>>) -> OutputElement;

    /// Checks if the network signaled that the result is ready.
    /// 
    /// # Parameters
    /// 
    /// * `information` - the [`Information`] representing the finished state
    fn is_finished(&self, information: Option<InformationType>) -> bool;

    /// Creates an random `OutputSensor`.
    fn random() -> Self;
}
