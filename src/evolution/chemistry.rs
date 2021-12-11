//! The `chemistry` module contains elementary actions
//! for the evolutionary network.
// extern crate serde;

use super::gene::CrossOver;
use super::helper::Iteration;
use serde::{de::DeserializeOwned, Serialize};
use std::fmt::Debug;

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
    /// as [`Iteration`](oben::evolution::helper::Iteration)
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
    /// as [`Iteration`](oben::evolution::helper::Iteration)
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

/// An `Input` represents an element that can convert external stimuli to substrates.
pub trait Input<
    InputElement: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InformationType: Information,
>: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned
{
    fn set_input(&mut self, input: InputElement);

    /// Returns the number of educts required to perform the reaction.
    fn get_input_substrates(&self) -> usize;

    /// Handles changes in the feedback substrates and returns if those changes
    /// entail a change in the input sensor.
    ///
    /// # Parameters
    ///
    /// *`changes` - the changed feedback substrate values if a change occurred
    fn handle_feedback_substrate_changes(&mut self, changes: Vec<Option<InformationType>>) -> bool;

    /// Creates an random `InputSensor`.
    fn random() -> Self;
}