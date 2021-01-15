//! The `chemistry` module contains elementary actions
//! for the evolutionary network.
// extern crate serde;

use super::gene::CrossOver;
use super::helper::Iteration;
use serde::{de::DeserializeOwned, Serialize};
use std::fmt::Debug;

/// An `Information` is a generic recombinable piece of data containing any informative value.
pub trait Information: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned {

}

/// A `State` is an elementary operation for comparing binary substrates.
pub trait State<T: Information>: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned {
    /// Compares a number of substrates for a logical property.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the substrates to perform the logical operation on
    ///
    /// # Panics
    ///
    /// If the number of substrates is not exactly equal to the required one.
    fn detect(&self, substrates: &[&T]) -> bool;

    /// Returns the number of substrates required to perform the logical comparison.
    fn get_substrate_number(&self) -> usize;

    /// Creates a random `Reaction`.
    fn random() -> Self;
}

/// A `Reaction` represents an elementary operation for modification of binary substrates.
pub trait Reaction<T: Information>: Clone + Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned {
    /// Performs the specified reaction and returns the products.
    ///
    /// # Parameters
    ///
    /// * `educts` - the educts to convert into products
    /// * `time_of_catalysis` - the timepoint at which the reaction happens
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
