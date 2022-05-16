//! The `parameter` module contains aliases and helper functions for working with mathematical
//! parameters.
extern crate serde;

use std::collections::HashMap;

use serde::de::DeserializeOwned;
use serde::Serialize;

use crate::evolution::chemistry::{Input, Output};
use crate::evolution::gene::Gene;

use super::super::chemistry::{Information, Reaction, State};
use super::super::gene::{CrossOver, Genome};
use super::super::helper::Iteration;

/// A type alias for the representation of an empty
/// [`Genome`](crate::evolution::gene::Genome)
/// without any function.
pub type NoOpGenome = Genome<
    NoOpReaction,
    NoOpState,
    NoOpSubstrate,
    NoOpInputElement,
    NoOpInputSensor,
    NoOpOutputElement,
    NoOpOutputSensor,
>;

/// A type alias for the representation of an empty
/// [`Gene`](crate::evolution::gene::Gene)
/// without any function.
pub type NoOpGene = Gene<NoOpReaction, NoOpState, NoOpSubstrate>;

/// A type alias for the underlying representation of a
/// [`Substrate`](crate::evolution::protein::Substrate)
/// containing no information or value.
pub type NoOpSubstrate = ();

/// A type alias for the underlying representation of a
/// [`Reaction`](crate::evolution::chemistry::Reaction)
/// having no function.
pub type NoOpReaction = ();

/// A type alias for the underlying representation of a
/// [`State`](crate::evolution::chemistry::State)
/// having no function.
pub type NoOpState = ();

/// A type alias for the underlying representation of an input
/// containing no information or value and having no function.
pub type NoOpInputElement = ();

/// A type alias for the underlying representation of an
/// [`Input`](crate::evolution::chemistry::Input)
/// containing no information or value and having no function.
pub type NoOpInputSensor = ();

/// A type alias for the underlying representation of an output
/// containing no information or value and having no function.
pub type NoOpOutputElement = ();

/// A type alias for the underlying representation of an
/// [`Output`](crate::evolution::chemistry::Output)
/// containing no information or value and having no function.
pub type NoOpOutputSensor = ();

impl Information for NoOpSubstrate {
    fn update_value(&mut self, _time_passed: i32) {}
}

impl CrossOver for NoOpSubstrate {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        ()
    }
}

impl<T> Reaction<T> for NoOpReaction
where
    T: Information,
{
    fn react(&self, _educts: &[&T], _reaction_time: Iteration) -> Vec<T> {
        Vec::new()
    }

    fn get_educt_number(&self) -> usize {
        0
    }

    fn get_product_number(&self) -> usize {
        0
    }

    fn random() -> Self {
        ()
    }
}

impl<T> State<T> for NoOpState
where
    T: Information,
{
    fn detect(&self, _substrates: &[&T], _detection_time: Iteration) -> bool {
        false
    }

    fn get_substrate_number(&self) -> usize {
        0
    }

    fn random() -> Self {
        ()
    }
}

impl<T, S> Input<T, S> for NoOpInputSensor
where
    T: Clone + Serialize + DeserializeOwned + Send + Sync + std::fmt::Debug + PartialEq,
    S: Information,
{
    fn set_input(&mut self, _input: T) -> Vec<S> {
        Vec::new()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        _changes: HashMap<usize, S>,
    ) -> std::option::Option<Vec<S>> {
        None
    }

    fn random() -> Self {
        ()
    }
}

impl<T, S> Output<T, S> for NoOpOutputSensor
where
    T: Clone + Serialize + DeserializeOwned + Send + Sync + std::fmt::Debug + PartialEq + Default,
    S: Information,
{
    fn get_output(&mut self, _information: Vec<Option<S>>) -> T {
        T::default()
    }

    fn is_finished(&self, _information: S) -> bool {
        false
    }

    fn random() -> Self {
        ()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        _changes: HashMap<usize, S>,
        _current_output_information: Vec<Option<S>>,
    ) -> Option<Vec<S>> {
        None
    }
}
