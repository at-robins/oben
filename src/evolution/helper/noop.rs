//! The `parameter` module contains aliases and helper functions for working with mathematical
//! parameters.
extern crate serde;

use crate::evolution::chemistry::Input;
use crate::evolution::gene::Gene;

use super::super::chemistry::{Information, Reaction, State};
use super::super::gene::{CrossOver, Genome};
use super::super::helper::Iteration;

/// A type alias for the representation of an empty
/// [`Genome`](crate::evolution::gene::Genome)
/// without any function.
pub type NoOpGenome =
    Genome<NoOpReaction, NoOpState, NoOpSubstrate, NoOpInputElement, NoOpInputSensor>;

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

impl Reaction<NoOpSubstrate> for NoOpReaction {
    fn react(&self, _educts: &[&NoOpSubstrate], _reaction_time: Iteration) -> Vec<NoOpSubstrate> {
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

impl State<NoOpSubstrate> for NoOpState {
    fn detect(&self, _substrates: &[&NoOpSubstrate], _detection_time: Iteration) -> bool {
        false
    }

    fn get_substrate_number(&self) -> usize {
        0
    }

    fn random() -> Self {
        ()
    }
}

impl Input<NoOpInputElement, NoOpSubstrate> for NoOpInputSensor {
    fn set_input(&mut self, _input: NoOpInputElement) {}

    fn handle_feedback_substrate_changes(&mut self, _changes: Vec<Option<NoOpSubstrate>>) -> bool {
        false
    }

    fn random() -> Self {
        ()
    }
}
