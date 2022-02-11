use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State, Output},
    gene::{CrossOver, Genome},
};

use super::noop::{NoOpInputElement, NoOpOutputElement};

/// A [`Genome`] for testing purposes.
pub type TestGenome = Genome<TestReaction, TestState, TestInformation, NoOpInputElement, TestInput, NoOpOutputElement, TestOutput>;

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// An [`Input`] for testing purposes.
pub struct TestInput {
    pub last_feedback_update: HashMap<usize, TestInformation>,
    pub last_set_input: NoOpInputElement,
    pub mocked_output_vector: Vec<TestInformation>,
}

impl Input<NoOpInputElement, TestInformation> for TestInput {
    fn set_input(&mut self, input: NoOpInputElement) -> Vec<TestInformation> {
        self.last_set_input = input;
        self.mocked_output_vector.clone()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        changes: HashMap<usize, TestInformation>,
    ) -> std::option::Option<Vec<TestInformation>> {
        self.last_feedback_update = changes.clone();
        if !changes.is_empty() {
            Some(self.mocked_output_vector.clone())
        } else {
            None
        }
    }

    fn random() -> Self {
        TestInput {
            last_feedback_update: HashMap::new(),
            last_set_input: (),
            mocked_output_vector: Vec::new(),
        }
    }
}

impl CrossOver for TestInput {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        Self::random()
    }
}

impl Default for TestInput {
    fn default() -> Self {
        Self {
            last_feedback_update: Default::default(),
            last_set_input: Default::default(),
            mocked_output_vector: Default::default(),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// An [`Output`] for testing purposes.
pub struct TestOutput {}

impl Output<NoOpOutputElement, TestInformation> for TestOutput {
    fn get_output(&self, _information: Vec<Option<TestInformation>>) -> NoOpOutputElement {
        ()
    }

    fn is_finished(&self, _information: TestInformation) -> bool {
        false
    }

    fn random() -> Self {
        TestOutput {}
    }
}

impl CrossOver for TestOutput {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        Self::random()
    }
}

impl Default for TestOutput {
    fn default() -> Self {
        Self {}
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// An [`Information`] for testing purposes.
pub struct TestInformation {
    pub value: usize,
}

impl Information for TestInformation {
    fn update_value(&mut self, _time_passed: i32) {}
}

impl CrossOver for TestInformation {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        self.clone()
    }
}

impl Default for TestInformation {
    fn default() -> Self {
        Self {
            value: Default::default(),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A [`Reaction`] for testing purposes.
pub struct TestReaction {}

impl CrossOver for TestReaction {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        TestReaction {}
    }
}

impl Reaction<TestInformation> for TestReaction {
    fn react(
        &self,
        _educts: &[&TestInformation],
        _reaction_time: crate::evolution::helper::Iteration,
    ) -> Vec<TestInformation> {
        Vec::new()
    }

    fn get_educt_number(&self) -> usize {
        0
    }

    fn get_product_number(&self) -> usize {
        0
    }

    fn random() -> Self {
        TestReaction {}
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A [`State`] for testing purposes.
pub struct TestState {}

impl CrossOver for TestState {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        TestState {}
    }
}

impl State<TestInformation> for TestState {
    fn detect(
        &self,
        _substrates: &[&TestInformation],
        _detection_time: crate::evolution::helper::Iteration,
    ) -> bool {
        false
    }

    fn get_substrate_number(&self) -> usize {
        0
    }

    fn random() -> Self {
        TestState {}
    }
}
