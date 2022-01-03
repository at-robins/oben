use serde::{Deserialize, Serialize};

use crate::evolution::{
    chemistry::{Information, Input, Reaction, State},
    gene::CrossOver,
};

use super::noop::NoOpInputElement;

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// An [`Input`] for testing purposes.
pub struct TestInput {
    pub last_feedback_update: Vec<Option<TestInformation>>,
    pub last_set_input: NoOpInputElement,
}

impl Input<NoOpInputElement, TestInformation> for TestInput {
    fn set_input(&mut self, input: NoOpInputElement) {
        self.last_set_input = input;
    }

    fn handle_feedback_substrate_changes(&mut self, changes: Vec<Option<TestInformation>>) -> bool {
        self.last_feedback_update = changes;
        true
    }

    fn random() -> Self {
        TestInput {
            last_feedback_update: Vec::new(),
            last_set_input: (),
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
