use serde::{Deserialize, Serialize};

use crate::evolution::{
    chemistry::Input,
    gene::CrossOver,
    helper::Nlbf64,
    neuron::{simple_neuron::POTENTIAL_HALFLIFE_TIME_MAXIMUM, SimpleNeuron},
};

/// The length of the read window of [`SimpleNeuronTextInputSensor`].
pub const READ_LENGTH: usize = 20;
/// The maximum position that can be accessed by an [`SimpleNeuronTextInputSensor`].
pub const MAX_POSITION: usize = 10_000_000 - READ_LENGTH;
/// The alphabet used by an [`SimpleNeuronTextInputSensor`].
pub const ALPHABET: [char; 73] = [
    'a', 'A', 'b', 'B', 'c', 'C', 'd', 'D', 'e', 'E', 'f', 'F', 'g', 'G', 'h', 'H', 'i', 'I', 'j',
    'J', 'k', 'K', 'l', 'L', 'm', 'M', 'n', 'N', 'o', 'O', 'p', 'P', 'q', 'Q', 'r', 'R', 's', 'S',
    't', 'T', 'u', 'U', 'v', 'V', 'w', 'W', 'x', 'X', 'y', 'Y', 'z', 'Z', 'ä', 'Ä', 'ü', 'Ü', 'ö',
    'Ö', ' ', '\t', ',', '.', '?', '!', ';', ':', '-', '\"', '\'', '\\', '/', '(', ')',
];
/// The character used to represent unknown or missing characters in a read window.
pub const UNKNOWN_CHAR: char = '�';

/// Converts a [`char`] to a [`Nlbf64`] where the numeric value represents its relative position in the alphabet.
/// If the character is not contained in the alphabet a value of `1.0` is returned.
///  
/// # Parameters
///
/// * `alphabet` - the alphabet to use for the conversion
/// * `character` - the character to convert
pub fn char_to_nlbf64(alphabet: &[char], character: char) -> Nlbf64 {
    let index = alphabet.binary_search(&character).unwrap_or(alphabet.len());
    let numeric: f64 = (index as f64) / (alphabet.len() as f64);
    numeric.into()
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Serialize, Deserialize)]
/// A text input sensor which utilises a fixed size read window to scan text, which can be controlled via a single feedback substrate.
pub struct SimpleNeuronTextInputSensor {
    input: Vec<char>,
    position: Nlbf64,
}

impl SimpleNeuronTextInputSensor {
    /// Creates a new `SimpleNeuronTextInputSensor`.
    pub fn new() -> Self {
        Self {
            input: Vec::new(),
            position: 0.0.into(),
        }
    }

    /// Returns the specified position as character index of an input string.
    ///
    /// # Parameters
    ///
    /// * `position` - the position to convert
    fn position_as_index(position: Nlbf64) -> usize {
        (position.value() * (MAX_POSITION as f64)) as usize
    }

    /// Returns the index of the start of the current read window.
    fn current_start_index(&self) -> usize {
        Self::position_as_index(self.position)
    }

    /// Returns the index of the end of the current read window.
    fn current_end_index(&self) -> usize {
        self.current_start_index() + READ_LENGTH - 1
    }

    /// Returns the currently selected read.
    pub fn current_read(&self) -> Vec<char> {
        let start_index = self.current_start_index();
        let mut end_index = self.current_end_index();
        if end_index >= self.input.len() {
            end_index = self.input.len() - 1;
        }
        self.input
            .iter()
            .enumerate()
            .filter(|(index, _)| *index >= start_index && *index <= end_index)
            .map(|(_, character)| *character)
            .collect()
    }

    /// Returns the currently selected read as neuronal representation.
    /// If the read is shorter than the input layer, it is padded with characters marked as unknown.
    fn current_read_as_neurons(&self) -> Vec<SimpleNeuron> {
        let mut current_read = self.current_read();
        // Pads the vector with unknown values if the end of the input is reached.
        if current_read.len() < READ_LENGTH {
            let difference = READ_LENGTH - current_read.len();
            let mut padding: Vec<char> = (0..difference).map(|_| UNKNOWN_CHAR).collect();
            current_read.append(&mut padding);
        }
        current_read
            .into_iter()
            .map(|character| char_to_nlbf64(&ALPHABET, character))
            .map(|value| SimpleNeuron::new(value, POTENTIAL_HALFLIFE_TIME_MAXIMUM))
            .collect()
    }
}

impl Input<String, SimpleNeuron> for SimpleNeuronTextInputSensor {
    fn set_input(&mut self, input: String) -> Vec<SimpleNeuron> {
        self.input = input.chars().collect();
        self.current_read_as_neurons()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        changes: std::collections::HashMap<usize, SimpleNeuron>,
    ) -> Option<Vec<SimpleNeuron>> {
        changes.get(&0).map(|position_neuron| {
            self.position = position_neuron.current_potential();
            self.current_read_as_neurons()
        })
    }

    fn random() -> Self {
        Self::new()
    }
}

impl CrossOver for SimpleNeuronTextInputSensor {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        Self::new()
    }
}
