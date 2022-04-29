use serde::{Deserialize, Serialize};

use crate::evolution::{
    chemistry::Output,
    gene::CrossOver,
    helper::Nlbf64,
    neuron::{simple_neuron::POTENTIAL_HALFLIFE_TIME_MAXIMUM, SimpleNeuron},
};

/// The length of the read window of [`SimpleNeuronAudioSixteenOutputSensor`].
pub const OUTPUT_WINDOW: usize = 100;
/// The maximum position that can be accessed by an [`SimpleNeuronAudioSixteenOutputSensor`].
pub const MAX_POSITION: usize = 10_000_000 - OUTPUT_WINDOW;

/// Converts a [`Nlbf64`] to a [`i16`] where the numeric value represents its relative position in the alphabet.
///  
/// # Parameters
///
/// * `value` - the neuron value to convert
pub fn nlbf64_to_audio_sixteen(value: Nlbf64) -> i16 {
    let range: f64 = (i16::MAX as f64) - (i16::MIN as f64);
    let converted: f64 = (i16::MIN as f64) + (value.value() * range);
    converted as i16
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Serialize, Deserialize)]
/// An sixteen bit audio output sensor which utilises a fixed size window to convert neuron values to an audio signal,
/// which can be controlled via a single feedback substrate.
pub struct SimpleNeuronAudioSixteenOutputSensor {
    audio: Vec<Nlbf64>,
    position: Nlbf64,
}

impl SimpleNeuronAudioSixteenOutputSensor {
    /// Creates a new `SimpleNeuronAudioSixteenOutputSensor`.
    pub fn new() -> Self {
        Self {
            audio: Vec::new(),
            position: 0.0.into(),
        }
    }

    /// Returns the specified position as character index of an output signal.
    ///
    /// # Parameters
    ///
    /// * `position` - the position to convert
    fn position_as_index(position: Nlbf64) -> usize {
        (position.value() * (MAX_POSITION as f64)) as usize
    }

    /// Checks if the audio vector is sufficient to handle the current window and
    /// expands the vector if necessary.
    fn check_audio_vector(&mut self) {
        let end = self.current_end_index();
        let difference: i128 = (end as i128) - ((self.audio.len() as i128) - 1);
        if difference > 0 {
            self.audio
                .extend((0..difference).map(|_| Nlbf64::from(0.0f64)));
        }
    }

    /// Returns the index of the start of the current window.
    fn current_start_index(&self) -> usize {
        Self::position_as_index(self.position)
    }

    /// Returns the index of the end of the current window.
    fn current_end_index(&self) -> usize {
        self.current_start_index() + OUTPUT_WINDOW - 1
    }

    /// Returns the currently selected read.
    pub fn current_read(&mut self) -> Vec<Nlbf64> {
        self.check_audio_vector();
        let start_index = self.current_start_index();
        let end_index = self.current_end_index();
        (&self.audio[start_index..=end_index])
            .iter()
            .map(|value| *value)
            .collect()
    }

    /// Writes the current output window to the underlying audio vector.
    ///
    /// # Parameters
    ///
    /// * `information` - the information vector to write to the current output window
    fn set_current_output_window(&mut self, information: Vec<Option<SimpleNeuron>>) {
        self.check_audio_vector();
        let start = self.current_start_index();
        for (index, info_option) in information.into_iter().enumerate() {
            if let Some(info) = info_option {
                self.audio[start + index] = info.current_potential();
            }
        }
    }
}

impl Output<Vec<i16>, SimpleNeuron> for SimpleNeuronAudioSixteenOutputSensor {
    fn get_output(&mut self, information: Vec<Option<SimpleNeuron>>) -> Vec<i16> {
        self.set_current_output_window(information);
        self.audio
            .iter()
            .map(|value| nlbf64_to_audio_sixteen(*value))
            .collect()
    }

    fn handle_feedback_substrate_changes(
        &mut self,
        changes: std::collections::HashMap<usize, SimpleNeuron>,
        current_output_information: Vec<Option<SimpleNeuron>>,
    ) -> Option<Vec<SimpleNeuron>> {
        let halflifes: Vec<f64> = current_output_information
            .iter()
            .map(|opt| opt.as_ref())
            .map(|neuron_option| {
                neuron_option
                    .map_or(POTENTIAL_HALFLIFE_TIME_MAXIMUM, SimpleNeuron::potential_halflife_time)
            })
            .collect();
        self.set_current_output_window(current_output_information);
        changes.get(&0).map(|position_neuron| {
            self.position = position_neuron.current_potential();
            self.current_read()
                .into_iter()
                .zip(halflifes.into_iter())
                .map(|(potential, halflife)| SimpleNeuron::new(potential, halflife))
                .collect()
        })
    }

    fn is_finished(&self, information: SimpleNeuron) -> bool {
        information.current_potential().value() > 0.9
    }

    fn random() -> Self {
        Self::new()
    }
}

impl CrossOver for SimpleNeuronAudioSixteenOutputSensor {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, _other: &Self) -> Self {
        Self::new()
    }
}
