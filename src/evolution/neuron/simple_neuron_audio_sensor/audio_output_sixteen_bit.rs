use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize, de::DeserializeOwned};

use crate::evolution::{chemistry::{Output, Input}, gene::{CrossOver, Genome}, helper::Nlbf64, neuron::{SimpleNeuron, SimpleDendriteActivationPotential, SimpleDendriteThreshold}};

const NUMBER_OF_PHONEMES: usize = 60;
const MAXIMUM_PHONEME_LENGTH: usize = 2000;
/// The maximum position that can be accessed by an [`SimpleNeuronAudioSixteenOutputSensor`].
pub const MAX_POSITION: usize = 10_000_000;

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
    phonemes: Vec<Vec<Nlbf64>>,
    audio: Vec<Nlbf64>,
    position: Nlbf64,
}

impl SimpleNeuronAudioSixteenOutputSensor {
    /// Creates a new `SimpleNeuronAudioSixteenOutputSensor`.
    ///
    /// # Parameters
    ///
    /// * `phonemes` - the audio building blocks
    pub fn new(phonemes: Vec<Vec<Nlbf64>>) -> Self {
        Self {
            phonemes,
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

    /// Returns the specified phoneme as index of the available phonemes.
    ///
    /// # Parameters
    ///
    /// * `phoneme_information` - the phoneme information to convert
    fn phoneme_as_index(&self, phoneme_information: Nlbf64) -> usize {
        let calculated_index =
            (phoneme_information.value() * (self.phonemes.len() as f64)).floor() as usize;
        // This fixes the fact that a value of exactly 1.0 corresponds to `self.phonemes.len()`.
        if calculated_index >= self.phonemes.len() {
            self.phonemes.len() - 1
        } else {
            calculated_index
        }
    }

    /// Returns the specifed phoneme.
    ///
    /// # Parameters
    ///
    /// * `phoneme_information` - the phoneme information to convert
    fn phoneme(&self, phoneme_information: Nlbf64) -> &Vec<Nlbf64> {
        self.phonemes
            .get(self.phoneme_as_index(phoneme_information))
            .expect("The phoneme bounds were checked.")
    }

    /// Checks if the audio vector is sufficient to handle the current window and
    /// expands the vector if necessary.
    ///
    /// # Parameters
    ///
    /// * `current_phoneme` - the currently selected phoneme
    fn check_audio_vector(&mut self, current_phoneme: Nlbf64) {
        let end = self.current_end_index(current_phoneme);
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
    ///
    /// # Parameters
    ///
    /// * `current_phoneme` - the currently selected phoneme
    fn current_end_index(&self, current_phoneme: Nlbf64) -> usize {
        self.current_start_index() + self.phoneme(current_phoneme).len() - 1
    }

    /// Writes the current output window to the underlying audio vector.
    ///
    /// # Parameters
    ///
    /// * `information` - the information vector to write to the current output window
    fn write_current_phoneme(&mut self, information: Vec<Option<SimpleNeuron>>) {
        let phoneme_index = information
            .get(0)
            .expect("There must be a phoneme output substrate.")
            .as_ref();
        // Mode of the phoneme interaction with the audio signal, e.g. addition or substraction.
        let mode = information
            .get(0)
            .expect("There must be a mode output substrate.")
            .as_ref();
        if phoneme_index.is_some() && mode.is_some() {
            let phoneme_index = phoneme_index.unwrap().current_potential();
            let mode = mode.unwrap().current_potential();
            self.check_audio_vector(phoneme_index);
            let start = self.current_start_index();
            let end = self.current_end_index(phoneme_index);
            for audio_index in start..=end {
                // The bounds were checked before, so the unwrap is fine.
                let new_value = if mode.value() < 0.5 {
                    (self.audio[audio_index].value()
                        + self
                            .phoneme(phoneme_index)
                            .get(audio_index - start)
                            .unwrap()
                            .value())
                        / 2.0
                } else {
                    (self.audio[audio_index].value() + 1.0
                        - self
                            .phoneme(phoneme_index)
                            .get(audio_index - start)
                            .unwrap()
                            .value())
                        / 2.0
                };
                self.audio[audio_index] = new_value.into();
            }
        }
    }
}

impl Output<Vec<i16>, SimpleNeuron> for SimpleNeuronAudioSixteenOutputSensor {
    fn get_output(&mut self, information: Vec<Option<SimpleNeuron>>) -> Vec<i16> {
        self.write_current_phoneme(information);
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
        self.write_current_phoneme(current_output_information);
        changes.get(&0).and_then(|position_neuron| {
            self.position = position_neuron.current_potential();
            None
        })
    }

    fn is_finished(&self, information: SimpleNeuron) -> bool {
        information.current_potential().value() > 0.9
    }

    fn random() -> Self {
        let phonemes = (0..NUMBER_OF_PHONEMES)
            .map(|_| {
                let phoneme_length = thread_rng().gen_range(1..MAXIMUM_PHONEME_LENGTH);
                (0..phoneme_length).map(|_| thread_rng().gen()).collect()
            })
            .collect();
        Self::new(phonemes)
    }
}

impl CrossOver for SimpleNeuronAudioSixteenOutputSensor {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        let crossover_phonemes = self
            .phonemes
            .iter()
            .zip(other.phonemes.iter())
            .map(|(a, b)| a.cross_over(b))
            .collect();

        Self::new(crossover_phonemes)
    }
}


pub fn mutate_phoneme<
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, SimpleNeuron>,
>(
    genome: &Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,
        Vec<i16>,
        SimpleNeuronAudioSixteenOutputSensor,
    >,
) -> Option<
    Genome<
        SimpleDendriteActivationPotential,
        SimpleDendriteThreshold,
        SimpleNeuron,
        InputElementType,
        InputSensorType,Vec<i16>,
        SimpleNeuronAudioSixteenOutputSensor,
    >,
> {
    let mut mutated_genome = genome.duplicate();
    let phoneme = mutated_genome
        .output_mut()
        .output_mut()
        .phonemes
        .get_mut(thread_rng().gen_range(0..NUMBER_OF_PHONEMES))
        .unwrap();
    let mutation_index: usize = thread_rng().gen_range(0..phoneme.len()); 
    let mutated_phoneme_value = Nlbf64::flip_random_bit(*(phoneme.get(mutation_index).unwrap()));
    phoneme[mutation_index] = mutated_phoneme_value;
    Some(mutated_genome)
}
