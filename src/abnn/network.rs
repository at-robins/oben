use std::{collections::VecDeque, num::NonZeroUsize, sync::Arc};

use rand::{thread_rng, Rng};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::evolution::{helper::Iteration, neuron};

use super::{
    dendrite::Dendrite,
    interactor::{self, Interactor},
    neuron::Neuron,
};

/// The maximum number of prediction errors that is remembered for
/// performance evaluation.
pub const NUMBER_OF_REMEMBERED_PREDICTION_ERRORS: usize = 100;

pub struct AssociationBasedNeuralNetwork {
    interactors: Vec<Arc<dyn Interactor + Send + Sync>>,
    neurons: Vec<Arc<Neuron>>,
    current_time: Iteration,
    prediction_errors: VecDeque<f64>,
}

impl AssociationBasedNeuralNetwork {
    pub fn new(
        number_of_neurons: NonZeroUsize,
        dendrites_per_neuron: NonZeroUsize,
        interactors: Vec<Arc<dyn Interactor + Send + Sync>>,
    ) -> Self {
        if dendrites_per_neuron >= number_of_neurons {
            panic!("The number of dendrites per neuron exceeds the total number of targetable neurons.");
        }
        // Creates empty neurons.
        let neurons: Vec<Arc<Neuron>> = (0..number_of_neurons.get())
            .map(|_| Arc::new(Neuron::new()))
            .collect();
        // Sets the input and output neurons at random.
        let mut available_neurons = neurons.clone();
        for interactor in interactors {
            let n_input = interactor.number_of_input_neurons();
            let n_output = interactor.number_of_output_neurons();
            if n_input > 0 {
                interactor.set_input_neurons(
                    (0..n_input)
                        .map(|_| remove_random(&mut available_neurons))
                        .collect(),
                );
            }
            if n_output > 0 {
                interactor.set_output_neurons(
                    (0..n_output)
                        .map(|_| remove_random(&mut available_neurons))
                        .collect(),
                );
            }
        }
        // Add the dendrites with random targets and weights.
        for neuron in neurons.iter() {
            let mut available_targets: Vec<Arc<Neuron>> = neurons
                .iter()
                .filter(|target| Arc::ptr_eq(neuron, target))
                .map(|target| *target)
                .collect();
            for _ in 0..dendrites_per_neuron.get() {
                neuron.add_dendrite(Arc::new(Dendrite::new(
                    remove_random(&mut available_targets),
                    thread_rng().gen_range(-1.0..=1.0),
                )));
            }
        }
        // Construct the network.
        Self {
            interactors,
            neurons,
            current_time: Iteration::new(),
            prediction_errors: VecDeque::with_capacity(NUMBER_OF_REMEMBERED_PREDICTION_ERRORS),
        }
    }

    pub fn run_until_evaluation(&mut self) {
        let mut dendrites_activated_since_last_evaluation: u64 = 0;
        let mut request_evaluation = false;
        while !request_evaluation {
            self.interactors
                .par_iter()
                .for_each(|interactor| interactor.update_iteration(self.current_time));
            let activated_dendrites: Vec<Arc<Dendrite>> = self
                .neurons
                .par_iter()
                .map(|neuron| neuron.try_trigger_action_potential(self.current_time))
                .flat_map(|dendrites| {
                    for dendrite in dendrites {
                        dendrite.trigger(self.current_time);
                    }
                    dendrites
                })
                .collect();
            dendrites_activated_since_last_evaluation += activated_dendrites.len() as u64;
            self.current_time.increment();
            request_evaluation = self
                .interactors
                .iter()
                .any(|interactor| interactor.request_evaluation());
        }

        todo!(); // Generate an prediction error.
        let prediction_error = 1.0;
        self.update_prediction_errors(prediction_error);

        let normalized_prediction_error = prediction_error - self.mean_prediction_error();

        if dendrites_activated_since_last_evaluation > 0 {
            // Resets the amount of times the dendrite was acitvated.
            self.neurons.par_iter().for_each(|neuron| {
                neuron.dendrites().into_iter().for_each(|dendrite| {
                    let scaling_factor = dendrite.times_activated() as f64
                        / dendrites_activated_since_last_evaluation as f64;
                    let updated_weight =
                        dendrite.weight() - scaling_factor * normalized_prediction_error;
                    dendrite.set_weight(updated_weight);
                    dendrite.set_times_activated(0)
                });
            });
        }
    }

    fn update_prediction_errors(&mut self, prediction_error: f64) {
        self.prediction_errors.push_back(prediction_error);
        if self.prediction_errors.len() > NUMBER_OF_REMEMBERED_PREDICTION_ERRORS {
            self.prediction_errors.pop_front();
        }
    }

    pub fn mean_prediction_error(&self) -> f64 {
        self.prediction_errors.into_iter().sum::<f64>() / (self.prediction_errors.len() as f64)
    }
}

fn remove_random<T>(values: &mut Vec<T>) -> T {
    if values.is_empty() {
        panic!("There is not enough elements for random removal.");
    }
    let random_index = thread_rng().gen_range(0..values.len());
    values.remove(random_index)
}
