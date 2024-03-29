use std::{collections::VecDeque, num::NonZeroUsize, sync::Arc};

use parking_lot::Mutex;
use rand::{thread_rng, Rng};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{abnn::dendrite::ACTIVATION_POTENTIAL_REINFORCEMENT, evolution::helper::Iteration};

use super::{dendrite::Dendrite, interactor::Interactor, neuron::Neuron};

/// The maximum number of prediction errors that is remembered for
/// performance evaluation.
pub const NUMBER_OF_REMEMBERED_PREDICTION_ERRORS: usize = 10000;

pub struct AssociationBasedNeuralNetwork {
    interactors: Vec<Arc<Mutex<dyn Interactor + Send + Sync>>>,
    neurons: Vec<Arc<Neuron>>,
    current_time: Iteration,
    last_evaluation_time: Iteration,
    prediction_errors: VecDeque<f64>,
}

impl AssociationBasedNeuralNetwork {
    pub fn new(
        number_of_neurons: NonZeroUsize,
        dendrites_per_neuron: NonZeroUsize,
        interactors: Vec<Arc<Mutex<dyn Interactor + Send + Sync>>>,
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
        for interactor in interactors.iter() {
            let n_input = interactor.lock().number_of_input_neurons();
            let n_output = interactor.lock().number_of_output_neurons();
            if n_input > 0 {
                interactor.lock().set_input_neurons(
                    (0..n_input)
                        .map(|_| remove_random(&mut available_neurons))
                        .collect(),
                );
            }
            if n_output > 0 {
                interactor.lock().set_output_neurons(
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
                .filter(|target| !Arc::ptr_eq(neuron, target))
                .map(|target| Arc::clone(target))
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
            last_evaluation_time: Iteration::new(),
            prediction_errors: VecDeque::with_capacity(NUMBER_OF_REMEMBERED_PREDICTION_ERRORS),
        }
    }

    pub fn run_until_evaluation(&mut self) {
        self.interactors
            .par_iter()
            .for_each(|interactor| interactor.lock().initialise_new_evaluation());
        let mut dendrites_activated_since_last_evaluation: u64 = 0;
        let mut request_evaluation = false;
        let mut last_activated_dendrites: Vec<Arc<Dendrite>> = Vec::new();
        while !request_evaluation && self.current_time - self.last_evaluation_time < 10000 {
            self.interactors
                .par_iter()
                .for_each(|interactor| interactor.lock().update_iteration(self.current_time));
            let activated_dendrites: Vec<Arc<Dendrite>> = self
                .neurons
                .par_iter()
                .map(|neuron| (neuron, neuron.value_at_timepoint(self.current_time), neuron.try_trigger_action_potential(self.current_time)))
                .flat_map(|(neuron, source_value, dendrites)| {
                    // Cross-reference learning and reinforcment of dendrite weights
                    // when an action potential is triggered.
                    if dendrites.len() > 0 {
                        for last_activated_dendrite in last_activated_dendrites.iter() {
                            if Arc::ptr_eq(&last_activated_dendrite.target(), neuron) {
                                last_activated_dendrite.set_weight(
                                    last_activated_dendrite.weight()
                                        + ACTIVATION_POTENTIAL_REINFORCEMENT,
                                );
                            }
                        }
                    }
                    for dendrite in dendrites.iter() {
                        dendrite.trigger(source_value, self.current_time);
                    }
                    dendrites
                })
                .collect();
            dendrites_activated_since_last_evaluation += activated_dendrites.len() as u64;
            last_activated_dendrites = activated_dendrites;
            self.current_time = self.current_time.increment();
            request_evaluation = self
                .interactors
                .iter()
                .any(|interactor| interactor.lock().request_evaluation());
            }
            
            let interactor_prediction_errors: Vec<f64> = self
            .interactors
            .iter()
            .filter_map(|interactor| interactor.lock().evalute_results())
            .collect();
            let prediction_error = interactor_prediction_errors.iter().sum::<f64>()
            / (interactor_prediction_errors.len() as f64);
        self.update_prediction_errors(prediction_error);
        
        let normalised_prediction_error =
        if let Some(mean_prediction_error) = self.mean_prediction_error() {
                prediction_error - mean_prediction_error
            } else {
                0.0
            };
        if dendrites_activated_since_last_evaluation > 0 {
            // Resets the amount of times the dendrite was acitvated.
            self.neurons.par_iter().for_each(|neuron| {
                neuron.dendrites().into_iter().for_each(|dendrite| {
                    let scaling_factor = dendrite.times_activated() as f64
                        / dendrites_activated_since_last_evaluation as f64;
                        let original_weight = dendrite.weight();
                        let updated_weight = if original_weight >= 0.0 {
                            - scaling_factor * normalised_prediction_error
                        } else {
                            scaling_factor * normalised_prediction_error
                        };
                        dendrite.add_weight_after_evaluation(updated_weight);
                        dendrite.set_times_activated(0)
                    });
                });
            }
            println!("{}", self.current_time - self.last_evaluation_time);
            self.last_evaluation_time = self.current_time;
    }

    fn update_prediction_errors(&mut self, prediction_error: f64) {
        self.prediction_errors.push_back(prediction_error);
        if self.prediction_errors.len() > NUMBER_OF_REMEMBERED_PREDICTION_ERRORS {
            self.prediction_errors.pop_front();
        }
    }

    pub fn mean_prediction_error(&self) -> Option<f64> {
        if self.prediction_errors.is_empty() {
            None
        } else {
            Some(self.prediction_errors.iter().sum::<f64>() / (self.prediction_errors.len() as f64))
        }
    }
}

fn remove_random<T>(values: &mut Vec<T>) -> T {
    if values.is_empty() {
        panic!("There is not enough elements for random removal.");
    }
    let random_index = thread_rng().gen_range(0..values.len());
    values.remove(random_index)
}
