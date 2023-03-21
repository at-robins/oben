use std::{collections::VecDeque, num::NonZeroUsize, sync::Arc};

use parking_lot::Mutex;
use rand::{thread_rng, Rng};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::evolution::helper::Iteration;

use super::{
    dendrite::Dendrite, interactor::Interactor, neuron::Neuron, parameters::ConfigurableParameters,
};

/// The maximum number of prediction errors that is remembered for
/// performance evaluation.
pub const NUMBER_OF_REMEMBERED_PREDICTION_ERRORS: usize = 200;

pub struct AssociationBasedNeuralNetwork {
    interactors: Vec<Arc<Mutex<dyn Interactor + Send + Sync>>>,
    neurons: Vec<Arc<Neuron>>,
    current_time: Iteration,
    last_evaluation_time: Iteration,
    prediction_errors: VecDeque<Vec<f64>>,
    configuration: Arc<ConfigurableParameters>,
    total_number_of_dendrites: usize,
}

impl AssociationBasedNeuralNetwork {
    pub fn new(
        number_of_neurons: NonZeroUsize,
        dendrites_per_neuron: NonZeroUsize,
        inhibitory_dendrites_per_neuron: usize,
        interactors: Vec<Arc<Mutex<dyn Interactor + Send + Sync>>>,
        configuration: ConfigurableParameters,
    ) -> Self {
        if dendrites_per_neuron.get() + inhibitory_dendrites_per_neuron >= number_of_neurons.get() {
            panic!("The number of dendrites per neuron exceeds the total number of targetable neurons.");
        }
        let shared_configuration = Arc::new(configuration);
        // Creates empty neurons.
        let neurons: Vec<Arc<Neuron>> = (0..number_of_neurons.get())
            .map(|_| Arc::new(Neuron::new(Arc::clone(&shared_configuration))))
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
                let dendrite = Arc::new(Dendrite::new(
                    remove_random(&mut available_targets),
                    Arc::clone(neuron),
                    thread_rng().gen_range(0.0..=1.0),
                    false,
                    Arc::clone(&shared_configuration),
                ));
                neuron.add_outgoing_dendrite(Arc::clone(&dendrite));
                dendrite.target().add_ingoing_dendrite(dendrite);
            }
            if inhibitory_dendrites_per_neuron > 0 {
                for _ in 0..inhibitory_dendrites_per_neuron {
                    let dendrite = Arc::new(Dendrite::new(
                        remove_random(&mut available_targets),
                        Arc::clone(neuron),
                        thread_rng().gen_range(0.0..=1.0),
                        true,
                        Arc::clone(&shared_configuration),
                    ));
                    neuron.add_outgoing_dendrite(Arc::clone(&dendrite));
                    dendrite.target().add_ingoing_dendrite(dendrite);
                }
            }
        }
        // Construct the network.
        Self {
            interactors,
            neurons,
            current_time: Iteration::new(),
            last_evaluation_time: Iteration::new(),
            prediction_errors: VecDeque::with_capacity(NUMBER_OF_REMEMBERED_PREDICTION_ERRORS),
            configuration: shared_configuration,
            total_number_of_dendrites: number_of_neurons.get() * dendrites_per_neuron.get(),
        }
    }

    pub fn run_until_evaluation(&mut self) {
        self.interactors
            .par_iter()
            .for_each(|interactor| interactor.lock().initialise_new_evaluation());
        let mut dendrites_activated_since_last_evaluation: u64 = 0;
        let mut request_evaluation = false;
        let mut last_activated_dendrites: Vec<Arc<Dendrite>> = Vec::new();
        let mut network_activity: f64 = self
            .configuration
            .dendrite_global_activity_regulation_midpoint();
        while !request_evaluation && self.current_time - self.last_evaluation_time < 10000 {
            self.interactors
                .par_iter()
                .for_each(|interactor| interactor.lock().update_iteration(self.current_time));
            let activated_dendrites: Vec<Arc<Dendrite>> = self
                .neurons
                .par_iter()
                .map(|neuron| (neuron, neuron.try_trigger_action_potential(self.current_time)))
                .flat_map(|(neuron, activated_dendrites)| {
                    // Cross-reference learning and reinforcment of dendrite weights
                    // when an action potential is triggered.
                    if activated_dendrites.len() > 0 {
                        for last_activated_dendrite in last_activated_dendrites.iter() {
                            if Arc::ptr_eq(&last_activated_dendrite.target(), neuron) {
                                last_activated_dendrite.set_weight(
                                    last_activated_dendrite.weight()
                                        + self
                                            .configuration
                                            .dendrite_activation_potential_reinforcement(),
                                );
                            }
                        }
                    }
                    for dendrite in activated_dendrites.iter() {
                        dendrite.trigger(network_activity, self.current_time);
                    }
                    activated_dendrites
                })
                .collect();
            network_activity = activated_dendrites.len() as f64
                / ((self.neurons.len() * self.total_number_of_dendrites) as f64);
            dendrites_activated_since_last_evaluation += activated_dendrites.len() as u64;
            last_activated_dendrites = activated_dendrites;
            self.current_time = self.current_time.increment();
            request_evaluation = self
                .interactors
                .iter()
                .any(|interactor| interactor.lock().request_evaluation());
        }

        // let interactor_prediction_errors: Vec<ErrorPropagationImpulse> = self
        //     .interactors
        //     .iter()
        //     .filter_map(|interactor| interactor.lock().evalute_results())
        //     .flat_map(|results| results)
        //     .collect();
        // let prediction_error = interactor_prediction_errors.iter().sum::<f64>()
        //     / (interactor_prediction_errors.len() as f64);
        // self.update_prediction_errors(prediction_error);

        // let normalised_prediction_error =
        //     if let Some(mean_prediction_error) = self.mean_prediction_error() {
        //         prediction_error - mean_prediction_error
        //     } else {
        //         0.0
        //     };
        if dendrites_activated_since_last_evaluation > 0 {
            // Resets the amount of times the dendrite was acitvated.
            // self.neurons.par_iter().for_each(|neuron| {
            //     neuron
            //         .outgoing_dendrites()
            //         .into_iter()
            //         .for_each(|dendrite| {
            //             let scaling_factor = dendrite.times_activated() as f64
            //                 / dendrites_activated_since_last_evaluation as f64;
            //             let updated_weight = if !dendrite.is_inhibitory() {
            //                 -scaling_factor * normalised_prediction_error
            //             } else {
            //                 scaling_factor * normalised_prediction_error
            //             };
            //             dendrite.add_weight_after_evaluation(updated_weight);
            //             dendrite.set_times_activated(0)
            //         });
            // });
            let prediction_error: Vec<f64> = self
                .interactors
                .par_iter()
                .filter_map(|interactor| interactor.lock().evalute_results())
                .flat_map(|results| results)
                .map(|(neuron, error)| {
                    neuron.propagate_error(error, dendrites_activated_since_last_evaluation);
                    println!("{:?}", error);
                    error.error()
                })
                .collect();
            self.update_prediction_errors(prediction_error);
        }
        self.neurons
            .par_iter()
            .map(|neuron| neuron.outgoing_dendrites())
            .for_each(|dendrites| {
                dendrites
                    .iter()
                    .for_each(|dendrite| dendrite.evaluate_error_propagation())
            });
        // println!("{}", self.current_time - self.last_evaluation_time);
        self.last_evaluation_time = self.current_time;
    }

    fn update_prediction_errors(&mut self, prediction_error: Vec<f64>) {
        self.prediction_errors.push_back(prediction_error);
        if self.prediction_errors.len() > NUMBER_OF_REMEMBERED_PREDICTION_ERRORS {
            self.prediction_errors.pop_front();
        }
    }

    pub fn mean_prediction_error(&self) -> Option<Vec<f64>> {
        if self.prediction_errors.is_empty() {
            None
        } else {
            Some(
                self.prediction_errors
                    .iter()
                    .fold(Vec::new(), |acc, x| {
                        if acc.is_empty() {
                            x.clone()
                        } else {
                            acc.iter().zip(x.iter()).map(|(a, b)| a.abs() + b.abs()).collect()
                        }
                    })
                    .iter()
                    .map(|error| error / (self.prediction_errors.len() as f64))
                    .collect(),
            )
        }
    }

    fn reset_neurons(&mut self) {
        for neuron in &self.neurons {
            neuron.set_value(0.2);
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

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ErrorPropagationImpulse {
    /// The ID of the output neuron from which the error originates.
    output_id: usize,
    /// The sum of consecutive dendrite weights.
    impact_factor: f64,
    /// The initial prediciton error.
    error: f64,
}

impl ErrorPropagationImpulse {
    pub fn new(output_id: usize, error: f64) -> Self {
        Self {
            output_id,
            impact_factor: 1.0,
            error,
        }
    }

    pub fn with_updated_impact_factor(&self, dendrite_weight: f64) -> Self {
        Self {
            output_id: self.output_id,
            impact_factor: dendrite_weight, //self.impact_factor *
            error: self.error,
        }
    }

    pub fn output_id(&self) -> usize {
        self.output_id
    }

    pub fn impact_factor(&self) -> f64 {
        self.impact_factor
    }

    pub fn error(&self) -> f64 {
        self.error
    }
}

unsafe impl Send for ErrorPropagationImpulse {}

unsafe impl Sync for ErrorPropagationImpulse {}
