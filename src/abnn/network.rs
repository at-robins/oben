use std::{
    collections::{HashSet, VecDeque},
    num::NonZeroUsize,
    sync::Arc,
};

use by_address::ByAddress;
use parking_lot::Mutex;
use rand::{thread_rng, Rng};
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IntoParallelIterator,
};

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
        let input_neurons: Vec<Arc<Neuron>> = interactors
            .iter()
            .flat_map(|interactor| interactor.lock().input_neurons())
            .collect();
        let output_neurons: Vec<Arc<Neuron>> = interactors
            .iter()
            .flat_map(|interactor| interactor.lock().output_neurons())
            .collect();
        // Add the dendrites with random targets and weights to every neuron except output neurons.
        for neuron in Self::neurons_exclude(&neurons, &output_neurons).iter() {
            let mut available_targets: Vec<Arc<Neuron>> =
                Self::neurons_exclude(&neurons, &input_neurons)
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

    fn neurons_exclude(list: &Vec<Arc<Neuron>>, exclude: &Vec<Arc<Neuron>>) -> Vec<Arc<Neuron>> {
        list.iter()
            .filter(|neuron| {
                let mut include = true;
                for exclude_neuron in exclude.iter() {
                    if Arc::ptr_eq(neuron, &exclude_neuron) {
                        include = false;
                        break;
                    }
                }
                include
            })
            .map(|neuron| Arc::clone(neuron))
            .collect()
    }

    pub fn run_until_evaluation(&mut self) {
        let mut update_neurons: HashSet<ByAddress<Arc<Neuron>>> = self
            .interactors
            .par_iter()
            .flat_map(|interactor| {
                let mut interactor_lock = interactor.lock();
                interactor_lock.initialise_new_evaluation(self.current_time);
                interactor_lock.input_neurons()
            })
            .flat_map(|neuron| neuron.targets())
            .map(|neuron| ByAddress(neuron))
            .collect();
        let mut request_evaluation = false;
        while !request_evaluation {
            let updated_input_neurons: Vec<ByAddress<Arc<Neuron>>> = self
                .interactors
                .par_iter()
                .flat_map(|interactor| {
                    interactor
                        .lock()
                        .update_iteration(self.current_time, update_neurons.is_empty())
                })
                .flat_map(|neuron| neuron.targets())
                .map(|neuron| ByAddress(neuron))
                .collect();
            update_neurons.extend(updated_input_neurons);
            update_neurons = update_neurons
                .par_iter()
                .flat_map(|neuron| {
                    if neuron.update_value() {
                        neuron.targets()
                    } else {
                        Vec::new()
                    }
                })
                .map(|neuron| ByAddress(neuron))
                .collect();
            self.current_time = self.current_time.increment();
            request_evaluation = self
                .interactors
                .iter()
                .any(|interactor| interactor.lock().request_evaluation());
        }
        // for neuron in self.neurons.iter() {
        //     let dendrite_values_out: Vec<f64> = neuron
        //         .outgoing_dendrites()
        //         .iter()
        //         .map(|a| a.weight())
        //         .collect();
        //     let dendrite_values_in: Vec<f64> = neuron
        //         .ingoing_dendrites()
        //         .iter()
        //         .map(|a| a.weight())
        //         .collect();
        //     println!(
        //         "Neuron: {} ; Dendrites out: {:?} ; Dendrites in: {:?}",
        //         neuron.value(),
        //         dendrite_values_out,
        //         dendrite_values_in
        //     );
        // }

        // self.neurons
        //     .par_iter()
        //     .flat_map(|neuron| neuron.outgoing_dendrites())
        //     .for_each(|dendrite| {
        //         let diff = (dendrite.source().value() - dendrite.target().value().powi(2));
        //         // dendrite.set_weight(dendrite.weight() * (2.0 - diff));
        //     });

        let prediction_error: Vec<f64> = self
            .interactors
            .par_iter()
            .filter_map(|interactor| {
                let mut interactor_lock = interactor.lock();
                interactor_lock
                    .evalute_results()
                    .map(|result| (interactor_lock.output_neurons(), result))
            })
            .map(|(neurons, error)| {
                println!("{:?}", error);
                // TODO: reasonable code
                for neuron in neurons.into_iter() {
                    neuron.propagate_error(error);
                }
                error.error()
            })
            .collect();
        self.update_prediction_errors(prediction_error);

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
                            acc.iter()
                                .zip(x.iter())
                                .map(|(a, b)| a.abs() + b.abs())
                                .collect()
                        }
                    })
                    .iter()
                    .map(|error| error / (self.prediction_errors.len() as f64))
                    .collect(),
            )
        }
    }

    // fn reset_neurons(&mut self) {
    //     for neuron in &self.neurons {
    //         neuron.set_value(0.2);
    //     }
    // }
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
    /// The distance from the original error source / output neuron.
    distance: u64,
    /// The initial prediciton error.
    error: f64,
}

impl ErrorPropagationImpulse {
    pub fn new(output_id: usize, error: f64) -> Self {
        Self {
            output_id,
            distance: 0,
            error,
        }
    }

    pub fn with_updated_distance(&self) -> Self {
        Self {
            output_id: self.output_id,
            distance: self.distance + 1,
            error: self.error,
        }
    }

    pub fn output_id(&self) -> usize {
        self.output_id
    }

    pub fn distance(&self) -> u64 {
        self.distance
    }

    pub fn error(&self) -> f64 {
        self.error
    }
}

unsafe impl Send for ErrorPropagationImpulse {}

unsafe impl Sync for ErrorPropagationImpulse {}
