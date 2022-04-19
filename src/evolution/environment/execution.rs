//! The `execution` module contains the executive setup of the evolutionary network.

use crate::evolution::chemistry::{Input, Output};
use crate::evolution::gene::CrossOver;
use crate::evolution::helper::ScalingFactor;
use rand::{thread_rng, Rng};
use rayon::prelude::*;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use super::super::chemistry::{Information, Reaction, State};
use super::super::gene::Genome;
use super::super::population::{Individual, Organism, OrganismInformation, Population};
use super::super::resource::Resource;
use super::configuration::Environment;
use super::MutationCompendium;
use uuid::Uuid;

/// An `EcologicalNiche` containing a [`Population`] and applying selective pressure.
///
/// [`Population`]: ../population/struct.Population.html
pub struct EcologicalNiche<
    SupplierResultInformationType,
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
    OutputElementType,
    OutputSensorType,
> {
    inner: Arc<
        InnerEcologicalNiche<
            SupplierResultInformationType,
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    >,
}

impl<
        SupplierResultInformationType: 'static,
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    EcologicalNiche<
        SupplierResultInformationType,
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Creates a new `EcologicalNiche` executing the evolutionary network by repeated
    /// mutagenesis, fitness evaluation and growth of a [`Population`].
    ///
    /// # Parameters
    ///
    /// * `environment` - the [`Environment`] that defines the basic properties of the network
    /// * `population` - the starting [`Population`] to alter during execution of the network
    /// * `supplier_function` - the function supplying the [`Population`] with test examples
    /// * `fitness_function` - the function evaluating the [`Population`]'s fitness based on
    /// the result obtained after supplying an example
    /// * `mutations` - a list of mutations that might occur during genome duplication
    ///
    /// [`Environment`]: ./struct.Environment.html
    /// [`Population`]: ../population/struct.Population.html
    pub fn new(
        environment: Environment,
        population: Population<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
        supplier_function: Box<
            dyn Fn() -> (InputElementType, SupplierResultInformationType) + Send + Sync + 'static,
        >,
        fitness_function: Box<
            dyn Fn(
                    Vec<OrganismInformation<SupplierResultInformationType, OutputElementType>>,
                    ScalingFactor,
                ) -> f64
                + Send
                + Sync
                + 'static,
        >,
        mutations: MutationCompendium<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    ) -> Self {
        EcologicalNiche {
            inner: Arc::new(InnerEcologicalNiche {
                environment: Arc::new(environment),
                population: Arc::new(Mutex::new(population)),
                supplier_function,
                fitness_function,
                mutations,
            }),
        }
    }

    /// Initialises the network.
    fn initialise(&self) {
        self.environment().initialise();
    }

    /// Returns the [`Environment`] of the network.
    ///
    /// [`Environment`]: ./struct.Environment.html
    fn environment(&self) -> Arc<Environment> {
        self.inner.environment.clone()
    }

    pub fn breath_life(&self) {
        // Initialise the environment.
        println!("Initialiasing the network...");
        self.initialise();
        // Start the network.
        println!("Starting execution...");
        let mut generation: u64 = 0;
        let mut fitness_scaling: ScalingFactor =
            self.environment().initial_fitness_scaling_factor();
        let mut start = Instant::now();
        loop {
            let spawn_counter = Arc::new(Mutex::new(0u32));
            let mating_counter = Arc::new(Mutex::new(0u32));
            generation += 1;
            println!("Generation {}", generation);
            // Age the population by a generation.
            self.inner.increment_age();
            // Challenge the organisms in the population.
            self.inner.individuals().par_iter().for_each(|individual| {
                Self::spawn_organism(self.inner.clone(), individual.clone(), fitness_scaling);
                *spawn_counter.lock().unwrap() += 1;
                if *spawn_counter.lock().unwrap() % 1000 == 0 {
                    println!("     Spawn Organism {}", *spawn_counter.lock().unwrap());
                }
            });
            // Distribute resources neccesarry for mating based on fitness.
            self.inner.distribute_resources();
            // Mate the organisms of the population and add offspring to the population.
            self.inner.individuals().par_iter().for_each(|individual| {
                Self::mate_organism(self.inner.clone(), individual.clone());
                *mating_counter.lock().unwrap() += 1;
                if *mating_counter.lock().unwrap() % 1000 == 0 {
                    println!("     Mate Organism {}", *mating_counter.lock().unwrap());
                }
            });
            // Kill individuals on statistical basis.
            self.inner
                .individuals()
                .par_iter()
                .filter(|individual| self.inner.died((*individual).clone()))
                .for_each(|individual| {
                    self.inner.remove_individual(individual.clone());
                });
            // Recycle resources.
            self.inner.recycle();
            // Print statistics
            let mut res: f64 = self
                .inner
                .individuals()
                .par_iter()
                .map(|a| {
                    InnerEcologicalNiche::<
                        SupplierResultInformationType,
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                        OutputElementType,
                        OutputSensorType,
                    >::get_accumulated_resources(a.clone())
                        + 1.0
                })
                .sum();
            res += self.inner.resources().total();
            let mean_fitness: f64 = self.inner.population_mean_fitness();
            println!("Size: {} : Bytes: {} ; Fitness: {} ; Fitness Scaling: {} ; Total Resources: {} ; Resources: {:?}",
                self.inner.population_size(),
                self.inner.population_mean_genome_size(),
                mean_fitness,
                fitness_scaling.exponent(),
                res,
                self.inner.resources());
            // Save the population in regular intervalls with a timestamp and print some information.
            if start.elapsed() >= self.environment().population_save_intervall() {
                self.save_population();
                start = Instant::now();
            }
            // Modify the fitness function scaling factor.
            if mean_fitness > 0.5 {
                fitness_scaling.decrement();
            } else if mean_fitness < 0.4 {
                fitness_scaling.increment();
            }
        }
    }

    /// Saves the current population.
    fn save_population(&self) {
        println!("Waiting to save...");
        let population_id = self.environment().generate_uuid();
        let save_path = self.environment().population_path(&population_id);
        self.inner.save_population(save_path);
        println!("Saved!");
    }

    /// Creates and tests the [`Organism`] and updates its evaluated fitness.
    ///
    /// # Parameters
    ///
    /// * `inner` - the [`Environment`] the [`Organism`] is living in
    /// * `individual` - the [`Individual`] describing the [`Organism`] to test
    ///
    /// [`Environment`]: ./struct.Environment.html
    /// [`Organism`]: ../population/struct.Organism.html
    /// [`Individual`]: ../population/struct.Individual.html
    fn spawn_organism(
        inner: Arc<
            InnerEcologicalNiche<
                SupplierResultInformationType,
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
        fitness_scaling: ScalingFactor,
    ) {
        let tested = inner.testing(individual.clone());
        if tested {
            // Transcribe / translate the genome and test the organism.
            Self::add_fitness(
                individual.clone(),
                Self::test_organism(inner.clone(), individual.clone(), fitness_scaling),
            );
        }
    }

    /// Mates the [`Organism`] and adds its offspring to the [`Population`].
    ///
    /// # Parameters
    ///
    /// * `inner` - the [`Environment`] the [`Organism`] is living in
    /// * `individual` - the [`Individual`] describing the [`Organism`] to test
    ///
    /// [`Environment`]: ./struct.Environment.html
    /// [`Organism`]: ../population/struct.Organism.html
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Population`]: ../population/struct.Population.html
    fn mate_organism(
        inner: Arc<
            InnerEcologicalNiche<
                SupplierResultInformationType,
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) {
        let offspring = Self::get_offspring(individual.clone(), inner.clone());
        // Add the mutated offspring to the population.
        inner.append_population(offspring);
    }

    /// Tests the [`Organism`] and returns the evaluated fitness.
    ///
    /// # Parameters
    ///
    /// * `inner` - the [`Environment`] the [`Organism`] is living in
    /// * `individual` - the [`Individual`] the [`Organism`] to test comes from
    ///
    /// [`Environment`]: ./struct.Environment.html
    /// [`Organism`]: ../population/struct.Organism.html
    /// [`Individual`]: ../population/struct.Individual.html
    fn test_organism(
        inner: Arc<
            InnerEcologicalNiche<
                SupplierResultInformationType,
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
        fitness_scaling: ScalingFactor,
    ) -> f64 {
        let mut organism = inner.load_organism(individual.clone());
        let mut organism_informations = Vec::new();
        // Repeatedly test the organism and supply all the testing information to the fitness
        // function.
        for _ in 0..inner.environment.testing_repetitions() {
            let (input, result_information) = (inner.supplier_function)();
            organism.set_input(input);
            let run_time = organism.live(&inner.environment);
            let output = organism.get_result();
            organism_informations.push(OrganismInformation::new(
                output,
                result_information,
                inner.get_bytes(individual.clone()) * 8,
                run_time,
                *(&inner.environment.lifespan()),
                inner.get_associated_inputs(individual.clone()),
                inner.get_associated_outputs(individual.clone()),
                organism.binary_size(),
                *(&inner.environment.max_organism_size()),
            ));
        }
        (inner.fitness_function)(organism_informations, fitness_scaling)
    }

    /// Adds the specified fitness to the specified [`Individual`].
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    /// * `fitness` - the fitness to add
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn add_fitness(
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
        fitness: f64,
    ) {
        let mut ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        ind.evaluate_new_fitness(fitness)
    }

    /// Generates offspring by sexual reproduction of the [`Individual`].
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    /// * `inner` - the inner environment
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn get_offspring(
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
        inner: Arc<
            InnerEcologicalNiche<
                SupplierResultInformationType,
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
    ) -> Vec<
        Individual<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    > {
        let mut offspring = Vec::new();
        // Use the accumulated resources to produce offspring.
        let number_of_offspring = inner.spend_resources_for_mating(individual.clone());
        for _ in 0..number_of_offspring {
            let partner = inner.get_random_genome();
            let ind = individual
                .lock()
                .expect("A thread paniced while holding the individual's lock.");
            offspring.push(ind.mate_and_mutate(partner, &inner.mutations, &inner.environment));
        }
        offspring
    }
}

struct InnerEcologicalNiche<
    SupplierResultInformationType,
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
    OutputElementType,
    OutputSensorType,
> {
    environment: Arc<Environment>,
    population: Arc<
        Mutex<
            Population<
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
    >,
    supplier_function:
        Box<dyn Fn() -> (InputElementType, SupplierResultInformationType) + Send + Sync + 'static>,
    fitness_function: Box<
        dyn Fn(
                Vec<OrganismInformation<SupplierResultInformationType, OutputElementType>>,
                ScalingFactor,
            ) -> f64
            + Send
            + Sync
            + 'static,
    >,
    mutations: MutationCompendium<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
}

impl<
        SupplierResultInformationType,
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    InnerEcologicalNiche<
        SupplierResultInformationType,
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Checks if the specified [`Individual`] died of age.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`] to check for extinction
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn died(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> bool {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        self.environment.death_chance(ind.age()) >= thread_rng().gen_range(0.0..=1.0)
    }

    /// Checks if the specified [`Individual`] is juvenil and still needs testing.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`] to check
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn is_juvenil(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> bool {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        self.environment
            .max_testing_age()
            .map_or(true, |max_age| ind.age() < max_age)
    }

    /// Checks if the specified [`Individual`] should be testedby chance.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`] to check
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn testing_by_chance(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> bool {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        thread_rng().gen_range(0.0..=1.0) <= self.environment.testing_chance(ind.times_tested())
    }

    /// Checks if the specified [`Individual`] should be tested.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`] to check
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn testing(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> bool {
        self.is_juvenil(individual.clone()) || self.testing_by_chance(individual)
    }

    /// Return the UUID of the specified [`Individual`].
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn get_uuid(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> Uuid {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        *ind.uuid()
    }

    /// Returns the [`Resource`]s accumulated by the specified [`Individual`] that can be
    /// used to generate offspring.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Resource`]: ../resource/struct.Resource.html
    fn get_accumulated_resources(
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> f64 {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        ind.resources()
    }

    /// Spends the maximum amount of [`Resource`]s possible to generate offspring and returns the
    /// number of offspring generated this way.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Resource`]: ../resource/struct.Resource.html
    fn spend_resources_for_mating(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> usize {
        individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.")
            .spend_resources_for_mating()
    }

    /// Return the size in bytes of the specified [`Individual`].
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn get_bytes(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> usize {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        ind.bytes()
    }

    /// Return the number of associated inputs for the specified [`Individual`].
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn get_associated_inputs(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> usize {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        ind.associated_inputs()
    }

    /// Return the number of associated outputs for the specified [`Individual`].
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn get_associated_outputs(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> usize {
        let ind = individual
            .lock()
            .expect("A thread paniced while holding the individual's lock.");
        ind.associated_outputs()
    }

    /// Returns the copy of a [`Genome`] of a random [`Individual`]
    /// where the chances for picking an individual correspond to its
    /// respective fitness.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock or the population is
    /// completely extinct.
    ///
    /// [`Individual`]: ./struct.Individual.html
    fn get_random_genome(
        &self,
    ) -> Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType, OutputElementType, OutputSensorType> {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .random_genome_fitness_based()
            // Unwrapping is safe here since we cannot call this function on empty populations.
            .unwrap()
            .duplicate()
    }

    /// Write a snapshot of the current [`Population`] to a JSON file.
    ///
    /// # Parameters
    ///
    /// * `save_path` - the JSON file the [`Population`] should be written to
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock or the file could not be created.
    ///
    /// [`Population`]: ../population/struct.Population.html
    fn save_population<P: AsRef<Path>>(&self, save_path: P) {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .snapshot_to_file(&save_path)
            .expect(&format!("The file {:?} could not be created.", save_path.as_ref()));
    }

    /// Returns all [`Individual`]s.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    fn individuals(
        &self,
    ) -> Vec<
        Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    > {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .individuals()
    }

    /// Returns the number of [`Individual`]s in the [`Population`].
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Population`]: ../population/struct.Population.html
    fn population_size(&self) -> usize {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .size()
    }

    /// Ages all [`Individual`]s in the [`Population`] by 1 generation.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Population`]: ../population/struct.Population.html
    fn increment_age(&self) {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .increment_age()
    }

    /// Distributes available [`Resource`]s among the [`Population`]
    /// based on the fitness of the [`Individual`]s.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Population`]: ../population/struct.Population.html
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn distribute_resources(&self) {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .distribute_resources()
    }

    /// Returns the mean fitness of all [`Individual`]s in the [`Population`].
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Population`]: ../population/struct.Population.html
    fn population_mean_fitness(&self) -> f64 {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .mean_fitness()
    }

    /// Returns the mean [`Genome`] size in byte of all [`Individual`]s in the [`Population`].
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock.
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Population`]: ../population/struct.Population.html
    fn population_mean_genome_size(&self) -> f64 {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .mean_genome_size()
    }

    /// Removes the specified [`Individual`] and repatriates its accumulated [`Resource`]s.
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`] to remove
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock or if the [`Individual`]
    /// could not be removed.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Resource`]: ../resource/struct.Resource.html
    fn remove_individual(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) {
        let uuid = self.get_uuid(individual.clone());
        // An individual consumes 1.0 resources when being born, so this has to be repatriated
        // additionally to the accumulated resources.
        let resources = Self::get_accumulated_resources(individual) + 1.0;
        {
            self.population
                .lock()
                .expect("A thread paniced while holding the population lock.")
                .repatriate_resources(resources);
        }
        {
            self.population
                .lock()
                .expect("A thread paniced while holding the population lock.")
                .remove(uuid)
                .expect("The individual could not be removed.");
        }
    }

    /// Recycles inavailable [`Resource`]s at the end of a generation.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    fn recycle(&self) {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .recycle();
    }

    /// Returns the [`Resource`]s.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    fn resources(&self) -> Resource {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .resources()
    }

    /// Load the [`Organism`] corresponding to the specified [`Individual`]
    ///
    /// # Parameters
    ///
    /// * `individual` - the [`Individual`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual lock.
    ///
    /// [`Organism`]: ../population/struct.population.html
    /// [`Individual`]: ../population/struct.Individual.html
    fn load_organism(
        &self,
        individual: Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    ) -> Organism<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    > {
        individual
            .lock()
            .expect("Another thread panicked while holding the individual lock.")
            .genome()
            .translate()
    }

    /// Add the [`Individual`]s to the [`Population`].
    ///
    /// # Parameters
    ///
    /// * `individuals` - the [`Individual`]s to add
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual lock.
    ///
    /// [`Population`]: ../population/struct.Population.html
    /// [`Individual`]: ../population/struct.Individual.html
    fn append_population(
        &self,
        individuals: Vec<
            Individual<
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
    ) {
        self.population
            .lock()
            .expect("A thread paniced while holding the population lock.")
            .append(individuals);
    }
}
