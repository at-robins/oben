//! The `execution` module contains the executive setup of the evolutionary network.

extern crate bitvec;
extern crate rand;
extern crate rayon;

use rand::{thread_rng, Rng};
use rayon::prelude::*;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use super::configuration::Environment;
use super::super::chemistry::{Information, Reaction, State};
use super::super::gene::{Genome, GenomeMutation};
use super::super::population::{Individual, Population, Organism, OrganismInformation};
use super::super::resource::Resource;
use uuid::Uuid;

/// An `EcologicalNiche` containing a [`Population`] and applying selective pressure.
///
/// [`Population`]: ../population/struct.Population.html
pub struct EcologicalNiche<I, M, R, S, T> {
    inner: Arc<InnerEcologicalNiche<I, M, R, S, T>>,
    phantom: PhantomData<M>,
}

impl<I: 'static, M: GenomeMutation<R, S, T>, R: Reaction<T>, S: State<T>, T: Information> EcologicalNiche<I, M, R, S, T> {
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
    ///
    /// [`Environment`]: ./struct.Environment.html
    /// [`Population`]: ../population/struct.Population.html
    pub fn new(
        environment: Environment<M, R, S, T>,
        population: Population<R, S, T>,
        supplier_function: Box<dyn Fn() -> (Vec<T>, I)  + Send + Sync + 'static>,
        fitness_function: Box<dyn Fn(Vec<OrganismInformation<I, T>>) -> f64 + Send + Sync + 'static>) -> Self {
        EcologicalNiche {
            inner: Arc::new(InnerEcologicalNiche {
                environment: Arc::new(environment),
                population: Arc::new(Mutex::new(population)),
                supplier_function,
                fitness_function,
                phantom: PhantomData,
            }),
            phantom: PhantomData,
        }
    }

    /// Initialises the network.
    fn initialise(&self) {
        self.environment().initialise();
    }

    /// Returns the [`Environment`] of the network.
    ///
    /// [`Environment`]: ./struct.Environment.html
    fn environment(&self) -> Arc<Environment<M, R, S, T>> {
        self.inner.environment.clone()
    }

    pub fn breath_life(&self) {
        // Initialise the environment.
        println!("Initialiasing the network...");
        self.initialise();
        // Start the network.
        println!("Starting execution...");
        let mut generation: u64 = 0;
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
                Self::spawn_organism(self.inner.clone(), individual.clone());
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
            self.inner.individuals().par_iter()
                .filter(|individual| self.inner.died((*individual).clone()))
                .for_each(|individual| {
                    self.inner.remove_individual(individual.clone());
                });
            // Recycle resources.
            self.inner.recycle();
            // Print statistics
            let mut res: f64 = self.inner.individuals().par_iter()
                .map(|a| InnerEcologicalNiche::<I, M, R, S, T>::get_accumulated_resources(a.clone()) + 1.0)
                .sum();
            res += self.inner.resources().total();
            println!("Size: {} : Bytes: {} ; Fitness: {} ; Total Resources: {} ; Resources: {:?}",
                self.inner.population_size(),
                self.inner.population_mean_genome_size(),
                self.inner.population_mean_fitness(),
                res,
                self.inner.resources());
            // Save the population in regular intervalls with a timestamp and print some information.
            if start.elapsed() >= self.environment().population_save_intervall() {
                self.save_population();
                start = Instant::now();
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
    fn spawn_organism(inner: Arc<InnerEcologicalNiche<I, M, R, S, T>>, individual: Arc<Mutex<Individual<R, S, T>>>) {
        let tested = inner.testing(individual.clone());
        if tested {
            // Transcribe / translate the genome and test the organism.
            Self::add_fitness(individual.clone(), Self::test_organism(inner.clone(), individual.clone()));
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
    fn mate_organism(inner: Arc<InnerEcologicalNiche<I, M, R, S, T>>, individual: Arc<Mutex<Individual<R, S, T>>>) {
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
    fn test_organism(inner: Arc<InnerEcologicalNiche<I, M, R, S, T>>, individual: Arc<Mutex<Individual<R, S, T>>>) -> f64 {
        let mut organism = inner.load_organism(individual.clone());
        let mut organism_informations = Vec::new();
        // Repeatedly test the organism and supply all the testing information to the fitness
        // function.
        for _ in 0..inner.environment.testing_repetitions() {
            let (input, result_information) = (inner.supplier_function)();
            organism.set_input(input);
            let run_time = organism.live(&inner.environment);
            let output = organism.get_result();
            organism_informations.push(
                OrganismInformation::new(
                    output,
                    result_information,
                    inner.get_bytes(individual.clone()) * 8,
                    run_time,
                    *(&inner.environment.lifespan()),
                    inner.get_associated_inputs(individual.clone()),
                    inner.get_associated_outputs(individual.clone()),
                    organism.binary_size(),
                    *(&inner.environment.max_organism_size())
                )
            );
        }
        (inner.fitness_function)(organism_informations)
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
    fn add_fitness(individual: Arc<Mutex<Individual<R, S, T>>>, fitness: f64) {
        let mut ind = individual.lock()
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
    fn get_offspring(individual: Arc<Mutex<Individual<R, S, T>>>, inner: Arc<InnerEcologicalNiche<I, M, R, S, T>>) -> Vec<Individual<R, S, T>> {
        let mut offspring = Vec::new();
        // Use the accumulated resources to produce offspring.
        let number_of_offspring = inner.spend_resources_for_mating(individual.clone());
        for _ in 0..number_of_offspring {
            let partner = inner.get_random_genome();
            let ind = individual.lock()
                .expect("A thread paniced while holding the individual's lock.");
            offspring.push(ind.mate_and_mutate(partner, &inner.environment));
        }
        offspring
    }
}

struct InnerEcologicalNiche<I, M, R, S, T> {
    environment: Arc<Environment<M, R, S, T>>,
    population: Arc<Mutex<Population<R, S, T>>>,
    supplier_function: Box<dyn Fn() -> (Vec<T>, I) + Send + Sync + 'static>,
    fitness_function: Box<dyn Fn(Vec<OrganismInformation<I, T>>) -> f64 + Send + Sync + 'static>,
    phantom: PhantomData<M>,
}

impl<I, M: GenomeMutation<R, S, T>, R: Reaction<T>, S: State<T>, T: Information> InnerEcologicalNiche<I, M, R, S, T> {
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
    fn died(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> bool {
        let ind = individual.lock()
            .expect("A thread paniced while holding the individual's lock.");
        self.environment.death_chance(ind.age()) >= thread_rng().gen_range(0.0, 1.0)
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
    fn is_juvenil(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> bool {
        let ind = individual.lock()
            .expect("A thread paniced while holding the individual's lock.");
        self.environment.max_testing_age().map_or(true, |max_age| ind.age() < max_age)
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
    fn testing_by_chance(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> bool {
        let ind = individual.lock()
            .expect("A thread paniced while holding the individual's lock.");
        thread_rng().gen_range(0.0, 1.0) <= self.environment.testing_chance(ind.times_tested())
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
    fn testing(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> bool {
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
    fn get_uuid(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> Uuid {
        let ind = individual.lock()
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
    fn get_accumulated_resources(individual: Arc<Mutex<Individual<R, S, T>>>) -> f64 {
        let ind = individual.lock()
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
    fn spend_resources_for_mating(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> usize {
        individual.lock()
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
    fn get_bytes(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> usize {
        let ind = individual.lock()
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
    fn get_associated_inputs(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> usize {
        let ind = individual.lock()
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
    fn get_associated_outputs(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> usize {
        let ind = individual.lock()
            .expect("A thread paniced while holding the individual's lock.");
        ind.associated_outputs()
    }

    /// Returns the copy of a [`Genome`] of a random [`Individual`].
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock or the population is
    /// completely extinct.
    ///
    /// [`Individual`]: ./struct.Individual.html
    fn get_random_genome(&self) -> Genome<R, S, T> {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .random_genome()
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
        self.population.lock()
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
    fn individuals(&self) -> Vec<Arc<Mutex<Individual<R, S, T>>>> {
        self.population.lock()
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
        self.population.lock()
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
        self.population.lock()
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
        self.population.lock()
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
        self.population.lock()
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
        self.population.lock()
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
    fn remove_individual(&self, individual: Arc<Mutex<Individual<R, S, T>>>) {
        let uuid = self.get_uuid(individual.clone());
        // An individual consumes 1.0 resources when being born, so this has to be repatriated
        // additionally to the accumulated resources.
        let resources = Self::get_accumulated_resources(individual) + 1.0;
        {
            &self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .repatriate_resources(resources);
        }
        {
            &self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .remove(uuid)
            .expect("The individual could not be removed.");
        }
    }

    /// Recycles inavailable [`Resource`]s at the end of a generation.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    fn recycle(&self) {
        &self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .recycle();
    }

    /// Returns the [`Resource`]s.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    fn resources(&self) -> Resource {
        self.population.lock()
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
    fn load_organism(&self, individual: Arc<Mutex<Individual<R, S, T>>>) -> Organism<R, S, T> {
        individual.lock()
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
    fn append_population(&self, individuals: Vec<Individual<R, S, T>>) {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .append(individuals);
    }

}
