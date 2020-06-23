//! The `environment` module contains the setup of the evolutionary network.

extern crate bitvec;
extern crate rayon;

use bitvec::{boxed::BitBox, order::Local};
use std::path::{Path, PathBuf};
use uuid::{Uuid, v1::Context, v1::Timestamp};
use super::population::{ClonalPopulation, Population, Organism, OrganismInformation};
use std::time::{Duration, Instant, SystemTime};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

/// The sub-folder in which genome files are stored.
const SUBFOLDER_GENOME: &str = "genomes/dummy";
/// The sub-folder in which genome files of extinct populations are stored.
const SUBFOLDER_GENOME_EXTINCT: &str = "extinct/dummy";
/// The sub-folder in which population snapshot files are stored.
const SUBFOLDER_POPULATION: &str = "populations/dummy";
/// The file extension of genome files.
const FILE_EXTENSION_GENOME: &str = "genome";
/// The file extension of population files.
const FILE_EXTENSION_POPULATION: &str = "population";
/// The context for UUID creation.
const UUID_CONTEXT: Context = Context::new(0);

/// An `Environment` specifing settings for an evolutionary network to develop in.
#[derive(Debug, PartialEq, Clone)]
pub struct Environment {
    working_directory: PathBuf,
    /// The chance of a single offspring to carry a mutation.
    mutation_rate: f64,
    /// The maximum size a [`ClonalPopulation`] can grow to.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    population_size: u32,
    /// The relative size of clonal sub-populations that can still be detected.
    /// Any population with a smaller size will go extinct.
    extinction_threshold: f64,
    /// The amount of time an [`Organism`] of a [`ClonalPopulation`] has to complete a task.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    /// [`Organism`]: ../population/struct.Organism.html
    lifespan: Duration,
    /// The time intervall in which to save the current [`Population`] to a file.
    ///
    /// [`Population`]: ../population/struct.Population.html
    population_save_intervall: Duration,
    /// The node for UUID creation.
    uuid_node: [u8; 6],
    /// The maximum age of a [`ClonalPopulation`] until testing and fitness determination is
    /// performed.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    max_testing_age: Option<u32>,
    /// The size of a newly founded, mutated [`ClonalPopulation`].
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    clonal_population_founding_size: f64,
}

impl Environment {
    // TODO: Implement an environment builder.
    pub fn new(working_directory: PathBuf,
        mutation_rate: f64,
        population_size: u32,
        extinction_threshold: f64,
        lifespan: Duration,
        population_save_intervall: Duration,
        uuid_node: [u8; 6],
        max_testing_age: Option<u32>,
        clonal_population_founding_size: f64,) -> Self {
        Environment{working_directory,
            mutation_rate,
            population_size,
            extinction_threshold,
            lifespan,
            population_save_intervall,
            uuid_node,
            max_testing_age,
            clonal_population_founding_size}
    }

    /// Returns the path to the working directory.
    pub fn working_directory(&self) -> &Path {
        Path::new(&self.working_directory)
    }

    /// Returns the chance of a single offspring carrying a mutation.
    pub fn mutation_rate(&self) -> f64 {
        self.mutation_rate
    }

    /// Returns the size in individuals a [`Population`] can grow to.
    ///
    /// [`Population`]: ../population/struct.Population.html
    pub fn population_size(&self) -> u32 {
        self.population_size
    }

    /// Returns the threshold of relative populaion size which is still detectable before
    /// a [`ClonalPopulation`] goes extinct.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    pub fn extinction_threshold(&self) -> f64 {
        self.extinction_threshold
    }

    /// Returns the node for UUID creation.
    pub fn uuid_node(&self) -> &[u8; 6] {
        &self.uuid_node
    }

    /// Returns the maximum age until an individual is tested to determine the mean fitness.
    pub fn max_testing_age(&self) -> Option<u32> {
        self.max_testing_age
    }

    /// Returns the starting size of a newly found clonal population.
    pub fn clonal_population_founding_size(&self) -> f64 {
        self.clonal_population_founding_size
    }

    /// Returns the timestamp for UUID creation.
    ///
    /// # Panics
    ///
    /// If the system clock was set to earlier than Unix time start.
    pub fn uuid_timestamp(&self) -> Timestamp {
        let now = SystemTime::now().duration_since(SystemTime::UNIX_EPOCH)
            .expect("The system clock is not set correctly, so no UUIDs can be created.");
        Timestamp::from_unix(UUID_CONTEXT, now.as_secs(), now.subsec_nanos())
    }

    /// Generates a UUID based on the `Environment`.
    ///
    /// # Panics
    ///
    /// If the system clock was set to earlier than Unix time start.
    pub fn generate_uuid(&self) -> Uuid {
        // It was ensured, that the UUID node is exactly 6 byte long,
        // so the creation cannot fail.
        Uuid::new_v1(self.uuid_timestamp(), self.uuid_node()).unwrap()
    }

    /// Returns the file path to the [`Genome`] with the specified UUID.
    ///
    /// # Parameters
    ///
    /// * `genome_uuid` - the UUID of the [`Genome`]
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn genome_path(&self, genome_uuid: &Uuid) -> PathBuf {
        let mut path_to_genome: PathBuf = self.working_directory().into();
        path_to_genome.push(SUBFOLDER_GENOME);
        path_to_genome.set_file_name(genome_uuid.to_string());
        path_to_genome.set_extension(FILE_EXTENSION_GENOME);
        path_to_genome
    }

    /// Returns the file path to the extinct [`Genome`] with the specified UUID.
    ///
    /// # Parameters
    ///
    /// * `genome_uuid` - the UUID of the extinct [`Genome`]
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn extinct_genome_path(&self, genome_uuid: &Uuid) -> PathBuf {
        let mut path_to_genome: PathBuf = self.working_directory().into();
        path_to_genome.push(SUBFOLDER_GENOME_EXTINCT);
        path_to_genome.set_file_name(genome_uuid.to_string());
        path_to_genome.set_extension(FILE_EXTENSION_GENOME);
        path_to_genome
    }

    /// Returns the file path to the [`Population`] with the specified UUID.
    ///
    /// # Parameters
    ///
    /// * `population_uuid` - the UUID of the [`Population`]
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn population_path(&self, population_uuid: &Uuid) -> PathBuf {
        let mut path_to_genome: PathBuf = self.working_directory().into();
        path_to_genome.push(SUBFOLDER_POPULATION);
        path_to_genome.set_file_name(population_uuid.to_string());
        path_to_genome.set_extension(FILE_EXTENSION_POPULATION);
        path_to_genome
    }

    /// Returns the amount of time an [`Organism`] of a [`ClonalPopulation`] has to complete a task.
    ///
    /// [`Organism`]: ../population/struct.Organism.html
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    pub fn lifespan(&self) -> Duration {
        self.lifespan
    }

    /// Returns the intervall in which the whole [`Population`] is saved to a file.
    ///
    /// [`Population`]: ../population/struct.Population.html
    pub fn population_save_intervall(&self) -> Duration {
        self.population_save_intervall
    }

    /// Initialises the environment.
    ///
    /// # Panics
    ///
    /// If the initialisation process fails.
    fn initialise(&self) {
        self.create_subfolder(SUBFOLDER_GENOME);
        self.create_subfolder(SUBFOLDER_POPULATION);
        self.create_subfolder(SUBFOLDER_GENOME_EXTINCT);
    }

    /// Creates the specified sub-folder based on the current working directory.
    ///
    /// # Parameters
    ///
    /// * `subfolder` - the sub-folder to create
    ///
    /// # Panics
    ///
    /// If creation of the sub-folder fails.
    fn create_subfolder(&self, subfolder: &str) {
        let mut folder_path: PathBuf = self.working_directory().into();
        folder_path.push(subfolder);
        if let Err(err) = std::fs::create_dir_all(&folder_path) {
            // Without those folders the whole network cannot work correctly,
            // so panicing is the preferred option here instead of error handling.
            panic!("The folder {:?} could not be created: {}", folder_path, err);
        }
    }
}

impl Default for Environment {
    fn default() -> Self {
        Environment {
            working_directory: PathBuf::from("./working_directory"),
            mutation_rate: 0.001,
            population_size: 1_000_000,
            extinction_threshold: 1.0 / 10_000_000.0,
            lifespan: Duration::from_secs(1),
            population_save_intervall: Duration::from_secs(1800),
            uuid_node: rand::random(),
            max_testing_age: None,
            clonal_population_founding_size: 1.0 / 1_000_000.0,
        }
    }
}

pub struct GlobalEnvironment<I> {
    inner: Arc<InnerGlobalEnvironment<I>>,
}

impl<I: 'static> GlobalEnvironment<I> {
    pub fn new(environment: Environment,
        population: Population,
        supplier_function: Box<dyn Fn() -> (Vec<BitBox<Local, u8>>, I)  + Send + Sync + 'static>,
        fitness_function: Box<dyn Fn(Vec<Option<BitBox<Local, u8>>>, I, OrganismInformation) -> f64 + Send + Sync + 'static>) -> Self {
        GlobalEnvironment {
            inner: Arc::new(InnerGlobalEnvironment {
                state: Arc::new(Mutex::new(GlobalEnvironmentState::Execution)),
                environment: Arc::new(environment),
                population: Arc::new(Mutex::new(population)),
                supplier_function,
                fitness_function,
            })
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
        let mut start = Instant::now();
        loop {
            generation += 1;
            println!("Generation {}", generation);
            let clonal_populations = self.inner.clonal_populations();
            // Challenge the organisms in the population and add the offspring to the population.
            clonal_populations.par_iter().for_each(|clonal_population| {
                Self::spawn_organism(self.inner.clone(), clonal_population.clone());
            });
            // Calculate the new relative amount of individuals per clonal population.
            let total_population_size: f64 = clonal_populations.par_iter()
                .map(|clonal_population| self.inner.get_relative_size(clonal_population.clone()))
                .sum();
            clonal_populations.par_iter().for_each(|clonal_population| {
                    self.inner.adjust_size_to_reference(clonal_population.clone(), total_population_size);
                });
            // Remove extinct populations that are below the dection threshold.
            clonal_populations.par_iter()
                .filter(|clonal_population| self.inner.is_extinct((*clonal_population).clone()))
                .for_each(|clonal_population| {
                    self.inner.remove_clonal_population(clonal_population.clone());
                });
            // Save the population in regular intervalls with a timestamp and print some information.
            if start.elapsed() >= self.environment().population_save_intervall() {
                println!("Waiting to save...");
                self.inner.set_state(GlobalEnvironmentState::Saving);
                let population_id = self.environment().generate_uuid();
                let save_path = self.environment().population_path(&population_id);
                self.inner.save_population(save_path);
                println!("Saved!");
                start = Instant::now();
                self.inner.set_state(GlobalEnvironmentState::Execution);
            }
            // if generation > 300 {
            //     break;
            // }
        }
    }



    fn spawn_organism(inner: Arc<InnerGlobalEnvironment<I>>, clonal_population: Arc<Mutex<ClonalPopulation>>) {
        let fitness;
        if inner.is_juvenil(clonal_population.clone()) {
            // Transcribe / translate the genome and test the organism.
            let (input, result_information) = (inner.supplier_function)();
            let organism = inner.load_organism(clonal_population.clone());
            organism.set_input(input);
            let run_time = organism.live(&inner.environment);
            let output = organism.get_result();
            let oi = OrganismInformation::new(clonal_population.lock().unwrap().bytes(), run_time, *(&inner.environment.lifespan));
            fitness = (inner.fitness_function)(output, result_information, oi);
        } else {
            fitness = clonal_population.lock().unwrap().fitness().expect("The sub-population should have been tested before.");
        }
        let mutated_offspring = Self::get_mutated_offspring(clonal_population.clone(), fitness, inner.clone());
        // Add the mutated offspring to the popuation.
        inner.append_population(mutated_offspring);
    }

    fn get_mutated_offspring(clonal_population: Arc<Mutex<ClonalPopulation>>, fitness: f64, inner: Arc<InnerGlobalEnvironment<I>>) -> Vec<ClonalPopulation> {
        let mutated_offspring;
        if inner.is_juvenil (clonal_population.clone()) {
            let mut cp = clonal_population.lock().unwrap();
            mutated_offspring = cp.evaluate_new_fitness(fitness, &inner.environment);
        } else {
            let mut cp = clonal_population.lock().unwrap();
            mutated_offspring = cp.grow(&inner.environment);
        }
        mutated_offspring.into_iter()
            .map(|genome| {
                let uuid = inner.environment.generate_uuid();
                ClonalPopulation::found(uuid, genome, &inner.environment)
            }).collect()
    }

}

struct InnerGlobalEnvironment<I> {
    environment: Arc<Environment>,
    population: Arc<Mutex<Population>>,
    supplier_function: Box<dyn Fn() -> (Vec<BitBox<Local, u8>>, I) + Send + Sync + 'static>,
    fitness_function: Box<dyn Fn(Vec<Option<BitBox<Local, u8>>>, I, OrganismInformation) -> f64 + Send + Sync + 'static>,
    state: Arc<Mutex<GlobalEnvironmentState>>
}

impl<I> InnerGlobalEnvironment<I> {
    /// Checks if the specified [`ClonalPopulation`] is already extinct.
    /// A populaion cannot go extinct before its fitness was determined.
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`] to check for extinction
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn is_extinct(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> bool {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        !cp.has_fitness() && self.environment.extinction_threshold() > cp.relative_size()
    }

    /// Checks if the specified [`ClonalPopulation`] is juvenil and still needs testing.
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`] to check
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn is_juvenil(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> bool {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        self.environment.max_testing_age().map_or(true, |max_age| cp.age() < max_age)
    }

    /// Return the UUID of the specified [`ClonalPopulation`].
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn get_uuid(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> Uuid {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        *cp.uuid()
    }

    /// Return the relative size of the specified [`ClonalPopulation`].
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn get_relative_size(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> f64 {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        cp.relative_size()
    }

    /// Adjusts the relative size of this `ClonalPopulation` to the size of the whole reference
    /// [`Population`].
    ///
    /// # Parameters
    ///
    /// * `reference` - the reference population's size
    /// * `clonal_population` - the [`ClonalPopulation`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`Population`]: ../population/struct.Population.html
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn adjust_size_to_reference(&self, clonal_population: Arc<Mutex<ClonalPopulation>>, reference: f64) {
        let mut cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        cp.adjust_size_to_reference(reference);
    }

    /// Sets the [`GlobalEnvironmentState`] the [`GlobalEnvironment`] as specified.
    ///
    /// # Parameters
    ///
    /// * `state` - the new state to set
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the state lock.
    ///
    /// [`GlobalEnvironmentState`]: ./enum.GlobalEnvironmentState.html
    /// [`GlobalEnvironment`]: ./struct.GlobalEnvironment.html
    fn set_state(&self, state: GlobalEnvironmentState) {
        *self.state.lock().expect("Could not set the state of the environment.") = state;
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

    /// Returns all [`ClonalPopulation`]s.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn clonal_populations(&self) -> Vec<Arc<Mutex<ClonalPopulation>>> {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .clonal_populations()
    }

    /// Removes the specified [`ClonalPopulation`].
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`] to remove
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock or if the [`ClonalPopulation`]
    /// could not be removed.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn remove_clonal_population(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) {
        let uuid = self.get_uuid(clonal_population);
        &self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .remove(uuid)
            .expect("Clonal population could not be removed.");
    }

    /// Load the [`Organism`] corresponding to the specified [`ClonalPopulation`]
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population lock.
    ///
    /// [`Organism`]: ../population/struct.population.html
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn load_organism(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> Organism {
        clonal_population.lock()
            .expect("Another thread panicked while holding the clonal popuation lock.")
            .genome()
            .translate()
    }

    /// Add the [`ClonalPopulation`]s to the [`Population`].
    ///
    /// # Parameters
    ///
    /// * `clonal_populations` - the [`ClonalPopulation`]s to add
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population lock.
    ///
    /// [`Population`]: ../population/struct.Population.html
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn append_population(&self, clonal_populations: Vec<ClonalPopulation>) -> Vec<Arc<Mutex<ClonalPopulation>>>{
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .append(clonal_populations)
    }

}

/// The state the [`GlobalEnvironment`] is in.
///
/// [`GlobalEnvironment`]: ./struct.GlobalEnvironment.html
enum GlobalEnvironmentState {
    /// Execution of the network.
    Execution,
    /// Saving the network.
    Saving,
}
