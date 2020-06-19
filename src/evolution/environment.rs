//! The `environment` module contains the setup of the evolutionary network.

extern crate bitvec;
extern crate rayon;

use bitvec::{boxed::BitBox, order::Local};
use std::path::{Path, PathBuf};
use uuid::{Uuid, v1::Context, v1::Timestamp};
use super::population::{ClonalPopulation, Population};
use super::gene::Genome;
use std::time::{Duration, Instant, SystemTime};
use rayon::{ThreadPool, ThreadPoolBuilder, ThreadPoolBuildError};
use std::sync::{Arc, Mutex};
use std::collections::HashMap;

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
/// The number of threads to handle the computation of individuals.
const THREAD_NUMBER: usize = 10;
/// The time the main loop sleeps before updating.
const UPDATE_MAIN_LOOP: Duration = Duration::from_secs(5);
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
    max_clonal_population_size: u32,
    /// The rate in individualts / second at which individuals of a population die.
    death_rate: f64,
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
    uuid_node: [u8; 6]
}

impl Environment {
    // TODO: Implement an environment builder.
    pub fn new(working_directory: PathBuf,
        mutation_rate: f64,
        max_clonal_population_size: u32,
        death_rate: f64,
        lifespan: Duration,
        population_save_intervall: Duration,
        uuid_node: [u8; 6]) -> Self {
        Environment{working_directory,
            mutation_rate,
            max_clonal_population_size,
            death_rate,
            lifespan,
            population_save_intervall,
            uuid_node}
    }

    /// Returns the path to the working directory.
    pub fn working_directory(&self) -> &Path {
        Path::new(&self.working_directory)
    }

    /// Returns the chance of a single offspring carrying a mutation.
    pub fn mutation_rate(&self) -> f64 {
        self.mutation_rate
    }

    /// Returns the maximum size a [`ClonalPopulation`] can grow to.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    pub fn max_clonal_population_size(&self) -> u32 {
        self.max_clonal_population_size
    }

    /// Returns the rate in individualts / second at which individuals of a [`ClonalPopulation`] die.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    pub fn death_rate(&self) -> f64 {
        self.death_rate
    }

    /// Returns the node for UUID creation.
    pub fn uuid_node(&self) -> &[u8; 6] {
        &self.uuid_node
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
            max_clonal_population_size: 1_000_000,
            death_rate: 1.0,
            lifespan: Duration::from_secs(1),
            population_save_intervall: Duration::from_secs(1800),
            uuid_node: rand::random(),
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
        fitness_function: Box<dyn Fn(Vec<Option<BitBox<Local, u8>>>, I) -> f64 + Send + Sync + 'static>,
        write_verification: HashMap<Uuid, bool>) -> Self {
        GlobalEnvironment {
            inner: Arc::new(InnerGlobalEnvironment {
                state: Arc::new(Mutex::new(GlobalEnvironmentState::Execution)),
                environment: Arc::new(environment),
                population: Arc::new(Mutex::new(population)),
                supplier_function,
                fitness_function,
                pool: ThreadPoolBuilder::new().num_threads(THREAD_NUMBER).build().unwrap(),
                death_timer: Arc::new(Mutex::new(HashMap::new())),
                write_verification
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
        for clonal_population in self.inner.clonal_populations() {
            Self::spawn_organism(self.inner.clone(), clonal_population.clone());
        }
        println!("Starting the save loop...");
        // Save the population in regular intervalls with a timestamp and print some information.
        let mut start = Instant::now();
        loop {
            if start.elapsed() >= self.environment().population_save_intervall() {
                println!("Waiting to save...");
                self.inner.set_state(GlobalEnvironmentState::Saving);
                let population_id = self.environment().generate_uuid();
                let save_path = self.environment().population_path(&population_id);
                self.inner.save_population(save_path);
                    println!("Saved!");
                start = Instant::now();
                self.inner.set_state(GlobalEnvironmentState::Execution);
            } else {
                std::thread::sleep(UPDATE_MAIN_LOOP);
            }
        }
    }



    fn spawn_organism(inner: Arc<InnerGlobalEnvironment<I>>, clonal_population: Arc<Mutex<ClonalPopulation>>) {
        //TODO: Replace all unwraps.
        inner.clone().pool.spawn(move || {
            // Stop execution while the network is saving. This might cause catastrophic death
            // events and may need to be altered.
            while !inner.is_executing() {
                std::thread::sleep(Duration::from_millis(50));
            }
            // Calculate the deaths of the population and stop execution if it went extinct.
            if inner.calculate_deaths(clonal_population.clone()) {
                return
            }
            // Transcribe / translate the genome and test the organism.
            let (input, result_information) = (inner.supplier_function)();
            let organism_genome = inner.load_genome(clonal_population.clone());
            let organism = organism_genome.translate();
            organism.set_input(input);
            organism.live(&inner.environment);
            let output = organism.get_result();
            let fitness = (inner.fitness_function)(output, result_information);
            let mutated_offspring = Self::get_mutated_offspring(clonal_population.clone(), fitness, inner.environment.clone());
            // Restart this population.
            let b = inner.clone();
            inner.pool.spawn(move || Self::spawn_organism(b, clonal_population));
            // Add the mutated offspring to the popuation and start execution.
            let mutated_offspring = inner.append_population(mutated_offspring);
            for offspring in mutated_offspring.into_iter() {
                let b = inner.clone();
                inner.pool.spawn(move || Self::spawn_organism(b, offspring));
            }
        });
    }

    fn get_mutated_offspring(clonal_population: Arc<Mutex<ClonalPopulation>>, fitness: f64, environment: Arc<Environment>) -> Vec<ClonalPopulation> {
        // fn get_mutated_offspring(clonal_population: &Arc<Mutex<ClonalPopulation>>, fitness: f64, write_verification: Arc<Mutex<HashMap<Uuid, bool>>>, environment: &Environment) -> Vec<ClonalPopulation> {
        let mutated_offspring;
        {
            let mut cp = clonal_population.lock().unwrap();
            // while write_verification.lock().unwrap().get(cp.uuid()).is_none() {
            //     std::thread::sleep(Duration::from_millis(100));
            // }
            mutated_offspring = cp.evaluate_new_fitness(fitness, &environment);
        }
        mutated_offspring.iter()
            .map(|genome| {
                let uuid = environment.generate_uuid();
                if let Err(err) = genome.write_to_file(environment.genome_path(&uuid)) {
                    // Panicing on an error is fine, since writing genomes to files is an
                    // integral part of the network.
                    panic!("Writing mutated genome {} to a file filed: {}", &uuid, err);
                }
                // write_verification.lock().unwrap().insert(uuid, true);
                ClonalPopulation::found(uuid)
            }).collect()
    }

}

struct InnerGlobalEnvironment<I> {
    environment: Arc<Environment>,
    population: Arc<Mutex<Population>>,
    supplier_function: Box<dyn Fn() -> (Vec<BitBox<Local, u8>>, I) + Send + Sync + 'static>,
    fitness_function: Box<dyn Fn(Vec<Option<BitBox<Local, u8>>>, I) -> f64 + Send + Sync + 'static>,
    pool: ThreadPool,
    death_timer: Arc<Mutex<HashMap<Uuid, Instant>>>,
    write_verification: HashMap<Uuid, bool>,
    state: Arc<Mutex<GlobalEnvironmentState>>
}

impl<I> InnerGlobalEnvironment<I> {
    /// Checks whether the [`GlobalEnvironment`] is currently executing the network.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the state lock.
    ///
    /// [`GlobalEnvironment`]: ./struct.GlobalEnvironment.html
    fn is_executing(&self) -> bool {
        self.state.lock().expect("Could not read the state of the environment.").is_executing()
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
    /// * `clonal_population` - the UUID of the [`ClonalPopulation`] to remove
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock or if the [`ClonalPopulation`]
    /// could not be removed.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn remove_clonal_population(&self, clonal_population: Uuid) {
        &self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .remove(clonal_population, &self.environment)
            .expect("Clonal population could not be removed.");
    }

    /// Sets the current death event of a [`ClonalPopulation`] and returns the last one.
    ///
    /// # Parameters
    ///
    /// * `clone_id` - the UUID of the clonal population
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the death event lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn set_death_event(&self, clone_id: Uuid) -> Option<Instant> {
        self.death_timer.lock()
            .expect("A thread paniced while holding the death event lock.")
            .insert(clone_id, Instant::now())
    }

    /// Calculates the deathhs in the specified [`ClonalPopulation`] and returns if the population
    /// went extinct.
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`]
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the population lock, the death timer lock or if the [`ClonalPopulation`]
    /// could not be removed.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn calculate_deaths(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> bool{
        // Calculate the deaths of the population and remove it if it went extinct.
        let mut cp = clonal_population.lock().unwrap();
        if let Some(last_death_event) = self.set_death_event(cp.uuid().clone()) {
            cp.death_event(last_death_event.elapsed().as_secs_f64() / self.environment.death_rate());
        }
        if cp.is_extinct() {
            self.remove_clonal_population(*cp.uuid());
            // self.write_verification.lock().unwrap().remove(cp.uuid());
            true
        } else {
            false
        }
    }

    /// Load the [`Genome`] of the specified [`ClonalPopulation`]
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`] to load the genome from
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population lock or the [`Genome`]
    /// could not be loaded.
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn load_genome(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> Genome {
        let clonal_population_uuid = clonal_population.lock()
            .expect("Another thread panicked while holding the clonal popuation lock.")
            .uuid()
            .clone();
        //         while self.write_verification.lock().unwrap().get(&cp_uuid).is_none() {
        //             std::thread::sleep(Duration::from_millis(100));
        //         }
        match Genome::load_from_file(self.environment.genome_path(&clonal_population_uuid)) {
            Ok(genome) => genome,
            Err(err) => panic!("Genome {} could not be loaded: {}", clonal_population_uuid, err),
        }
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

impl GlobalEnvironmentState {
    /// Checks whether the [`GlobalEnvironment`] is currently executing the network.
    ///
    /// [`GlobalEnvironment`]: ./struct.GlobalEnvironment.html
    fn is_executing(&self) -> bool {
        match self {
            GlobalEnvironmentState::Execution => true,
            _ => false
        }
    }
}
