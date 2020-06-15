//! The `environment` module contains the setup of the evolutionary network.

extern crate bitvec;
extern crate rayon;

use bitvec::boxed::BitBox;
use std::path::{Path, PathBuf};
use uuid::{Uuid, v1::Context, v1::Timestamp};
use super::population::Population;
use std::time::{Duration, Instant, SystemTime};
use rayon::{ThreadPool};
use std::sync::{Arc, Mutex};

/// The sub-folder in which genome files are stored.
const SUBFOLDER_GENOME: &str = "genomes";
/// The sub-folder in which genome files of extinct populations are stored.
const SUBFOLDER_GENOME_EXTINCT: &str = "extinct";
/// The sub-folder in which population snapshot files are stored.
const SUBFOLDER_POPULATION: &str = "populations";
/// The file extension of genome files.
const FILE_EXTENSION_GENOME: &str = "genome";
/// The file extension of population files.
const FILE_EXTENSION_POPULATION: &str = "population";
/// The number of threads to handle the computation of individuals.
const THREAD_NUMBER: usize = 1000;
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
    /// The amount of time an [`Individual`] of a [`ClonalPopulation`] has to complete a task.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    /// [`Individual`]: ../population/struct.Individual.html
    lifespan: Duration,
    /// The time intervall in which to save the current [`Population`] to a file.
    ///
    /// [`Population`]: ../population/struct.Population.html
    population_save_intervall: Duration,
    /// The node for UUID creation.
    uuid_node: [u8; 6]
}

impl Environment {
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
    pub fn genome_path(&self, genome_uuid: Uuid) -> PathBuf {
        let mut path_to_genome: PathBuf = self.working_directory().into();
        path_to_genome.push(SUBFOLDER_GENOME);
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
    pub fn population_path(&self, population_uuid: Uuid) -> PathBuf {
        let mut path_to_genome: PathBuf = self.working_directory().into();
        path_to_genome.push(SUBFOLDER_POPULATION);
        path_to_genome.set_file_name(population_uuid.to_string());
        path_to_genome.set_extension(FILE_EXTENSION_POPULATION);
        path_to_genome
    }

    /// Returns the amount of time an [`Individual`] of a [`ClonalPopulation`] has to complete a task.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
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

pub struct GlobalEnvironment {
    environment: Environment,
    population: Arc<Mutex<Population>>,
    supplier_function: Box<dyn Fn() -> Vec<BitBox>>,
    fitness_function: Box<dyn Fn(Vec<Option<BitBox>>) -> f64>,
    pool: Arc<Mutex<ThreadPool>>,
}

impl GlobalEnvironment {
    pub fn breath_life(&mut self) {
        self.environment.initialise();
        let pool = self.pool.lock().expect("A thread panicked while having access to the complete threadpool.");
        //for clonal_population in
        // Save the population in regular intervalls with a timestamp.
        let mut start = Instant::now();
        loop {
            if start.elapsed() >= self.environment.population_save_intervall() {
                let population_id = self.environment.generate_uuid();
                let save_path = self.environment.population_path(population_id);
                self.population.lock()
                    .expect("A thread paniced while holding the population lock.")
                    .snapshot_to_file(&save_path)
                    .expect(&format!("The file {:?} could not be created.", save_path));
                start = Instant::now();
            } else {
                std::thread::sleep(UPDATE_MAIN_LOOP);
            }
        }
        // TODO: Push all ClonalPopulations to threadpool
    }
}
