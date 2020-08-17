//! The `configuration` module contains the configuration setup of the evolutionary network.

extern crate rand;

use rand::{thread_rng, Rng};
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};
use super::super::resource::Resource;
use uuid::{Uuid, v1::Context, v1::Timestamp};

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

/// An `EnvironmentBuilder` specifing settings for an evolutionary network to develop in and
/// returning the corresponding [`Environment`].
///
/// [`Environment`]: ./struct.Environment.html
#[derive(Debug, PartialEq, Clone)]
pub struct EnvironmentBuilder {
    working_directory: PathBuf,
    /// The chance of a single offspring to carry a mutation.
    mutation_rate: Option<f64>,
    /// The maximum size a [`Individual`] can grow to.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    population_size: u32,
    /// The half life time of resources that are released by the death of individuals
    /// before becomming available again.
    resource_half_life: f64,
    /// The amount of time an [`Organism`] of a [`Individual`] has to complete a task.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Organism`]: ../population/struct.Organism.html
    lifespan: Duration,
    /// The time intervall in which to save the current [`Population`] to a file.
    ///
    /// [`Population`]: ../population/struct.Population.html
    population_save_intervall: Duration,
    /// The node for UUID creation.
    uuid_node: [u8; 6],
    /// The maximum age of a [`Individual`] until testing and fitness determination is
    /// always performed. After this age testing is performed by chance.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    max_testing_age: Option<u32>,
    /// The midpoint of the sigmoid determining the chance of death depending on the age of a
    /// [`Individual`].
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    death_age_sigmoid_midpoint: f64,
    /// The chance per growth event that a lateral gene transfer between two
    /// sub-populations happens.
    lateral_gene_transfer_chance: Option<f64>,
    /// The 50 percent midpoint of test cycles of the chance determining sigmoid for testing a [`Individual`].
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    testing_chance_sigmoid_midpoint: f64,
    /// The number of times a [`Individual`] is tested per test cycle.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    testing_repetitions: u32,
    /// The maximum size a single [`Organism`] can grow to during testing in bit.
    ///
    /// [`Organism`]: ../population/struct.Organism.html
    max_organism_size: usize,
}

impl EnvironmentBuilder {

    pub fn new() -> Self {
        EnvironmentBuilder {
            working_directory: PathBuf::from("./working_directory"),
            mutation_rate: None,
            population_size: 1_000_000,
            resource_half_life: 3.0,
            lifespan: Duration::from_secs(1),
            population_save_intervall: Duration::from_secs(1800),
            uuid_node: rand::random(),
            max_testing_age: None,
            death_age_sigmoid_midpoint: 3.0,
            lateral_gene_transfer_chance: None,
            testing_chance_sigmoid_midpoint: 50.0,
            testing_repetitions: 1,
            max_organism_size: 8 * 1024 * 1024 * 50
        }
    }

    /// Build an [`Environment`] to specify run time properties of an evolutionary network.
    ///
    /// [`Environment`]: ./struct.Environment.html
    pub fn build(&self) -> Environment {
        Environment {
            working_directory: self.working_directory.clone(),
            mutation_rate: self.mutation_rate_or_default(),
            population_size: self.population_size,
            resource_half_life: self.resource_half_life,
            lifespan: self.lifespan,
            population_save_intervall: self.population_save_intervall,
            uuid_node: self.uuid_node,
            max_testing_age: self.max_testing_age,
            death_age_sigmoid_midpoint: self.death_age_sigmoid_midpoint,
            lateral_gene_transfer_chance: self.lateral_gene_transfer_chance_or_default(),
            testing_chance_sigmoid_midpoint: self.testing_chance_sigmoid_midpoint,
            testing_repetitions: self.testing_repetitions,
            max_organism_size: self.max_organism_size,
            uuid_context: Context::new(0),
        }
    }

    /// Sets the working directory as specified.
    ///
    /// # Parameters
    ///
    /// * `working_directory` - the working directory
    pub fn working_directory<P: AsRef<Path>>(&mut self, working_directory: P) -> &mut Self {
        self.working_directory = working_directory.as_ref().into();
        self
    }

    /// Sets the mutation rate as specified.
    ///
    /// # Parameters
    ///
    /// * `mutation_rate` - the chance of a single mutation occuring
    pub fn mutation_rate(&mut self, mutation_rate: f64) -> &mut Self {
        self.mutation_rate = Some(mutation_rate);
        self
    }

    /// Sets the population size as specified.
    ///
    /// # Parameters
    ///
    /// * `population_size` - the maximum size of the total population
    pub fn population_size(&mut self, population_size: u32) -> &mut Self {
        self.population_size = population_size;
        self
    }

    /// Sets the half life time of [`Resource`]s before becomming available again.
    ///
    /// # Parameters
    ///
    /// * `half_life` - the half life time
    ///
    /// ../resource/struct.Resource.html
    pub fn resource_half_life(&mut self, half_life: f64) -> &mut Self {
        self.resource_half_life = half_life;
        self
    }

    /// Sets the maximum run time of a sub-population per task as specified.
    ///
    /// # Parameters
    ///
    /// * `lifespan` - the maximum run time per tested task
    pub fn lifespan(&mut self, lifespan: Duration) -> &mut Self {
        self.lifespan = lifespan;
        self
    }

    /// Sets the intervall in which the total population is saved to disk as specified.
    /// The actual saving intervall might be longer as saving is only performed after each
    /// generation.
    ///
    /// # Parameters
    ///
    /// * `population_save_intervall` - the save intervall
    pub fn population_save_intervall(&mut self, population_save_intervall: Duration) -> &mut Self {
        self.population_save_intervall = population_save_intervall;
        self
    }

    /// Sets the UUID node as specified.
    ///
    /// # Parameters
    ///
    /// * `uuid_node` - the node for UUID generation
    pub fn uuid_node(&mut self, uuid_node: [u8; 6]) -> &mut Self {
        self.uuid_node = uuid_node;
        self
    }

    /// Sets the number of times testing of a sub-population is mandatorily performed as specified.
    /// After this threshold is reached, testing is performed on a statistical basis.
    /// If `None` every sub-population is tested in every generation.
    ///
    /// # Parameters
    ///
    /// * `max_testing_age` - the maximum mandatory testing age
    pub fn max_testing_age(&mut self, max_testing_age: Option<u32>) -> &mut Self {
        self.max_testing_age = max_testing_age;
        self
    }

    /// Sets the maximum size a single [`Organism`] can grow to during testing in bit.
    ///
    /// # Parameters
    ///
    /// * `max_testing_age` - the maximum size
    ///
    /// [`Organism`]: ../population/struct.Organism.html
    pub fn max_organism_size(&mut self, max_organism_size: usize) -> &mut Self {
        self.max_organism_size = max_organism_size;
        self
    }

    /// Sets the midpoint of the sigmoid determining the chance of death for an individual
    /// depending on its age.
    ///
    /// # Parameters
    ///
    /// * `death_age_sigmoid_midpoint` - the midpoint of the sigmoid
    pub fn death_age_sigmoid_midpoint(&mut self, death_age_sigmoid_midpoint: f64) -> &mut Self {
        self.death_age_sigmoid_midpoint = death_age_sigmoid_midpoint;
        self
    }

    /// Sets the chance of lateral gene transfer as specified.
    ///
    /// # Parameters
    ///
    /// * `lateral_gene_transfer_chance` - the chance of lateral gene transfer
    pub fn lateral_gene_transfer_chance(&mut self, lateral_gene_transfer_chance: f64) -> &mut Self {
        self.lateral_gene_transfer_chance = Some(lateral_gene_transfer_chance);
        self
    }

    /// Sets the midpoint of the sigmoid determining the chance of statistical testing as specified.
    /// This is the number of tests for a sub-populations where the chance of testing
    /// per generation is 50 percent.
    ///
    /// # Parameters
    ///
    /// * `testing_chance_sigmoid_midpoint` - the midpoint of the testing chance sigmoid
    pub fn testing_chance_sigmoid_midpoint(&mut self, testing_chance_sigmoid_midpoint: f64) -> &mut Self {
        self.testing_chance_sigmoid_midpoint = testing_chance_sigmoid_midpoint;
        self
    }

    /// Sets the number of repetitions per testing cycle. The fitness function is called
    /// that many times for evalution.
    ///
    /// # Parameters
    ///
    /// * `testing_repetitions` - the number of repetitions per testing cycle
    pub fn testing_repetitions(&mut self, testing_repetitions: u32) -> &mut Self {
        self.testing_repetitions = testing_repetitions;
        self
    }

    /// Returns the mutation rate if set. Otherwise defaults to a population size dependent value.
    fn mutation_rate_or_default(&self) -> f64 {
        if let Some(mutation_rate) = self.mutation_rate {
            mutation_rate
        } else {
            // Defaults to around 1000 mutations per generation.
            1000.0 / self.population_size as f64
        }
    }

    /// Returns the chance of lateral gene transfer if set.
    /// Otherwise defaults to a population size dependent value.
    fn lateral_gene_transfer_chance_or_default(&self) -> f64 {
        if let Some(lateral_gene_transfer_chance) = self.lateral_gene_transfer_chance {
            lateral_gene_transfer_chance
        } else {
            // Defaults to around 10 events per generation.
            10.0 / self.population_size as f64
        }
    }
}

impl Default for EnvironmentBuilder {
    fn default() -> Self {
        EnvironmentBuilder::new()
    }
}

/// An `Environment` specifing settings for an evolutionary network to develop in.
#[derive(Debug)]
pub struct Environment {
    working_directory: PathBuf,
    /// The chance of a single offspring to carry a mutation.
    mutation_rate: f64,
    /// The maximum size a [`Individual`] can grow to.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    population_size: u32,
    /// The half life time of resources that are released by the death of individuals
    /// before becomming available again.
    resource_half_life: f64,
    /// The amount of time an [`Organism`] of a [`Individual`] has to complete a task.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Organism`]: ../population/struct.Organism.html
    lifespan: Duration,
    /// The time intervall in which to save the current [`Population`] to a file.
    ///
    /// [`Population`]: ../population/struct.Population.html
    population_save_intervall: Duration,
    /// The node for UUID creation.
    uuid_node: [u8; 6],
    /// The maximum age of a [`Individual`] until testing and fitness determination is
    /// always performed. After this age testing is performed by chance.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    max_testing_age: Option<u32>,
    /// The midpoint of the sigmoid determining the chance of death depending on the age of a
    /// [`Individual`].
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    death_age_sigmoid_midpoint: f64,
    /// The chance per growth event that a lateral gene transfer between two
    /// sub-populations happens.
    lateral_gene_transfer_chance: f64,
    /// The 50 percent midpoint of test cycles of the chance determining sigmoid for testing a [`Individual`].
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    testing_chance_sigmoid_midpoint: f64,
    /// The number of times a [`Individual`] is tested per test cycle.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    testing_repetitions: u32,
    /// The maximum size a single [`Organism`] can grow to during testing in bit.
    ///
    /// [`Organism`]: ../population/struct.Organism.html
    max_organism_size: usize,
    /// The context for UUID creation.
    uuid_context: Context,
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

    /// Returns the size in individuals a [`Population`] can grow to.
    ///
    /// [`Population`]: ../population/struct.Population.html
    pub fn population_size(&self) -> u32 {
        self.population_size
    }

    /// Returns the half life time of [`Resource`]s before becomming available again.
    /// [`Resource`]s are consumed by generation of offspring and released by death of
    /// [`Individual`]s.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn resource_half_life(&self) -> f64 {
        self.resource_half_life
    }

    /// Creates [`Resource`]s based on population size and half life time.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn generate_resources(&self) -> Resource {
        Resource::new(self.population_size() as f64, self.resource_half_life())
    }

    /// Returns the node for UUID creation.
    pub fn uuid_node(&self) -> &[u8; 6] {
        &self.uuid_node
    }

    /// Returns the maximum age until an individual is tested to determine the mean fitness.
    pub fn max_testing_age(&self) -> Option<u32> {
        self.max_testing_age
    }

    /// Returns the maximum size an individual can grow to during testing in bit.
    pub fn max_organism_size(&self) -> usize {
        self.max_organism_size
    }

    /// The midpoint of the testing chance sigmoid in test cycles.
    pub fn testing_chance_sigmoid_midpoint(&self) -> f64 {
        self.testing_chance_sigmoid_midpoint
    }

    /// Returns the chance for testing a sub-population.
    ///
    /// # Parameters
    ///
    /// * `age` - the number of times the sub-population was already tested
    pub fn testing_chance(&self, age: u32) -> f64 {
        // Asymptote.
        let a = 1.0;
        // Midpoint.
        let m = self.testing_chance_sigmoid_midpoint();
        // Steepness.
        let s = self.testing_chance_sigmoid_midpoint() * 0.2;
        a / (1.0 + (((age as f64) - m) / s).exp())
    }

    /// Returns the chance of death for an inividual of the specified age.
    ///
    /// # Parameters
    ///
    /// * `age` - the number of times the sub-population was already tested
    pub fn death_chance(&self, age: u32) -> f64 {
        // Asymptote.
        let a = 1.0;
        // Midpoint.
        let m = self.death_age_sigmoid_midpoint();
        // Steepness.
        let s = self.testing_chance_sigmoid_midpoint() * 0.3;
        a / (1.0 + ((m - (age as f64)) / s).exp())
    }

    /// Returns the midpoint of the sigmoid determing the chance of death depending on the age
    /// of an individual.
    pub fn death_age_sigmoid_midpoint(&self) -> f64 {
        self.death_age_sigmoid_midpoint
    }

    /// Checks if lateral gene transfer happend on a statistical basis.
    pub fn lateral_gene_transfer(&self) -> bool {
        thread_rng().gen_range(0.0, 1.0) <= self.lateral_gene_transfer_chance
    }

    /// Returns the number of repetitions per testing cycle.
    pub fn testing_repetitions(&self) -> u32 {
        self.testing_repetitions
    }

    /// Returns the timestamp for UUID creation.
    ///
    /// # Panics
    ///
    /// If the system clock was set to earlier than Unix time start.
    pub fn uuid_timestamp(&self) -> Timestamp {
        let now = SystemTime::now().duration_since(SystemTime::UNIX_EPOCH)
            .expect("The system clock is not set correctly, so no UUIDs can be created.");
        Timestamp::from_unix(&self.uuid_context, now.as_secs(), now.subsec_nanos())
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

    /// Returns the amount of time an [`Organism`] of a [`Individual`] has to complete a task.
    ///
    /// [`Organism`]: ../population/struct.Organism.html
    /// [`Individual`]: ../population/struct.Individual.html
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
    pub fn initialise(&self) {
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
        EnvironmentBuilder::new().build()
    }
}

impl<E> From<E> for Environment where E: AsRef<EnvironmentBuilder>{
    fn from(builder: E) -> Self {
        builder.as_ref().build()
    }
}
