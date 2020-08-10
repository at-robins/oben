//! The `environment` module contains the setup of the evolutionary network.

extern crate bitvec;
extern crate rayon;
extern crate rand;

use super::binary::BinarySubstrate;
use std::path::{Path, PathBuf};
use uuid::{Uuid, v1::Context, v1::Timestamp};
use super::population::{Individual, Population, Organism, OrganismInformation};
use super::gene::Genome;
use super::resource::Resource;
use std::time::{Duration, Instant, SystemTime};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use rand::{thread_rng, Rng};

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
    /// The maximum offspring a [`Individual`] can produce per generation.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    max_offspring: usize,
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
            max_offspring: 10,
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
            max_offspring: self.max_offspring,
            death_age_sigmoid_midpoint: self.death_age_sigmoid_midpoint,
            lateral_gene_transfer_chance: self.lateral_gene_transfer_chance_or_default(),
            testing_chance_sigmoid_midpoint: self.testing_chance_sigmoid_midpoint,
            testing_repetitions: self.testing_repetitions,
            max_organism_size: self.max_organism_size
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

    /// Sets the maximum number of offspring an individual can produce per generation as specified.
    ///
    /// # Parameters
    ///
    /// * `offspring` - the limit of offspring per generation
    pub fn max_offspring(&mut self, offspring: usize) -> &mut Self {
        self.max_offspring = offspring;
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
#[derive(Debug, PartialEq, Clone)]
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
    /// The maximum number of offspring a [`Individual`] can produce per generation.
    ///
    /// [`Individual`]: ../population/struct.Individual.html
    max_offspring: usize,
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
        let s = self.testing_chance_sigmoid_midpoint() * 0.1;
        a / (1.0 + ((m - (age as f64)) / s).exp())
    }

    /// Returns the maximum offspring per individual per generation.
    pub fn max_offspring(&self) -> usize {
        self.max_offspring
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
        EnvironmentBuilder::new().build()
    }
}

impl<E> From<E> for Environment where E: AsRef<EnvironmentBuilder>{
    fn from(builder: E) -> Self {
        builder.as_ref().build()
    }
}

/// A `GlobalEnvironment` containing a [`Population`] and applying selective pressure.
///
/// [`Population`]: ../population/struct.Population.html
pub struct GlobalEnvironment<I> {
    inner: Arc<InnerGlobalEnvironment<I>>,
}

impl<I: 'static> GlobalEnvironment<I> {
    /// Creates a new `GlobalEnvironment` executing the evolutionary network by repeated
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
        environment: Environment,
        population: Population,
        supplier_function: Box<dyn Fn() -> (Vec<BinarySubstrate>, I)  + Send + Sync + 'static>,
        fitness_function: Box<dyn Fn(Vec<OrganismInformation<I>>) -> f64 + Send + Sync + 'static>) -> Self {
        GlobalEnvironment {
            inner: Arc::new(InnerGlobalEnvironment {
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
            let counter = Arc::new(Mutex::new(0u32));
            generation += 1;
            println!("Generation {}", generation);
            // Challenge the organisms in the population.
            self.inner.individuals().par_iter().for_each(|individual| {
                Self::spawn_organism(self.inner.clone(), individual.clone());
                *counter.lock().unwrap() += 1;
                if *counter.lock().unwrap() % 1000 == 0 {
                    println!("     Spawn Organism {}", *counter.lock().unwrap());
                }
            });
            // Distribute resources neccesarry for mating based on fitness.
            self.inner.distribute_resources();
            // Mate the organisms of the population and add offspring to the population.
            self.inner.individuals().par_iter().for_each(|individual| {
                Self::mate_organism(self.inner.clone(), individual.clone());
                *counter.lock().unwrap() += 1;
                if *counter.lock().unwrap() % 1000 == 0 {
                    println!("     Mate Organism {}", *counter.lock().unwrap());
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
            // Print statistics.
            println!("Size: {} : Bytes: {} ; Fitness: {}",
                self.inner.population_size(),
                self.inner.population_mean_genome_size(),
                self.inner.population_mean_fitness());
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
    fn spawn_organism(inner: Arc<InnerGlobalEnvironment<I>>, individual: Arc<Mutex<Individual>>) {
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
    fn mate_organism(inner: Arc<InnerGlobalEnvironment<I>>, individual: Arc<Mutex<Individual>>) {
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
    fn test_organism(inner: Arc<InnerGlobalEnvironment<I>>, individual: Arc<Mutex<Individual>>) -> f64 {
        let organism = inner.load_organism(individual.clone());
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
    fn add_fitness(individual: Arc<Mutex<Individual>>, fitness: f64) {
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
    fn get_offspring(individual: Arc<Mutex<Individual>>, inner: Arc<InnerGlobalEnvironment<I>>) -> Vec<Individual> {
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

struct InnerGlobalEnvironment<I> {
    environment: Arc<Environment>,
    population: Arc<Mutex<Population>>,
    supplier_function: Box<dyn Fn() -> (Vec<BinarySubstrate>, I) + Send + Sync + 'static>,
    fitness_function: Box<dyn Fn(Vec<OrganismInformation<I>>) -> f64 + Send + Sync + 'static>,
}

impl<I> InnerGlobalEnvironment<I> {
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
    fn died(&self, individual: Arc<Mutex<Individual>>) -> bool {
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
    fn is_juvenil(&self, individual: Arc<Mutex<Individual>>) -> bool {
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
    fn testing_by_chance(&self, individual: Arc<Mutex<Individual>>) -> bool {
        let ind = individual.lock()
            .expect("A thread paniced while holding the individual's lock.");
        thread_rng().gen_range(0.0, 1.0) <= self.environment.testing_chance(ind.age())
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
    fn testing(&self, individual: Arc<Mutex<Individual>>) -> bool {
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
    fn get_uuid(&self, individual: Arc<Mutex<Individual>>) -> Uuid {
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
    fn get_accumulated_resources(individual: Arc<Mutex<Individual>>) -> f64 {
        let ind = individual.lock()
            .expect("A thread paniced while holding the individual's lock.");
        // An individual consumes 1.0 resources when being born, so this has to be repatriated
        // additionally to the accumulated resources.
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
    fn spend_resources_for_mating(&self, individual: Arc<Mutex<Individual>>) -> usize {
        individual.lock()
            .expect("A thread paniced while holding the individual's lock.")
            .spend_resources_for_mating(&self.environment)
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
    fn get_bytes(&self, individual: Arc<Mutex<Individual>>) -> usize {
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
    fn get_associated_inputs(&self, individual: Arc<Mutex<Individual>>) -> usize {
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
    fn get_associated_outputs(&self, individual: Arc<Mutex<Individual>>) -> usize {
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
    fn get_random_genome(&self) -> Genome {
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
    fn individuals(&self) -> Vec<Arc<Mutex<Individual>>> {
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
    fn remove_individual(&self, individual: Arc<Mutex<Individual>>) {
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
    fn load_organism(&self, individual: Arc<Mutex<Individual>>) -> Organism {
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
    fn append_population(&self, individuals: Vec<Individual>) {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .append(individuals);
    }

}
