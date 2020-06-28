//! The `environment` module contains the setup of the evolutionary network.

extern crate bitvec;
extern crate rayon;
extern crate rand;

use bitvec::{boxed::BitBox, order::Local};
use std::path::{Path, PathBuf};
use uuid::{Uuid, v1::Context, v1::Timestamp};
use super::population::{ClonalPopulation, Population, Organism, OrganismInformation};
use super::gene::{Gene, GenomeMutation};
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
    /// The maximum size a [`ClonalPopulation`] can grow to.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    population_size: u32,
    /// The relative size of clonal sub-populations that can still be detected.
    /// Any population with a smaller size will go extinct.
    extinction_threshold: Option<f64>,
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
    /// always performed. After this age testing is performed by chance.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    max_testing_age: Option<u32>,
    /// The size of a newly founded, mutated [`ClonalPopulation`].
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    clonal_population_founding_size: Option<f64>,
    /// The standard deviation of the fitness dependent growth of a [`ClonalPopulation`].
    /// It is specified in percentage of the sub-population's size.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    clonal_population_growth_sd: f64,
    /// The chance per growth event that a lateral gene transfer between two
    /// sub-populations happens.
    lateral_gene_transfer_chance: Option<f64>,
    /// The 50 percent midpoint of test cycles of the chance determining sigmoid for testing a [`ClonalPopulation`].
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    testing_chance_sigmoid_midpoint: f64,
    /// The number of times a [`ClonalPopulation`] is tested per test cycle.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
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
            extinction_threshold: None,
            lifespan: Duration::from_secs(1),
            population_save_intervall: Duration::from_secs(1800),
            uuid_node: rand::random(),
            max_testing_age: None,
            clonal_population_founding_size: None,
            clonal_population_growth_sd: 0.1,
            lateral_gene_transfer_chance: None,
            testing_chance_sigmoid_midpoint: 50.0,
            testing_repetitions: 1,
            max_organism_size: 8 * 1024 * 1024 * 50
        }
    }

    /// Build an [`Environment`] to specify run time properties of an evolutionary network.
    pub fn build(&self) -> Environment {
        Environment {
            working_directory: self.working_directory.clone(),
            mutation_rate: self.mutation_rate_or_default(),
            population_size: self.population_size,
            extinction_threshold: self.extinction_threshold_or_default(),
            lifespan: self.lifespan,
            population_save_intervall: self.population_save_intervall,
            uuid_node: self.uuid_node,
            max_testing_age: self.max_testing_age,
            clonal_population_founding_size: self.clonal_population_founding_size_or_default(),
            clonal_population_growth_sd: self.clonal_population_growth_sd,
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

    /// Sets the extinction thresholdas as specified.
    ///
    /// # Parameters
    ///
    /// * `extinction_threshold` - the minimum relative size of a sub-population before being
    /// no longer detectable and thereby going extinct
    pub fn extinction_threshold(&mut self, extinction_threshold: f64) -> &mut Self {
        self.extinction_threshold = Some(extinction_threshold);
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

    /// Sets the relative size of newly found populations as specified.
    ///
    /// # Parameters
    ///
    /// * `clonal_population_founding_size` - the initial size of sub-populations
    pub fn clonal_population_founding_size(&mut self, clonal_population_founding_size: f64) -> &mut Self {
        self.clonal_population_founding_size = Some(clonal_population_founding_size);
        self
    }

    /// Sets the standard deviation for growth of sub-populations as specified.
    /// The standard deviation is specified in percentage of relative sub-population size.
    ///
    /// # Parameters
    ///
    /// * `clonal_population_growth_sd` - the growth standard deviation
    pub fn clonal_population_growth_sd(&mut self, clonal_population_growth_sd: f64) -> &mut Self {
        self.clonal_population_growth_sd = clonal_population_growth_sd;
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

    /// Returns the extinction threshold if set. Otherwise defaults to a population size dependent value.
    fn extinction_threshold_or_default(&self) -> f64 {
        if let Some(extinction_threshold) = self.extinction_threshold {
            extinction_threshold
        } else {
            // Defaults to a bit less than one individual.
            0.3 / self.population_size as f64
        }
    }

    /// Returns the clonal population founding size if set.
    /// Otherwise defaults to a population size dependent value.
    fn clonal_population_founding_size_or_default(&self) -> f64 {
        if let Some(clonal_population_founding_size) = self.clonal_population_founding_size {
            clonal_population_founding_size
        } else {
            // Defaults to a bit more than one individual.
            3.0 / self.population_size as f64
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
    /// always performed. After this age testing is performed by chance.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    max_testing_age: Option<u32>,
    /// The size of a newly founded, mutated [`ClonalPopulation`].
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    clonal_population_founding_size: f64,
    /// The standard deviation of the fitness dependent growth of a [`ClonalPopulation`].
    /// It is specified in percentage of the sub-population's size.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    clonal_population_growth_sd: f64,
    /// The chance per growth event that a lateral gene transfer between two
    /// sub-populations happens.
    lateral_gene_transfer_chance: f64,
    /// The 50 percent midpoint of test cycles of the chance determining sigmoid for testing a [`ClonalPopulation`].
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    testing_chance_sigmoid_midpoint: f64,
    /// The number of times a [`ClonalPopulation`] is tested per test cycle.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
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

    /// Returns the threshold of relative population size which is still detectable before
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
    /// * `size` - the relative size of the sub-population
    /// * `age` - the number of times the sub-population was already tested
    pub fn testing_chance(&self, size: f64, age: u32) -> f64 {
        // Asymptote.
        let a = 1.0;
        // Midpoint.
        let m = self.testing_chance_sigmoid_midpoint();
        // Steepness.
        let s = self.testing_chance_sigmoid_midpoint() * 0.2;
        0.5 * size + 0.5 * (a / (1.0 + (((age as f64) - m) / s).exp()))
    }

    /// Returns the starting size of a newly found clonal population.
    pub fn clonal_population_founding_size(&self) -> f64 {
        self.clonal_population_founding_size
    }

    /// Returns the growth standard deviation of clonal populations in percentage of their size.
    pub fn clonal_population_growth_sd(&self) -> f64 {
        self.clonal_population_growth_sd
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
        EnvironmentBuilder::new().build()
    }
}

impl<E> From<E> for Environment where E: AsRef<EnvironmentBuilder>{
    fn from(builder: E) -> Self {
        builder.as_ref().build()
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
            // Challenge the organisms in the population and add the offspring to the population.
            self.inner.clonal_populations().par_iter().for_each(|clonal_population| {
                Self::spawn_organism(self.inner.clone(), clonal_population.clone());
                *counter.lock().unwrap() += 1;
                if *counter.lock().unwrap() % 1000 == 0 {
                    println!("     Organism {}", *counter.lock().unwrap());
                }
            });
            // Calculate the new relative amount of individuals per clonal population.
            let clonal_populations = self.inner.clonal_populations();
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
            let mut bytes = 0usize;
            let mut fitness = 0.0;
            for cp in self.inner.population.lock().unwrap().clonal_populations() {
                let c = cp.lock().unwrap();
                bytes += c.bytes();
                fitness += c.fitness().unwrap_or(0.0);
            }
            fitness /= self.inner.population.lock().unwrap().clonal_populations().len() as f64;
            println!("Size: {} : Bytes: {} ; Fitness: {}", self.inner.population.lock().unwrap().clonal_populations().len(), bytes, fitness);
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

    fn spawn_organism(inner: Arc<InnerGlobalEnvironment<I>>, clonal_population: Arc<Mutex<ClonalPopulation>>) {
        let tested = inner.testing(clonal_population.clone());
        if tested {
            // Transcribe / translate the genome and test the organism.
            let organism = inner.load_organism(clonal_population.clone());
            for _ in 0..inner.environment.testing_repetitions() {
                Self::add_fitness(clonal_population.clone(), Self::test_organism(inner.clone(), &organism, clonal_population.clone()));
            }
        }
        let mutated_offspring = Self::get_mutated_offspring(clonal_population.clone(), inner.clone());
        // Add the mutated offspring to the population.
        inner.append_population(mutated_offspring);
    }

    /// Tests the [`Organism`] and returns the evaluated fitness.
    ///
    /// # Parameters
    ///
    /// * `inner` - the [`Environment`] the [`Organism`] is living in
    /// * `organism` - the [`Organism`] to test
    /// * `clonal_population` - the [`ClonalPopulation`] the [`Organism`] comes from
    ///
    /// [`Environment`]: ./struct.Environment.html
    /// [`Organism`]: ../population/struct.Organism.html
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn test_organism(inner: Arc<InnerGlobalEnvironment<I>>, organism: &Organism, clonal_population: Arc<Mutex<ClonalPopulation>>) -> f64 {
        let (input, result_information) = (inner.supplier_function)();
        organism.set_input(input);
        let run_time = organism.live(&inner.environment);
        let output = organism.get_result();
        let oi = OrganismInformation::new(
            inner.get_bytes(clonal_population.clone()) * 8,
            run_time,
            *(&inner.environment.lifespan()),
            inner.get_associated_inputs(clonal_population.clone()),
            organism.binary_size(),
            *(&inner.environment.max_organism_size())
        );
        (inner.fitness_function)(output, result_information, oi)
    }

    /// Adds the specified fitness to the specified [`ClonalPopulation`].
    ///
    /// # Parameters
    ///
    /// * `clonal_population` - the [`ClonalPopulation`]
    /// * `fitness` - the fitness to add
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    fn add_fitness(clonal_population: Arc<Mutex<ClonalPopulation>>, fitness: f64) {
        let mut cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        cp.evaluate_new_fitness(fitness)
    }

    fn get_mutated_offspring(clonal_population: Arc<Mutex<ClonalPopulation>>, inner: Arc<InnerGlobalEnvironment<I>>) -> Vec<ClonalPopulation> {
        let mut mutated_offspring;
        {
            // This must be scoped since lateral
            // gene transfer from the same clonal population is possible. This specifically
            // would result into a dead lock.
            let mut cp = clonal_population.lock()
                .expect("A thread paniced while holding the clonal population's lock.");
            mutated_offspring = cp.grow(&inner.environment);
        }
        // Perform lateral gene transfer by chance.
        if inner.environment.lateral_gene_transfer() {
            // This must be called before aquiring the clonal population's lock since lateral
            // gene transfer from the same clonal population is possible. This specifically
            // would result into a dead lock.
            let gene = inner.get_random_gene();
            let cp = clonal_population.lock().unwrap();
            if let Some(mutated_genome) = GenomeMutation::lateral_gene_transfer(cp.genome(), gene) {
                mutated_offspring.push(mutated_genome);
            }
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
}

impl<I> InnerGlobalEnvironment<I> {
    /// Checks if the specified [`ClonalPopulation`] is already extinct.
    /// A population cannot go extinct before its fitness was determined.
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
        cp.has_fitness() && self.environment.extinction_threshold() > cp.relative_size()
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

    /// Checks if the specified [`ClonalPopulation`] should be testedby chance.
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
    fn testing_by_chance(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> bool {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        thread_rng().gen_range(0.0, 1.0) <= self.environment.testing_chance(cp.relative_size(), cp.age())
    }

    /// Checks if the specified [`ClonalPopulation`] should be tested.
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
    fn testing(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> bool {
        self.is_juvenil(clonal_population.clone()) || self.testing_by_chance(clonal_population)
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

    /// Return the size in bytes of the specified [`ClonalPopulation`].
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
    fn get_bytes(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> usize {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        cp.bytes()
    }

    /// Return the number of associated inputs for the specified [`ClonalPopulation`].
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
    fn get_associated_inputs(&self, clonal_population: Arc<Mutex<ClonalPopulation>>) -> usize {
        let cp = clonal_population.lock()
            .expect("A thread paniced while holding the clonal population's lock.");
        cp.associated_inputs()
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

    /// Returns the copy of a random gene in a random [`ClonalPopulation`].
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock or the population is
    /// completely extinct.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    fn get_random_gene(&self) -> Gene {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .random_gene()
            // Unwrapping is safe here since we cannot call this function on empty populations.
            .unwrap()
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
            .expect("Another thread panicked while holding the clonal population lock.")
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
    fn append_population(&self, clonal_populations: Vec<ClonalPopulation>) {
        self.population.lock()
            .expect("A thread paniced while holding the population lock.")
            .append(clonal_populations);
    }

}
