//! The `population` module contains the administrative part of the evolutionary network.
extern crate bitvec;
extern crate rand;
extern crate rand_distr;
extern crate rmp_serde;
extern crate uuid;
extern crate serde;

use std::collections::VecDeque;
use std::time::{Duration, Instant};
use std::rc::Rc;
use std::cell::RefCell;
use bitvec::{boxed::BitBox, order::Local};
use super::protein::{Substrate, Receptor};
use super::gene::{Genome, GenomeMutation};
use super::environment::Environment;
use uuid::Uuid;
use serde::{Serialize, Deserialize};
use rand_distr::{Binomial, Distribution};
use std::sync::{Arc, Mutex};
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::collections::HashMap;
use std::io::{Read, Write};

pub struct Organism {
    substrates: Vec<Rc<RefCell<Substrate>>>,
    input: Vec<Option<Rc<RefCell<Substrate>>>>,
    output: Vec<Option<Rc<RefCell<Substrate>>>>,
}

impl Organism {
    /// Creates a new `Organism` from the specified [`Substrate`]s.
    ///
    /// # Parameters
    ///
    /// * `substrates` - all [`Substrate`]s that make up the organism
    /// * `input` - the [`Substrate`]s linked to sensorical input
    /// * `output` - the [`Substrate`]s linked to the output
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn new(substrates: Vec<Rc<RefCell<Substrate>>>,
        input: Vec<Option<Rc<RefCell<Substrate>>>>,
        output: Vec<Option<Rc<RefCell<Substrate>>>>) -> Self {
        Organism{substrates, input, output}
    }

    /// Starts activity of all [`Receptor`]s and [`CatalyticCentre`]s linked to the [`Substrate`]s
    /// of this `Organism`. Execution will be aborted if the execution takes longer than the
    /// lifespan defined by the [`Environment`]. Returns the execution time.
    ///
    /// # Parameters
    ///
    /// * `environment` - the [`Environment`] the `Organism` lives in
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`CatalyticCentre`]: ../protein/struct.CatalyticCentre.html
    /// [`Receptor`]: ../protein/struct.Receptor.html
    /// [`Environment`]: ../environment/struct.Environment.html
    pub fn live(&self, environment: &Environment) -> Duration {
        let birth = Instant::now();
        let mut actions = VecDeque::<Rc<Receptor>>::new();
        // Add all receptors detecting changes to the input.
        for input_substrate in self.input.iter().filter_map(|val| val.as_ref()) {
            for receptor in input_substrate.borrow().receptors() {
                actions.push_back(receptor);
            }
        }
        // Run all receptors and subsequently add receptors detecting substrates,
        // which were modified during the run.
        // If the task takes longer than the specified threshold,
        // the run will be aborted.
        while !actions.is_empty() && birth.elapsed() < environment.lifespan() {
            for cascading_receptor in actions.pop_front().unwrap().detect() {
                actions.push_back(cascading_receptor);
            }
        }
        // Return the time it took to perform the task, but never more than the lifespan.
        let run_time = birth.elapsed();
        if run_time < environment.lifespan() {
            run_time
        } else {
            environment.lifespan()
        }
    }

    /// Returns the number of [`Substrate`]s this `Organism` consists of.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn number_of_substrates(&self) -> usize {
        self.substrates.len()
    }

    /// Returns the values of output [`Substrate`]s if any.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_result(&self) -> Vec<Option<BitBox<Local, u8>>> {
        self.output.iter().map(|sub| sub.as_ref().and_then(|some| Some(some.borrow().value().clone()))).collect()
    }

    /// Sets the input of the `Organism` to the specified values.
    ///
    /// # Prameters
    ///
    /// * `input` - the input values to set
    ///
    /// # Panics
    ///
    /// If the length of `input` is not equal to the number of inputs specified by the organism.
    pub fn set_input(&self, input: Vec<BitBox<Local, u8>>) {
        if input.len() != self.input.len() {
            panic!("The organism needs {} inputs, but {} were specified.", self.input.len(), input.len());
        }
        for ip in input.into_iter().enumerate() {
            if let Some(substrate) = self.input[ip.0].as_ref() {
                substrate.borrow_mut().set_value(ip.1);
            }
        }
    }
}

/// Generates the specified number of randomly mutated versions of the specified [`Genome`].
///
/// # Parameters
///
/// * `number_of_mutated_genomes` - the number of mutated [`Genome`]s to produce
/// * `genome` - the source [`Genome`] to mutate
///
/// [`Genome`]: ../gene/struct.Genome.html
fn produce_mutated_genomes(number_of_mutated_genomes: u32, genome: &Genome) -> Vec<Genome>{
    (0..number_of_mutated_genomes).map(|_| rand::random::<GenomeMutation>())
        .filter_map(|mutation| mutation.mutate(genome))
        .collect()
}

#[derive(Debug, PartialEq, Clone)]
/// Information of about an [`Organism`]s performance.
///
/// [`Organism`]: ./struct.Organism.html
pub struct OrganismInformation {
    bytes: usize,
    run_time: Duration,
    max_run_time: Duration,
}

impl OrganismInformation {
    /// Creates a new `OrganismInformation` containig information about the performance of an
    /// [`Organism`] during a specific task.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn new(bytes: usize, run_time: Duration, max_run_time: Duration) -> Self {
        OrganismInformation{bytes, run_time, max_run_time}
    }

    /// The size of the [`Organism`]s [`Genome`] in bytes.
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`Organism`]: ./struct.Organism.html
    pub fn bytes(&self) -> usize {
        self.bytes
    }

    /// The time it took the [`Organism`] to perform a task.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn run_time(&self) -> Duration {
        self.run_time
    }

    /// The maximal time a [`Organism`] is allowed to spend before timing out.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn max_runtime(&self) -> Duration {
        self.max_run_time
    }
}


#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A `ClonalPopulation` is a population of individuals with exactly the same [`Genome`].
///
/// [`Genome`]: ../gene/struct.Genome.html
pub struct ClonalPopulation {
    source: Uuid,
    genome: Genome,
    size: f64,
    fitness: Option<f64>,
    bytes: usize,
}

impl ClonalPopulation {
    /// Founds a new `ClonalPopulation` with a single individual.
    ///
    /// # Parameters
    ///
    /// * `uuid` - the UUID of the `ClonalPopulation`
    /// * `genome` - the [`Genome`] of the population
    /// * `environment` - the [`Environment`] the population is living in
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`Environment`]: ../environment/struct.Environment.html
    pub fn found(uuid: Uuid, genome: Genome, environment: &Environment) -> Self {
        let bytes = genome.binary_size();
        ClonalPopulation{
            source: uuid,
            genome,
            size: 1.0 / environment.population_size() as f64,
            fitness: None,
            bytes,
        }
    }

    /// Returns the UUID of this `ClonalPopulation`.
    pub fn uuid(&self) -> &Uuid {
        &self.source
    }

    /// Returns the number of bytes this `ClonalPopulation` contains.
    pub fn bytes(&self) -> usize {
        self.bytes
    }

    /// Returns the relative size of this `ClonalPopulation`.
    pub fn relative_size(&self) -> f64 {
        self.size
    }

    /// Checks whether the fitness value of this `ClonalPopulation` was already set.
    pub fn has_fitness(&self) -> bool {
        self.fitness.is_some()
    }

    /// Adjusts the relative size of this `ClonalPopulation` to the size of the whole reference
    /// [`Population`].
    ///
    /// # Parameters
    ///
    /// * `reference` - the reference population's size
    ///
    /// [`Population`]: ./struct.Population.html
    pub fn adjust_size_to_reference(&mut self, reference: f64) {
        self.fitness = self.fitness.map(|f| f / reference);
    }

    /// Adds the specified fitness value by setting the fitness to the mean of old and
    /// and new fitness. If no fitness was set prevously, sets the specified fitness
    /// as new fitness.
    ///
    /// # Parameters
    ///
    /// * `fitness` - the new fitness to add to the current fitness value
    fn add_fitness(&mut self, fitness: f64) {
        if let Some(old_fitness) = self.fitness {
            self.fitness = Some((fitness + old_fitness) / 2.0);
        } else {
            self.fitness = Some(fitness);
        }
    }

    /// Generates mutated [`Genome`]s and updates the population size based on the evaluated fitness.
    ///
    /// # Parameters
    ///
    /// * `fitness` - the newly evaluated fitness of the population
    /// * `environment` - the [`Environment`] the population is growing in
    ///
    /// # Panics
    ///
    /// If internal loading the source [`Genome`] from its associated file failed,
    /// while generating mutants.
    ///
    /// [`Environment`]: ../environment/struct.Environment.html
    pub fn evaluate_new_fitness(&mut self, fitness: f64, environment: &Environment) -> Vec<Genome> {
        self.add_fitness(fitness);
        // The fitness was just set, so the unwrap call must succeed.
        let relative_total_offspring = self.size * self.fitness.unwrap();
        // Determine the number of genetically different offspring.
        let absolute_total_offspring = relative_total_offspring.floor() as u64;
        let absolute_mutated_offspring = Binomial::new(absolute_total_offspring, environment.mutation_rate()).expect("The mutation rate is not set between 0 and 1.").sample(&mut rand::thread_rng());
        // Grow the population by the non-mutated offspring.
        self.size += relative_total_offspring - (absolute_mutated_offspring as f64 / environment.population_size() as f64);
        // Generate the mutations.
        produce_mutated_genomes(absolute_mutated_offspring as u32, &self.genome)
    }

    /// Returns the [`Genome`] of this [`ClonalPopulation`].
    ///
    /// # Panics
    ///
    /// If loading the [`Genome`] from its associated file failed.
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn genome(&self) -> &Genome {
        &self.genome
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A `SerialisablePopulation` is a simplified version of a [`Population`] that can easily
/// be serialised and deserialised.
///
/// [`Population`]: ./struct.Population.html
struct SerialisablePopulation {
    clonal_populations: Vec<ClonalPopulation>,
}

impl SerialisablePopulation {
    /// Returns the UUID and fitness of the fittest individual of the population.
    fn fittest_individual(&self) -> (Option<Uuid>, Option<f64>) {
        let (mut uuid, mut fitness) = (None, None);
        for individual in &self.clonal_populations {
            match fitness {
                Some(fit) if individual.fitness.is_some() && fit >= individual.fitness.unwrap() => {},
                _ => { uuid = Some(*individual.uuid());
                    fitness = individual.fitness;
                }
            }
        }
        (uuid, fitness)
    }
}

impl From<&Population> for SerialisablePopulation {
    fn from(pop: &Population) -> Self {
        let mut serialisable_population = SerialisablePopulation{clonal_populations: vec!()};
        for clonal_population in pop.clonal_populations.values() {
            match clonal_population.lock() {
                Ok(cp) => serialisable_population.clonal_populations.push((*cp).clone()),
                Err(err) => {
                    // A clonal population cannot end up in an invalid state, so serialising it
                    // does no harm.
                    let poisoned: ClonalPopulation = err.into_inner().clone();
                    serialisable_population.clonal_populations.push(poisoned);
                },
            }
        }
        serialisable_population
    }
}

#[derive(Debug, Clone)]
/// A `Population` is a population of clonal sub-populations with different [`Genome`]s.
///
/// [`Genome`]: ../gene/struct.Genome.html
pub struct Population {
    clonal_populations: HashMap<Uuid, Arc<Mutex<ClonalPopulation>>>,
}

impl Population {
    /// Creates a new `Population` from the specified clonal sub-populations.
    ///
    /// # Parameters
    ///
    /// `sub_populations` - the [`ClonalPopulation`]s to group into a `Population`
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn new(sub_populations: Vec<ClonalPopulation>) -> Self {
        let mut clonal_populations = HashMap::new();
        for cp in sub_populations.into_iter() {
            clonal_populations.insert(*cp.uuid(), Arc::new(Mutex::new(cp)));
        }
        Population{clonal_populations}
    }

    /// Write this `Population` to a JSON file if possible.
    /// An error will be returned if writing to the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file the `Population` should be written to
    pub fn snapshot_to_file<P>(&self, path_to_file: P) -> Result<(), Box<dyn Error + 'static>> where P: AsRef<Path> {
        let mut file = File::create(&path_to_file)?;
        let serialisable_population: SerialisablePopulation = self.into();
        // TODO: Remove debug print statement and return an PopulationInformation struct instead.
        let (id, fitness) = serialisable_population.fittest_individual();
        println!("Population: {:?}\nID: {:?}\nFitness: {:?}", path_to_file.as_ref(), id, fitness);
        let ser = rmp_serde::to_vec(&serialisable_population)?;
        file.write_all(&ser)?;
        Ok(file.sync_all()?)
    }

    /// Load a `Population` from a JSON file if possible.
    /// An error will be returned if parsing the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file from which the `Population` should be loaded
    pub fn load_from_file<P>(path_to_file: P) -> Result<Self, Box<dyn Error>> where P: AsRef<Path> {
        let mut file = File::open(&path_to_file)?;
        let mut file_content = Vec::new();
        file.read_to_end(&mut file_content)?;
        let serialisable_population: SerialisablePopulation = rmp_serde::from_read_ref(&file_content)?;
        Ok(serialisable_population.into())
    }

    /// Returns the [`ClonalPopulation`]s that are part of this `Population`.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn clonal_populations(&self) -> Vec<Arc<Mutex<ClonalPopulation>>> {
        self.clonal_populations.values().map(|val| val.clone()).collect()
    }

    /// Inserts all the [`ClonalPopulation`]s into the `Population` and returns a reference
    /// to them.
    ///
    /// # Parameters
    ///
    /// * `clonal_populations` - the [`ClonalPopulation`]s to append
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn append(&mut self, clonal_populations: Vec<ClonalPopulation>) -> Vec<Arc<Mutex<ClonalPopulation>>> {
        clonal_populations.into_iter()
            .map(|clonal_population| (clonal_population.uuid().clone(), Arc::new(Mutex::new(clonal_population))))
            .map(|clonal_population| {
                self.clonal_populations.insert(clonal_population.0, clonal_population.1.clone());
                clonal_population.1
            }).collect()
    }

    /// Remove the [`ClonalPopulation`] from the `Population`.
    /// An error will be returned if there is no [`ClonalPopulation`] with the specified UUID or
    /// if the correspondig [`Genome`] file could not be moved to the extinct sub-folder.
    ///
    /// # Parameters
    ///
    /// * `clonal_population_uuid` - the UUID of the [`ClonalPopulation`] to remove
    /// * `environment` - the [`Environment`] the population is growing in
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`Environment`]: ../environment/struct.Environment.html
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn remove(&mut self, clonal_population_uuid: Uuid, environment: &Environment) -> Result<Arc<Mutex<ClonalPopulation>>, Box<dyn Error>> {
        let removed = self.clonal_populations.remove(&clonal_population_uuid)
            .ok_or::<RemoveError>(RemoveError::new(&clonal_population_uuid))?;
        // Move the extinct genome to a backup folder.
        std::fs::rename(environment.genome_path(&clonal_population_uuid), environment.extinct_genome_path(&clonal_population_uuid))?;
        Ok(removed)
    }
}

impl From<SerialisablePopulation> for Population {
    fn from(serial: SerialisablePopulation) -> Self {
        Self::new(serial.clonal_populations)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
/// A `RemoveError` is returned when a UUID with no matching [`ClonalPopulation`] is flagged for
/// removal.
pub struct RemoveError {
    description: String,
}

impl RemoveError {
    /// Creates a `RemoveError` from the specified UUID.
    ///
    /// # Parameters
    ///
    /// * `uuid` - the UUID flagged for removal
    pub fn new(uuid: &Uuid) -> Self {
        RemoveError {description: format!("No clonal population with UUID {} is present in the population.", uuid)}
    }
}

impl std::fmt::Display for RemoveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description)
    }
}

impl Error for RemoveError {
    fn description(&self) -> &str {
        &self.description
    }
}
