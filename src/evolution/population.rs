//! The `population` module contains the administrative part of the evolutionary network.
extern crate rand;
extern crate rand_distr;
extern crate rmp_serde;
extern crate uuid;
extern crate serde;

use std::collections::VecDeque;
use std::time::{Duration, Instant};
use std::rc::Rc;
use std::cell::RefCell;
use super::binary::BinarySubstrate;
use super::protein::{Substrate, Receptor};
use super::gene::{CrossOver, Gene, Genome, GenomeMutation};
use super::environment::Environment;
use uuid::Uuid;
use serde::{Serialize, Deserialize};
// use rand_distr::{Binomial, Distribution};
use std::sync::{Arc, Mutex};
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::collections::HashMap;
use std::io::{Read, Write};
use rand::{thread_rng, Rng};

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
        while !actions.is_empty() && birth.elapsed() < environment.lifespan() && self.binary_size() < environment.max_organism_size() {
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

    /// Returns the number of bits of all [`Substrate`]s that are part of the `Organism`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn binary_size(&self) -> usize {
        self.substrates.iter()
            .map(|s| s.borrow().binary_size())
            .sum()
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
    pub fn get_result(&self) -> Vec<Option<BinarySubstrate>> {
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
    pub fn set_input(&self, input: Vec<BinarySubstrate>) {
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

#[derive(Debug, PartialEq, Clone)]
/// Information of about an [`Organism`]s performance.
///
/// [`Organism`]: ./struct.Organism.html
pub struct OrganismInformation<I> {
    result: Vec<Option<BinarySubstrate>>,
    result_info: I,
    genome_size: usize,
    run_time: Duration,
    max_run_time: Duration,
    associated_inputs: usize,
    associated_outputs: usize,
    organism_size: usize,
    max_organism_size: usize,
}

impl<I> OrganismInformation<I> {
    /// Creates a new `OrganismInformation` containig information about the performance of an
    /// [`Organism`] during a specific task.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn new(result: Vec<Option<BinarySubstrate>>, result_info: I, genome_size: usize, run_time: Duration, max_run_time: Duration, associated_inputs: usize, associated_outputs: usize, organism_size: usize, max_organism_size: usize) -> Self {
        OrganismInformation{
            result,
            result_info,
            genome_size,
            run_time,
            max_run_time,
            associated_inputs,
            associated_outputs,
            organism_size,
            max_organism_size
        }
    }

    /// Returns the result of one testing run as calculated by the [`Organism`].
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn result(&self) -> &Vec<Option<BinarySubstrate>> {
        &self.result
    }

    /// Returns the result related information about the testing run supplied by the
    /// supplier function; e.g. the real result.
    pub fn result_info(&self) -> &I {
        &self.result_info
    }

    /// The size of the [`Organism`]s [`Genome`] in bit.
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`Organism`]: ./struct.Organism.html
    pub fn genome_size(&self) -> usize {
        self.genome_size
    }

    /// The size of the [`Organism`]s in bit.
    /// This currently only takes [`Substrate`] values into account.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Organism`]: ./struct.Organism.html
    pub fn organism_size(&self) -> usize {
        self.organism_size
    }

    /// The maximum size of an [`Organism`]s in bit.
    /// This currently only takes [`Substrate`] values into account.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Organism`]: ./struct.Organism.html
    pub fn max_organism_size(&self) -> usize {
        self.max_organism_size
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

    /// Returns the execution time relative to the maximal run time allowed for execution
    /// of a task.
    pub fn relative_run_time(&self) -> f64 {
        self.run_time().as_secs_f64() / self.max_runtime().as_secs_f64()
    }

    /// Returns the size of an [`Organism`] relative to the maximal size.
    /// This can be more than 1 since the maximum size only is a cap for further task execution.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn relative_organism_size(&self) -> f64 {
        self.organism_size() as f64 / self.max_organism_size() as f64
    }

    /// The number of inputs of the [`Organism`] that are already associated.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn associated_inputs(&self) -> usize {
        self.associated_inputs
    }

    /// The number of outputs of the [`Organism`] that are already associated.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn associated_outputs(&self) -> usize {
        self.associated_outputs
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A `ClonalPopulation` is a population of individuals with exactly the same [`Genome`].
///
/// [`Genome`]: ../gene/struct.Genome.html
pub struct ClonalPopulation {
    source: Uuid,
    genome: Genome,
    fitness: Option<f64>,
    bytes: usize,
    age: u32
}

impl ClonalPopulation {
    /// Founds a new `ClonalPopulation` with a single individual.
    ///
    /// # Parameters
    ///
    /// * `uuid` - the UUID of the `ClonalPopulation`
    /// * `genome` - the [`Genome`] of the population
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn found(uuid: Uuid, genome: Genome) -> Self {
        let bytes = genome.binary_size();
        ClonalPopulation{
            source: uuid,
            genome,
            fitness: None,
            bytes,
            age: 0,
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

    /// Returns the age of this `ClonalPopulation`.
    pub fn age(&self) -> u32 {
        self.age
    }

    /// Returns the fitness of this `ClonalPopulation` if any.
    pub fn fitness(&self) -> Option<f64> {
        self.fitness
    }

    /// Checks whether the fitness value of this `ClonalPopulation` was already set.
    pub fn has_fitness(&self) -> bool {
        self.fitness.is_some()
    }

    /// Returns the number of associated inputs for this `ClonalPopulation` contains.
    pub fn associated_inputs(&self) -> usize {
        self.genome().number_of_associated_inputs()
    }

    /// Returns the number of associated outputs for this `ClonalPopulation` contains.
    pub fn associated_outputs(&self) -> usize {
        self.genome().number_of_associated_outputs()
    }

    /// Returns the values of all outputs for this `ClonalPopulation`.
    pub fn output_values(&self) -> Vec<Option<BinarySubstrate>> {
        self.genome().get_output_values()
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
            // Calculate the mean of the current and all previous fitness values.
            let f_old = old_fitness * (self.age as f64);
            self.fitness = Some((fitness  + f_old) / ((self.age + 1) as f64));
        } else {
            self.fitness = Some(fitness);
        }
        self.age += 1;
    }

    /// Updates the `Population`'s fitness.
    ///
    /// # Parameters
    ///
    /// * `fitness` - the newly evaluated fitness of the population
    pub fn evaluate_new_fitness(&mut self, fitness: f64) {
        self.add_fitness(fitness);
    }

    /// Recombine the [`Genome`] of this `Individual` and its mating partner and return the
    /// resulting [`Genome`].
    ///
    /// # Parameters
    ///
    /// * `partner` - the mating partner's [`Genome`]
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn mate(&self, partner: Genome) -> Genome {
        self.genome().cross_over(&partner)
    }

    /// Recombine the [`Genome`] of this `Individual` and its mating partner and return the
    /// resulting `Individual`.
    ///
    /// # Parameters
    ///
    /// * `partner` - the mating partner's [`Genome`]
    /// * `environment` - the [`Environment`] the `Individual` is living in
    ///
    /// [`Environment`]: ../environment/struct.Environment.html
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn mate_and_mutate(&self, partner: Genome, environment: &Environment) -> ClonalPopulation {
        let offspring_genome = self.mate(partner);
        let random_chance: f64 = rand::thread_rng().gen_range(0.0, 1.0);
        // Calculate the number of mutations corresponding to the generated uniform random percentage.
        //    P("n mutations in a single genome") = "mutation rate" ^ n
        // => n = log(base: "mutation rate", value: P)
        let number_of_mutations = random_chance.log(environment.mutation_rate()).floor() as usize;
        if let Some(mutated_offspring_genome) = GenomeMutation::mutate_n_times(number_of_mutations, &offspring_genome) {
            // If the mutation was successful, return the mutated individual.
            ClonalPopulation::found(environment.generate_uuid(), mutated_offspring_genome)
        } else {
            // If the mutation was not successful, return the recombined individual without any
            // mutations.
            ClonalPopulation::found(environment.generate_uuid(), offspring_genome)
        }
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
    fn fittest_individual(&self) -> (Option<Uuid>, Option<f64>, Option<ClonalPopulation>, Option<Vec<Option<BinarySubstrate>>>) {
        let (mut uuid, mut fitness, mut cp, mut out) = (None, None, None, None);
        let filter: Vec<&ClonalPopulation> = self.clonal_populations.iter().filter(|i| i.age() >= 1).collect();
        for individual in filter {
            match fitness {
                Some(fit) if (individual.fitness.is_none() || fit >= individual.fitness.unwrap()) => {},
                _ => { uuid = Some(*individual.uuid());
                    fitness = individual.fitness;
                    cp = Some(individual.clone());
                    out = Some(individual.output_values())
                }
            }
        }
        (uuid, fitness, cp, out)
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
        let (id, fitness, cp, out) = serialisable_population.fittest_individual();
        println!("\nFittest:\nPopulation: {:?}\nID: {:?}\nFitness: {:?}\n\n{:#?}\nValues: {:?}\n\n",
            path_to_file.as_ref(), id, fitness, cp, out);
        let ser = rmp_serde::to_vec(&serialisable_population)?;
        Ok(file.write_all(&ser)?)
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

    /// Inserts all the [`ClonalPopulation`]s into the `Population`.
    ///
    /// # Parameters
    ///
    /// * `clonal_populations` - the [`ClonalPopulation`]s to append
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn append(&mut self, clonal_populations: Vec<ClonalPopulation>) {
        for cp in clonal_populations.into_iter() {
            self.clonal_populations.insert(*cp.uuid(), Arc::new(Mutex::new(cp)));
        }
    }

    /// Returns the copy of a random gene in a random [`ClonalPopulation`] if there are any.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn random_gene(&self) -> Option<Gene> {
        if self.clonal_populations.len() == 0 {
            None
        } else {
            let random_population_index = thread_rng().gen_range(0, self.clonal_populations.len());
            for (index, value) in self.clonal_populations.values().enumerate() {
                if random_population_index == index {
                    return Some(value.lock()
                        .expect("A thread paniced while holding the clonal population's lock.")
                        .genome()
                        .duplicate_random_gene())
                }
            }
            // This None is unreachable as the value must have been set before.
            None
        }
    }

    /// Returns a random [`Genome`] of a random [`ClonalPopulation`] if there are any.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn random_genome(&self) -> Option<Genome> {
        if self.clonal_populations.len() == 0 {
            None
        } else {
            let random_population_index = thread_rng().gen_range(0, self.clonal_populations.len());
            self.clonal_populations.values()
                .nth(random_population_index)
                .and_then(|value| Some(value.lock()
                    .expect("A thread paniced while holding the clonal population's lock.")
                    .genome()
                    .duplicate())
                )
        }
    }

    /// Returns a random [`ClonalPopulation`] if there is any.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the clonal population's lock.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn random_individual(&self) -> Option<Arc<Mutex<ClonalPopulation>>> {
        if self.clonal_populations.len() == 0 {
            None
        } else {
            let random_population_index = thread_rng().gen_range(0, self.clonal_populations.len());
            self.clonal_populations.values()
                .nth(random_population_index)
                .and_then(|value| Some(value.clone()))
        }
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
    pub fn remove(&mut self, clonal_population_uuid: Uuid) -> Result<Arc<Mutex<ClonalPopulation>>, Box<dyn Error>> {
        let removed = self.clonal_populations.remove(&clonal_population_uuid)
            .ok_or::<RemoveError>(RemoveError::new(&clonal_population_uuid))?;
        Ok(removed)
    }

    /// Resets the fitness and age of all [`ClonalPopulation`] as if they were never tested.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn reset_fitness(&self) {
        for clonal_population in self.clonal_populations.values() {
            let mut cp = clonal_population.lock()
                .expect("A thread paniced while holding the clonal population's lock.");
            cp.fitness = None;
            cp.age = 0;
        }
    }

    /// Returns the number of [`ClonalPopulation`]s that are part of this `Population`.
    pub fn size(&self) -> usize {
        self.clonal_populations.len()
    }

    /// Calculates the mean fitness of the [`ClonalPopulation`]s that are part of this `Population`.
    pub fn mean_fitness(&self) -> f64 {
        if self.size() > 0 {
            let mut fitness = 0.0;
            for clonal_population in self.clonal_populations.values() {
                let cp = clonal_population.lock()
                    .expect("A thread paniced while holding the clonal population's lock.");
                if let Some(f) = cp.fitness() {
                    fitness += f;
                }
            }
            fitness / (self.size() as f64)
        } else {
            // If the population is empty, the fitness is zero.
            0.0
        }
    }

    /// Calculates the mean genome size in bytes of the [`ClonalPopulation`]s
    /// that are part of this `Population`.
    pub fn mean_genome_size(&self) -> f64 {
        if self.size() > 0 {
            let mut g_size = 0.0;
            for clonal_population in self.clonal_populations.values() {
                let cp = clonal_population.lock()
                    .expect("A thread paniced while holding the clonal population's lock.");
                g_size += cp.bytes() as f64;
            }
            g_size / (self.size() as f64)
        } else {
            // If the population is empty, the genome size is zero.
            0.0
        }
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
