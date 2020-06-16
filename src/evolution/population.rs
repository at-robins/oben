//! The `population` module contains the administrative part of the evolutionary network.
extern crate bitvec;
extern crate rand;
extern crate rand_distr;
extern crate serde_json;
extern crate uuid;
extern crate serde;

use std::collections::VecDeque;
use std::time::Instant;
use std::rc::Rc;
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

pub struct Organism<'a> {
    substrates: Vec<Substrate>,
    input: Vec<&'a Substrate>,
    output: Vec<&'a Substrate>,
}

impl<'a> Organism<'a> {
    /*pub fn new(input_substrates: u32, output_substrates: u32) -> Self {

    }*/
    pub fn live(&self, environment: &Environment) {
        let birth = Instant::now();
        let mut actions = VecDeque::<Rc<Receptor>>::new();
        // Add all receptors detecting changes to the input.
        for input_substrate in &self.input {
            for receptor in input_substrate.receptors() {
                actions.push_back(receptor);
            }
        }
        // Run all receptors and subsequently add receptors detecting substrates,
        // which were modified during the run.
        // If the task takes longer than the specified threshold,
        // the run will be aborted.
        while !actions.is_empty() && birth.elapsed() <= environment.lifespan() {
            for cascading_receptor in actions.pop_front().unwrap().detect() {
                actions.push_back(cascading_receptor);
            }
        }
    }

    pub fn get_result(&self) -> Vec<BitBox<Local, u8>> {
        self.output.iter().map(|sub| sub.value().clone()).collect()
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

#[derive(Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub struct ClonalPopulation {
    source: Uuid,
    size: u32,
    fitness: Option<f64>,
    death_counter: f64,
}

impl ClonalPopulation {
    /// Founds a new `ClonalPopulation` with a single individual.
    ///
    /// # Parameters
    ///
    /// * `uuid` - the UUID of the `ClonalPopulation`
    pub fn found(uuid: Uuid) -> Self {
        ClonalPopulation{
            source: uuid,
            size: 1,
            fitness: None,
            death_counter: 0.0,
        }
    }

    /// Returns the UUID of this `ClonalPopulation`.
    pub fn uuid(&self) -> &Uuid {
        &self.source
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
            self.fitness = Some((fitness + old_fitness) / 2f64);
        } else {
            self.fitness = Some(fitness);
        }
    }

    /// Sets adds the specified value to the size of the `ClonalPopulation` if it does not
    /// exceed the internally specified maximum population size.
    ///
    /// # Parameters
    ///
    /// * `size` - the population size by which to grow the `ClonalPopulation`
    /// * `environment` - the [`Environment`] the population is growing in
    ///
    /// [`Environment`]: ../environment/struct.Environment.html
    fn grow(&mut self, size: u32, environment: &Environment) {
        let new_size = self.size + size;
        if new_size > environment.max_clonal_population_size() {
            self.size = environment.max_clonal_population_size();
        } else {
            self.size = new_size;
        }
    }

    /// Kills the specified amount of individuals. The remainder is stored and taken into account
    /// on subsequent death events.
    ///
    /// # Parameters
    ///
    /// * `death_count` - the number of individuals that died
    ///
    /// # Panics
    ///
    /// If `death_count` is negative.
    pub fn death_event(&mut self, death_count: f64) {
        if death_count < 0.0 {
            panic!("Death count cannot be negative: {}", death_count);
        } else {
            self.death_counter += death_count;
            let actual_deaths = self.death_counter.floor();
            self.size -= actual_deaths as u32;
            self.death_counter -= actual_deaths;
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
        let total_offspring = (self.size as f64 * self.fitness.unwrap()) as u64;
        // Determine the number of genetically different offspring.
        let mutated_offspring = Binomial::new(total_offspring, environment.mutation_rate()).expect("The mutation rate is not set between 0 and 1.").sample(&mut rand::thread_rng());
        // Grow the population by the non-mutated offspring.
        self.grow((total_offspring - mutated_offspring) as u32, environment);
        // Generate the mutations.
        let genome = self.get_genome(environment);
        produce_mutated_genomes(mutated_offspring as u32, &genome)
    }

    /// Loads the [`Genome`] from its associated file.
    ///
    /// # Parameter
    ///
    /// * `environment` - the [`Environment`] the population is growing in
    ///
    /// # Panics
    ///
    /// If loading the [`Genome`] from its associated file failed.
    ///
    /// [`Environment`]: ../environment/struct.Environment.html
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn get_genome(&self, environment: &Environment) -> Genome {
        match Genome::load_from_file(environment.genome_path(&self.source)) {
            Ok(genome) => genome,
            // I/O is a subsistantial part of the system, so its advised to panic here instead
            // of performing error handling.
            Err(err) => panic!("Loading of genome {} resulted in error: {}", self.source, err),
        }
    }

    /// Checks if this `ClonalPopulation` is already extinct.
    pub fn is_extinct(&self) -> bool {
        self.size == 0
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct SerialisablePopulation {
    clonal_populations: Vec<ClonalPopulation>,
}

#[derive(Debug, Clone)]
pub struct Population {
    clonal_populations: HashMap<Uuid, Arc<Mutex<ClonalPopulation>>>,
}

impl Population {
    /// Create an immutable, serialisable snapshot of the current population.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn snapshot(&self) -> SerialisablePopulation {
        let mut serialisable_population = SerialisablePopulation{clonal_populations: vec!()};
        for clonal_population in self.clonal_populations.values() {
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

    /// Write this `Population` to a JSON file if possible.
    /// An error will be returned if writing to the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file the `Population` should be written to
    pub fn snapshot_to_file<P>(&self, path_to_file: P) -> Result<(), Box<dyn Error + 'static>> where P: AsRef<Path> {
        let file = File::create(path_to_file)?;
        let serialisable_population = self.snapshot();
        let write_success = serde_json::to_writer(file, &serialisable_population)?;
        Ok(write_success)
    }

    /// Returns the [`ClonalPopulation`]s that are part of this `Population`.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    pub fn clonal_populations(&self) -> Vec<&Arc<Mutex<ClonalPopulation>>> {
        self.clonal_populations.values().collect()
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
    pub fn remove(&mut self, clonal_population_uuid: &Uuid, environment: &Environment) -> Result<Arc<Mutex<ClonalPopulation>>, Box<dyn Error>> {
        let removed = self.clonal_populations.remove(&clonal_population_uuid)
            .ok_or::<RemoveError>(RemoveError::new(clonal_population_uuid))?;
        // Move the extinct genome to a backup folder.
        {
            // We do not care for the integrity of the clonal population since the UUID is immuatble.
            let clonal_population = removed.lock().map_or_else(|val| val.into_inner(), |val| val);
            std::fs::rename(environment.genome_path(clonal_population.uuid()), environment.extinct_genome_path(clonal_population.uuid()))?;
        }
        Ok(removed)
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
