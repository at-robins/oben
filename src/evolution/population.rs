//! The `population` module contains the administrative part of the evolutionary network.
extern crate bitvec;
extern crate rand;
extern crate rand_distr;
extern crate serde_json;
extern crate uuid;
extern crate serde;

use std::collections::VecDeque;
use std::time::{Duration, Instant};
use std::rc::Rc;
use bitvec::{boxed::BitBox, order::Local};
use super::protein::{Substrate, Receptor};
use super::gene::{Genome, GenomeMutation};
use super::environment::Environment;
use uuid::Uuid;
use serde::{Serialize, Deserialize};
use rand_distr::{Binomial, Distribution};

pub struct Individual<'a> {
    substrates: Vec<Substrate>,
    input: Vec<&'a Substrate>,
    output: Vec<&'a Substrate>,
    lifespan: Duration,
}

impl<'a> Individual<'a> {
    /*pub fn new(input_substrates: u32, output_substrates: u32) -> Self {

    }*/
    pub fn live(&self) {
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
        while !actions.is_empty() && birth.elapsed() <= self.lifespan {
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
        match Genome::load_from_file(environment.genome_path(self.source)) {
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

pub struct Population {

}
