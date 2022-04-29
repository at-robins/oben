//! The `population` module contains the administrative part of the evolutionary network.
extern crate rand;
extern crate rand_distr;
extern crate rmp_serde;
extern crate serde;
extern crate uuid;

use super::chemistry::{Information, Input, Output, Reaction, State};
use super::environment::{Environment, MutationCompendium};
use super::gene::{CrossOver, Gene, Genome};
use super::helper::{ActionChain, Iteration};
use super::protein::{InputSensor, OutputSensor, Receptor, Substrate};
use super::resource::Resource;
use rand::{thread_rng, Rng};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{Read, Write};
use std::marker::PhantomData;
use std::path::Path;
use std::rc::Rc;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};
use uuid::Uuid;

pub struct Organism<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
    OutputElementType,
    OutputSensorType,
> {
    substrates: Vec<Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
    input: InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>,
    output:
        OutputSensor<ReactionType, StateType, InformationType, OutputElementType, OutputSensorType>,
    time_alive: Iteration,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    Organism<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Creates a new `Organism` from the specified [`Substrate`]s.
    ///
    /// # Parameters
    ///
    /// * `substrates` - all [`Substrate`]s that make up the organism
    /// * `input` - the [`Substrate`]s linked to sensorical input
    /// * `output` - the [`Substrate`]s linked to the output
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn new(
        substrates: Vec<Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        input: InputSensor<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
        output: OutputSensor<
            ReactionType,
            StateType,
            InformationType,
            OutputElementType,
            OutputSensorType,
        >,
    ) -> Self {
        Organism {
            substrates,
            input,
            output,
            time_alive: Iteration::new(),
        }
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
    pub fn live(&mut self, environment: &Environment) -> Duration {
        let birth = Instant::now();
        // Set the current time to the last update, so the organsim can be reused (set_input()).
        let mut actions: ActionChain<Rc<Receptor<ReactionType, StateType, InformationType>>> =
            self.time_alive().into();
        let mut finished = false;
        // Add all receptors detecting changes to the input.
        for receptor in self.input.cascading_receptors() {
            actions.push_action(receptor);
        }
        // Run all receptors and subsequently add receptors detecting substrates,
        // which were modified during the run.
        // If the task takes longer than the specified threshold,
        // the run will be aborted.
        while !actions.is_empty()
            && !finished
            && birth.elapsed() < environment.lifespan()
            && self.binary_size() < environment.max_organism_size()
        {
            let time = actions.current_iteration();
            let mut input_feedback_changes: HashMap<usize, InformationType> = HashMap::new();
            let mut output_feedback_changes: HashMap<usize, InformationType> = HashMap::new();
            actions
                .pop_actions()
                .iter()
                .filter(|receptor| receptor.detect(time))
                .flat_map(|receptor| {
                    let result = receptor.catalyse(time);
                    // Fill the input feedback association map.
                    for (associations, value) in result.input_feedback_associations.into_iter() {
                        for association in associations {
                            input_feedback_changes.insert(association, value.clone());
                        }
                    }
                    // Fill the output feedback association map.
                    for (associations, value) in result.output_feedback_associations.into_iter() {
                        for association in associations {
                            output_feedback_changes.insert(association, value.clone());
                        }
                    }
                    // Set the finished flag if necessary.
                    finished = finished || self.output.is_finished(result.output_finish_substrate);
                    // Return the receptors to be processed in the subsequent step.
                    result.cascading_receptors
                })
                .for_each(|cascading_receptor| {
                    actions.push_action(cascading_receptor);
                });
            // Update input feedback substrates and add all receptors if changes to the input were detected.
            for receptor in self.input.feedback_update(input_feedback_changes) {
                actions.push_action(receptor);
            }
            // Update output feedback substrates and add all receptors if changes to the output were detected.
            for receptor in self.output.feedback_update(output_feedback_changes, time) {
                actions.push_action(receptor);
            }
        }
        self.time_alive = actions.current_iteration();
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
        self.substrates
            .iter()
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
    pub fn get_result(&mut self) -> OutputElementType {
        self.output.output_as_element(self.time_alive())
    }

    /// Returns the time the `Organism` is alive.
    pub fn time_alive(&self) -> Iteration {
        self.time_alive
    }

    /// Sets the input of the `Organism` to the specified values.
    ///
    /// # Prameters
    ///
    /// * `input` - the input values to set
    pub fn set_input(&mut self, input: InputElementType) {
        self.input.set_input(input);
    }
}

#[derive(Debug, PartialEq, Clone)]
/// Information of about an [`Organism`]s performance.
///
/// [`Organism`]: ./struct.Organism.html
pub struct OrganismInformation<SupplierResultInformationType, ExecutionResultElementType> {
    result: ExecutionResultElementType,
    result_info: SupplierResultInformationType,
    genome_size: usize,
    run_time: Duration,
    max_run_time: Duration,
    associated_inputs: usize,
    associated_outputs: usize,
    organism_size: usize,
    max_organism_size: usize,
}

impl<SupplierResultInformationType, ExecutionResultElementType>
    OrganismInformation<SupplierResultInformationType, ExecutionResultElementType>
{
    /// Creates a new `OrganismInformation` containig information about the performance of an
    /// [`Organism`] during a specific task.
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn new(
        result: ExecutionResultElementType,
        result_info: SupplierResultInformationType,
        genome_size: usize,
        run_time: Duration,
        max_run_time: Duration,
        associated_inputs: usize,
        associated_outputs: usize,
        organism_size: usize,
        max_organism_size: usize,
    ) -> Self {
        OrganismInformation {
            result,
            result_info,
            genome_size,
            run_time,
            max_run_time,
            associated_inputs,
            associated_outputs,
            organism_size,
            max_organism_size,
        }
    }

    /// Returns the result of one testing run as calculated by the [`Organism`].
    ///
    /// [`Organism`]: ./struct.Organism.html
    pub fn result(&self) -> &ExecutionResultElementType {
        &self.result
    }

    /// Returns the result related information about the testing run supplied by the
    /// supplier function; e.g. the real result.
    pub fn result_info(&self) -> &SupplierResultInformationType {
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
/// A `Individual` is describes an [`Organism`] with a defined [`Genome`].
///
/// [`Genome`]: ../gene/struct.Genome.html
/// [`Organism`]: ./struct.Organism.html
pub struct Individual<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
    OutputElementType,
    OutputSensorType,
> {
    phantom_r: PhantomData<ReactionType>,
    phantom_s: PhantomData<StateType>,
    phantom_t: PhantomData<InformationType>,
    phantom_input_element: PhantomData<InputElementType>,
    phantom_input_sensor: PhantomData<InputSensorType>,
    phantom_output_element: PhantomData<OutputElementType>,
    phantom_output_sensor: PhantomData<OutputSensorType>,
    source: Uuid,
    genome: Genome<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >,
    fitness: Option<f64>,
    bytes: usize,
    age: u32,
    tested: u32,
    resources: f64,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    Individual<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Creates a new `Individual`.
    ///
    /// # Parameters
    ///
    /// * `uuid` - the UUID of the `Individual`
    /// * `genome` - the [`Genome`] of the `Individual`
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn new(
        uuid: Uuid,
        genome: Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    ) -> Self {
        let bytes = genome.binary_size();
        Individual {
            phantom_r: PhantomData,
            phantom_s: PhantomData,
            phantom_t: PhantomData,
            phantom_input_element: PhantomData,
            phantom_input_sensor: PhantomData,
            phantom_output_element: PhantomData,
            phantom_output_sensor: PhantomData,
            source: uuid,
            genome,
            fitness: None,
            bytes,
            age: 0,
            tested: 0,
            resources: 0.0,
        }
    }

    /// Returns the UUID of this `Individual`.
    pub fn uuid(&self) -> &Uuid {
        &self.source
    }

    /// Returns the number of bytes this `Individual` contains.
    pub fn bytes(&self) -> usize {
        self.bytes
    }

    /// Returns the age of this `Individual`.
    pub fn age(&self) -> u32 {
        self.age
    }

    /// Returns how often this `Individual` was already tested.
    pub fn times_tested(&self) -> u32 {
        self.tested
    }

    /// Returns the [`Resource`]s this `Individual` has accumulated.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn resources(&self) -> f64 {
        self.resources
    }

    /// Returns the fitness of this `Individual` if any.
    pub fn fitness(&self) -> Option<f64> {
        self.fitness
    }

    /// Checks whether the fitness value of this `Individual` was already set.
    pub fn has_fitness(&self) -> bool {
        self.fitness.is_some()
    }

    /// Returns the number of associated inputs for this `Individual` contains.
    pub fn associated_inputs(&self) -> usize {
        self.genome().input().number_of_associated_inputs()
    }

    /// Returns the number of associated outputs for this `Individual` contains.
    pub fn associated_outputs(&self) -> usize {
        self.genome().output().number_of_associated_outputs()
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
            let f_old = old_fitness * (self.tested as f64);
            self.fitness = Some((fitness + f_old) / ((self.tested + 1) as f64));
        } else {
            self.fitness = Some(fitness);
        }
        self.tested += 1;
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
    pub fn mate(
        &self,
        partner: Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    ) -> Genome<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    > {
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
    pub fn mate_and_mutate(
        &self,
        partner: Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
        mutations: &MutationCompendium<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
        environment: &Environment,
    ) -> Individual<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    > {
        let offspring_genome = self.mate(partner);
        if let Some(mutated_offspring_genome) = mutations.mutate(&offspring_genome) {
            // If the mutation was successful, return the mutated individual.
            Individual::new(environment.generate_uuid(), mutated_offspring_genome)
        } else {
            // If the mutation was not successful, return the recombined individual without any
            // mutations.
            Individual::new(environment.generate_uuid(), offspring_genome)
        }
    }

    /// Returns the [`Genome`] of this [`Individual`].
    ///
    /// # Panics
    ///
    /// If loading the [`Genome`] from its associated file failed.
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn genome(
        &self,
    ) -> &Genome<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    > {
        &self.genome
    }

    /// Spends the maximum amount of [`Resource`]s possible to generate offspring and returns the
    /// number of offspring generated this way.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn spend_resources_for_mating(&mut self) -> usize {
        let number_of_offspring = self.resources().floor() as usize;
        self.resources -= number_of_offspring as f64;
        number_of_offspring
    }

    /// Adds the specified amount of [`Resource`]s to the `Individual`.
    ///
    /// # Parameters
    ///
    /// * `amount` - the amount of [`Resource`]s to add
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    fn aquire_resources(&mut self, amount: f64) {
        self.resources += amount;
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A `SerialisablePopulation` is a simplified version of a [`Population`] that can easily
/// be serialised and deserialised.
///
/// [`Population`]: ./struct.Population.html
struct SerialisablePopulation<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
    OutputElementType,
    OutputSensorType,
> {
    individuals: Vec<
        Individual<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    >,
    resources: Resource,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    SerialisablePopulation<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Returns the UUID and fitness of the fittest individual of the population.
    fn fittest_individual(
        &self,
    ) -> (
        Option<Uuid>,
        Option<f64>,
        Option<
            Individual<
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
    ) {
        let (mut uuid, mut fitness, mut ind) = (None, None, None);
        let filter: Vec<
            &Individual<
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        > = self.individuals.iter().filter(|i| i.age() >= 1).collect();
        for individual in filter {
            match fitness {
                Some(fit)
                    if (individual.fitness.is_none() || fit >= individual.fitness.unwrap()) => {},
                _ => {
                    uuid = Some(*individual.uuid());
                    fitness = individual.fitness;
                    ind = Some(individual.clone());
                },
            }
        }
        (uuid, fitness, ind)
    }
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    From<
        &Population<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    >
    for SerialisablePopulation<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    fn from(
        pop: &Population<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    ) -> Self {
        let mut serialisable_population = SerialisablePopulation {
            individuals: vec![],
            resources: pop.resources,
        };
        for individual in pop.individuals.values() {
            match individual.lock() {
                Ok(ind) => serialisable_population.individuals.push((*ind).clone()),
                Err(err) => {
                    // A individual cannot end up in an invalid state, so serialising it
                    // does no harm.
                    let poisoned: Individual<
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                        OutputElementType,
                        OutputSensorType,
                    > = err.into_inner().clone();
                    serialisable_population.individuals.push(poisoned);
                },
            }
        }
        serialisable_population
    }
}

#[derive(Debug, Clone)]
/// A `Population` is a population of individuals with different [`Genome`]s.
///
/// [`Genome`]: ../gene/struct.Genome.html
pub struct Population<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
    OutputElementType,
    OutputSensorType,
> {
    individuals: HashMap<
        Uuid,
        Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    >,
    resources: Resource,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    Population<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    /// Creates a new `Population` from the specified individuals.
    ///
    /// # Parameters
    ///
    /// `founding_individuals` - the [`Individual`]s to group into a `Population`
    /// `resources` - the [`Resource`]s the `Population` has access to
    ///
    /// [`Individual`]: ./struct.Individual.html
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn new(
        founding_individuals: Vec<
            Individual<
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
        resources: Resource,
    ) -> Self {
        let mut individuals = HashMap::new();
        for ind in founding_individuals.into_iter() {
            individuals.insert(*ind.uuid(), Arc::new(Mutex::new(ind)));
        }
        Population {
            individuals,
            resources,
        }
    }

    /// Write this `Population` to a JSON file if possible.
    /// An error will be returned if writing to the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file the `Population` should be written to
    pub fn snapshot_to_file<P>(&self, path_to_file: P) -> Result<(), Box<dyn Error + 'static>>
    where
        P: AsRef<Path>,
    {
        let mut file = File::create(&path_to_file)?;
        let serialisable_population: SerialisablePopulation<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        > = self.into();
        // TODO: Remove debug print statement and return an PopulationInformation struct instead.
        let (id, fitness, ind) = serialisable_population.fittest_individual();
        println!(
            "\nFittest:\nPopulation: {:?}\nID: {:?}\nFitness: {:?}\n\n{:#?}\n\n",
            path_to_file.as_ref(),
            id,
            fitness,
            ind
        );
        let ser = rmp_serde::to_vec(&serialisable_population)?;
        Ok(file.write_all(&ser)?)
    }

    /// Load a `Population` from a JSON file if possible.
    /// An error will be returned if parsing the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file from which the `Population` should be loaded
    pub fn load_from_file<P>(path_to_file: P) -> Result<Self, Box<dyn Error>>
    where
        P: AsRef<Path>,
    {
        let mut file = File::open(&path_to_file)?;
        let mut file_content = Vec::new();
        file.read_to_end(&mut file_content)?;
        let serialisable_population: SerialisablePopulation<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        > = rmp_serde::from_read_ref(&file_content)?;
        Ok(serialisable_population.into())
    }

    /// Returns the [`Individual`]s that are part of this `Population`.
    ///
    /// [`Individual`]: ./struct.Individual.html
    pub fn individuals(
        &self,
    ) -> Vec<
        Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    > {
        self.individuals.values().map(|val| val.clone()).collect()
    }

    /// Inserts all the [`Individual`]s into the `Population`.
    ///
    /// # Parameters
    ///
    /// * `individuals` - the [`Individual`]s to append
    ///
    /// [`Individual`]: ./struct.Individual.html
    pub fn append(
        &mut self,
        individuals: Vec<
            Individual<
                ReactionType,
                StateType,
                InformationType,
                InputElementType,
                InputSensorType,
                OutputElementType,
                OutputSensorType,
            >,
        >,
    ) {
        for ind in individuals.into_iter() {
            if let Some(ele) = self
                .individuals
                .insert(*ind.uuid(), Arc::new(Mutex::new(ind)))
            {
                // Should - for some reason - a duplicate UUID arise, repatriate the resources of
                // the removed individual.
                let removed_resources = ele
                    .lock()
                    .expect("A thread paniced while holding the individual's lock.")
                    .resources()
                    + 1.0;
                self.repatriate_resources(removed_resources);
            }
        }
    }

    /// Returns the copy of a random gene in a random [`Individual`] if there are any.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ./struct.Individual.html
    pub fn random_gene(&self) -> Option<Gene<ReactionType, StateType, InformationType>> {
        if self.individuals.len() == 0 {
            None
        } else {
            let random_population_index = thread_rng().gen_range(0..self.individuals.len());
            for (index, value) in self.individuals.values().enumerate() {
                if random_population_index == index {
                    return Some(
                        value
                            .lock()
                            .expect("A thread paniced while holding the individual's lock.")
                            .genome()
                            .duplicate_random_gene(),
                    );
                }
            }
            // This None is unreachable as the value must have been set before.
            None
        }
    }

    /// Returns a random [`Genome`] of a random [`Individual`] if there are any.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ./struct.Individual.html
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn random_genome(
        &self,
    ) -> Option<
        Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    > {
        if self.individuals.len() == 0 {
            None
        } else {
            let random_population_index = thread_rng().gen_range(0..self.individuals.len());
            self.individuals
                .values()
                .nth(random_population_index)
                .and_then(|value| {
                    Some(
                        value
                            .lock()
                            .expect("A thread paniced while holding the individual's lock.")
                            .genome()
                            .duplicate(),
                    )
                })
        }
    }

    /// Returns a random [`Genome`] of a random [`Individual`] if there are any.
    /// The chance of selecting an individual is based on its respective fitness.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ./struct.Individual.html
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn random_genome_fitness_based(
        &self,
    ) -> Option<
        Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    > {
        if self.individuals.len() == 0 {
            None
        } else {
            let individuals: Vec<
                Arc<
                    Mutex<
                        Individual<
                            ReactionType,
                            StateType,
                            InformationType,
                            InputElementType,
                            InputSensorType,
                            OutputElementType,
                            OutputSensorType,
                        >,
                    >,
                >,
            > = self
                .individuals
                .values()
                .map(|arc| Arc::clone(arc))
                .collect();
            let fitness_values: Vec<f64> = individuals
                .iter()
                .map(|ind| {
                    ind.lock()
                        .expect("A thread paniced while holding the individual's lock.")
                        .fitness()
                        .unwrap_or(0.0f64)
                })
                .collect();
            let random_fitness_sum: f64 = thread_rng().gen_range(0.0..=fitness_values.iter().sum());
            let mut random_population_index: usize = 0;
            let mut fitness_sum: f64 = 0.0;
            for fitness in fitness_values {
                fitness_sum += fitness;
                if random_fitness_sum <= fitness_sum {
                    break;
                } else {
                    random_population_index += 1;
                }
            }

            individuals.get(random_population_index).and_then(|value| {
                Some(
                    value
                        .lock()
                        .expect("A thread paniced while holding the individual's lock.")
                        .genome()
                        .duplicate(),
                )
            })
        }
    }

    /// Returns a random [`Individual`] if there is any.
    ///
    /// # Panics
    ///
    /// If another thread paniced while holding the individual's lock.
    ///
    /// [`Individual`]: ./struct.Individual.html
    pub fn random_individual(
        &self,
    ) -> Option<
        Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    > {
        if self.individuals.len() == 0 {
            None
        } else {
            let random_population_index = thread_rng().gen_range(0..self.individuals.len());
            self.individuals
                .values()
                .nth(random_population_index)
                .and_then(|value| Some(value.clone()))
        }
    }

    /// Returns the fittest [`Individual`] of the population that has a minimum age as specified.
    /// If the population is empty or no [`Individual`] meets the age criterium, `None` is returned.
    /// 
    /// # Parameters
    /// 
    /// * `minimum_age` - the minimum age to take the [`Individual`] into account
    pub fn fittest_individual(
        &self,
        minimum_age: u32,
    ) -> Option<
        Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
    > {
        let individuals_of_age: Vec<
            Arc<
                Mutex<
                    Individual<
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                        OutputElementType,
                        OutputSensorType,
                    >,
                >,
            >,
        > = self
            .individuals
            .values()
            .filter(|ind| {
                ind.lock()
                    .expect("A thread paniced while holding the individual's lock.")
                    .age()
                    >= minimum_age
            })
            .map(|arc| Arc::clone(arc))
            .collect();
        let mut fittest_individual = None;
        let mut max_fitness: f64 = 0.0;
        for individual in individuals_of_age.into_iter() {
            let current_fitness: f64 = individual
                .lock()
                .expect("A thread paniced while holding the individual's lock.")
                .fitness()
                .unwrap_or(0.0);
            if fittest_individual.is_none() || current_fitness > max_fitness {
                max_fitness = current_fitness;
                fittest_individual = Some(individual);
            }
        }
        fittest_individual
    }

    /// Remove the [`Individual`] from the `Population`.
    /// An error will be returned if there is no [`Individual`] with the specified UUID or
    /// if the correspondig [`Genome`] file could not be moved to the extinct sub-folder.
    ///
    /// # Parameters
    ///
    /// * `individual_uuid` - the UUID of the [`Individual`] to remove
    /// * `environment` - the [`Environment`] the population is growing in
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    /// [`Environment`]: ../environment/struct.Environment.html
    /// [`Individual`]: ./struct.Individual.html
    pub fn remove(
        &mut self,
        individual_uuid: Uuid,
    ) -> Result<
        Arc<
            Mutex<
                Individual<
                    ReactionType,
                    StateType,
                    InformationType,
                    InputElementType,
                    InputSensorType,
                    OutputElementType,
                    OutputSensorType,
                >,
            >,
        >,
        Box<dyn Error>,
    > {
        let removed = self
            .individuals
            .remove(&individual_uuid)
            .ok_or::<RemoveError>(RemoveError::new(&individual_uuid))?;
        Ok(removed)
    }

    /// Resets the fitness and age of all [`Individual`] as if they were never tested.
    ///
    /// [`Individual`]: ./struct.Individual.html
    pub fn reset_fitness(&self) {
        for individual in self.individuals.values() {
            let mut ind = individual
                .lock()
                .expect("A thread paniced while holding the individual's lock.");
            ind.fitness = None;
            ind.age = 0;
        }
    }

    /// Increments the age of every [`Individual`] by 1 generation.
    ///
    /// [`Individual`]: ./struct.Individual.html
    pub fn increment_age(&self) {
        for individual in self.individuals.values() {
            individual
                .lock()
                .expect("A thread paniced while holding the individual's lock.")
                .age += 1;
        }
    }

    /// Returns the number of [`Individual`]s that are part of this `Population`.
    pub fn size(&self) -> usize {
        self.individuals.len()
    }

    /// Calculates the mean fitness of the [`Individual`]s that are part of this `Population`.
    pub fn mean_fitness(&self) -> f64 {
        let mut count = 0.0;
        if self.size() > 0 {
            let mut fitness = 0.0;
            for individual in self.individuals.values() {
                let ind = individual
                    .lock()
                    .expect("A thread paniced while holding the individual's lock.");
                if let Some(f) = ind.fitness() {
                    fitness += f;
                    count += 1.0;
                }
            }
            fitness / count
        } else {
            // If the population is empty, the fitness is zero.
            0.0
        }
    }

    /// Calculates the mean genome size in bytes of the [`Individual`]s
    /// that are part of this `Population`.
    pub fn mean_genome_size(&self) -> f64 {
        if self.size() > 0 {
            let mut g_size = 0.0;
            for individual in self.individuals.values() {
                let ind = individual
                    .lock()
                    .expect("A thread paniced while holding the individual's lock.");
                g_size += ind.bytes() as f64;
            }
            g_size / (self.size() as f64)
        } else {
            // If the population is empty, the genome size is zero.
            0.0
        }
    }

    /// Add the specified amount of [`Resource`]s.
    ///
    /// # Parameters
    ///
    /// * `amount` - the amount of [`Resource`]s to add
    ///
    /// # Panics
    ///
    /// If the specified `amount` is not a valid positive number.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn repatriate_resources(&mut self, amount: f64) {
        self.resources.repatriate_resources(amount);
    }

    /// Recycles inavailable [`Resource`]s at the end of a generation.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn recycle(&mut self) {
        self.resources.recycle();
    }

    /// Returns the [`Resource`]s of the `Population`.
    ///
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn resources(&self) -> Resource {
        self.resources
    }

    /// Distributes available [`Resource`]s based on the fitness of the [`Individual`]s.
    ///
    /// [`Individual`]: ./struct.Individual.html
    /// [`Resource`]: ../resource/struct.Resource.html
    pub fn distribute_resources(&mut self) {
        let mut requests: Vec<(
            Arc<
                Mutex<
                    Individual<
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                        OutputElementType,
                        OutputSensorType,
                    >,
                >,
            >,
            f64,
        )> = Vec::with_capacity(self.individuals.capacity());
        let mut total_request = 0.0;
        // Calculate the maximum resources per individual that might be aquired.
        for individual in self.individuals() {
            if let Some(fitness) = individual
                .lock()
                .expect("Another thread panicked while holding the individual lock.")
                .fitness()
            {
                // Aquire resources. The higher the fitness, the slighter the difference needed
                // for significant resource advantage.
                let request = if fitness < 0.5 {
                    1.0 + 2.0 * fitness
                } else {
                    196.0 * fitness - 96.0
                };
                total_request += request;
                requests.push((individual.clone(), request));
            }
        }
        // Distribute resources based on available resources and claims.
        if total_request > 0.0 {
            let aquired_resources = self.resources.claim_resources(total_request);
            for (individual, request) in requests {
                let share = aquired_resources * (request / total_request);
                individual
                    .lock()
                    .expect("Another thread panicked while holding the individual lock.")
                    .aquire_resources(share);
            }
        }
    }
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        OutputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + Serialize
            + DeserializeOwned,
        OutputSensorType: Output<OutputElementType, InformationType>,
    >
    From<
        SerialisablePopulation<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    >
    for Population<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
        OutputElementType,
        OutputSensorType,
    >
{
    fn from(
        serial: SerialisablePopulation<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
            OutputElementType,
            OutputSensorType,
        >,
    ) -> Self {
        Self::new(serial.individuals, serial.resources)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
/// A `RemoveError` is returned when a UUID with no matching [`Individual`] is flagged for
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
        RemoveError {
            description: format!("No individual with UUID {} is present in the population.", uuid),
        }
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
