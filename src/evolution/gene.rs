//! The `gene` module contains information structures to describe,
//! replicate and mutate the evolutionary network.
extern crate rand;
extern crate rmp_serde;
extern crate serde;

use super::chemistry::{Information, Input, Reaction, State};
use super::helper::{a_or_b, do_a_or_b};
use super::population::Organism;
use super::protein::{CatalyticCentre, InputSensor, Receptor, Substrate};
use rand::{thread_rng, Rng};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{Read, Write};
use std::marker::PhantomData;
use std::num::NonZeroUsize;
use std::path::Path;
use std::rc::{Rc, Weak};

/// A `Genome` is a collection of individual [`Gene`]s and associations between them.
/// A `Genome` is required to consist of 1 or more genes.
///
/// [`Gene`]: ./struct.Gene.html
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType> {
    phantom_reaction: PhantomData<ReactionType>,
    phantom_state: PhantomData<StateType>,
    phantom_information: PhantomData<InformationType>,
    phantom_input_element: PhantomData<InputElementType>,
    phantom_input_sensor: PhantomData<InputSensorType>,
    input: GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    >,
    output: Vec<Option<GeneSubstrate>>,
    genes: Vec<Gene<ReactionType, StateType, InformationType>>,
    associations: Vec<GeneAssociation<InformationType>>,
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
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
    > Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates a new `Genome` containing the specified inputs, outputs and [`Gene`]s.
    /// A `Genome` needs at least 1 [`Gene`].
    ///
    /// # Parameters
    ///
    /// * `input` - the input substrates
    /// * `output` - the output substrates
    /// * `genes` - the initial [`Gene`]s
    ///
    /// # Panics
    ///
    /// If the vector of `genes` is empty.
    pub fn new(
        input: GenomicInputSensor<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
        output: Vec<Option<GeneSubstrate>>,
        genes: Vec<Gene<ReactionType, StateType, InformationType>>,
    ) -> Self {
        if genes.is_empty() {
            panic!("A genome needs at least one gene.");
        }
        Genome {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_input_element: PhantomData,
            phantom_input_sensor: PhantomData,
            input,
            output,
            genes,
            associations: Vec::new(),
        }
    }

    /// Get the number of [`Gene`]s in this `Genome`.
    /// A `Genome` must encode 1 or more [`Gene`]s.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn number_of_genes(&self) -> NonZeroUsize {
        NonZeroUsize::new(self.genes.len())
            .expect("No gene is encoded by this genome. This is forbidden by the contract.")
    }

    /// Get the number of associated outputs [`Substrate`]s in this `Genome`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn number_of_associated_outputs(&self) -> usize {
        self.output.iter().filter(|i| i.is_some()).count()
    }

    /// Get the number of output [`Substrate`]s in this `Genome`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn number_of_outputs(&self) -> usize {
        self.output.len()
    }

    /// Get the number of [`GeneAssociation`]s in this `Genome`.
    ///
    /// [`GeneAssociation`]: ./struct.GeneAssociation.html
    pub fn number_of_associations(&self) -> usize {
        self.associations.len()
    }

    /// Returns the input.
    pub fn input(
        &self,
    ) -> &GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    > {
        &self.input
    }

    /// Returns the input as mutable.
    pub fn input_mut(
        &mut self,
    ) -> &mut GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    > {
        &mut self.input
    }

    /// Returns the output at the specified index if any.
    pub fn output(&self, output_index: usize) -> Option<Option<GeneSubstrate>> {
        self.output.get(output_index).and_then(|inner| Some(*inner))
    }

    /// Returns the association at the specified index if any.
    pub fn association(
        &self,
        association_index: usize,
    ) -> Option<&GeneAssociation<InformationType>> {
        self.associations.get(association_index)
    }

    /// Returns the association at the specified index as mutable if any.
    pub fn association_mut(
        &mut self,
        association_index: usize,
    ) -> Option<&mut GeneAssociation<InformationType>> {
        self.associations.get_mut(association_index)
    }

    /// Sets the output association at the specified index to the specified value
    /// and returns the previous value.
    ///
    /// # Parameters
    ///
    /// * `output_index` - the index of the output association to changes
    /// * `output_value` - the new value for the specified output association
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    pub fn set_output(
        &mut self,
        output_index: usize,
        output_value: Option<GeneSubstrate>,
    ) -> Option<GeneSubstrate> {
        if output_index >= self.output.len() {
            panic!("The output vector of this genome is of length {}, but insertion of element {:#?} at index {} was attempted.", self.number_of_outputs(), output_value, output_index);
        }
        let old_value = self.output[output_index];
        self.output[output_index] = output_value;
        old_value
    }

    /// Returns a copy of the output [`Substrate`]s.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_output_values(&self) -> Vec<Option<InformationType>> {
        self.output
            .iter()
            .map(|substrate| {
                substrate
                    .and_then(|a| self.get_substrate(a))
                    .and_then(|inner| Some(inner.clone()))
            })
            .collect()
    }

    /// Returns the specified [`Substrate`] if contained within the `Genome`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_substrate(&self, substrate: GeneSubstrate) -> Option<&InformationType> {
        self.genes
            .get(substrate.gene)
            .and_then(|gene| gene.substrates.get(substrate.substrate))
    }

    /// Returns the specified [`Substrate`] as mutable if contained within the `Genome`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_substrate_mut(&mut self, substrate: GeneSubstrate) -> Option<&mut InformationType> {
        self.genes
            .get_mut(substrate.gene)
            .and_then(|gene| gene.substrates.get_mut(substrate.substrate))
    }

    /// Adds a [`Gene`] to the `Genome` if possible and returns the index of the
    /// new gene.
    /// This function will fail if the underlying vector would overflow due to the addition.
    ///
    /// # Parameters
    ///
    /// * `gene` - the [`Gene`] to add to the `Genome`
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn add_gene(
        &mut self,
        gene: Gene<ReactionType, StateType, InformationType>,
    ) -> Option<usize> {
        if let Some(new_index) = self.genes.len().checked_add(1) {
            self.genes.push(gene);
            Some(new_index)
        } else {
            None
        }
    }

    /// Adds a substrate to a [`Gene`] if possible and returns a reference to the
    /// new substrate.
    ///
    /// # Parameters
    ///
    /// * `gene` - the [`Gene`] to add to the substrate to
    /// * `substrate` - the substrate to add
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn add_substrate_to_gene(
        &mut self,
        gene: usize,
        substrate: InformationType,
    ) -> Option<GeneSubstrate> {
        if self.number_of_genes().get() > gene {
            self.genes[gene]
                .add_substrate(substrate)
                .and_then(|sub_index| Some(GeneSubstrate::new(gene, sub_index)))
        } else {
            None
        }
    }

    /// Removes the [`Gene`] at the specified index and returns it.
    /// All input-, output- and substrate-associations with the specified gene
    /// are cleared.
    /// The number of [`Gene`]s in the `Genome` cannot be reduced to zero.
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene to remove
    ///
    /// # Panics
    ///
    /// If the index is out of bounds or the `Genome` contains only a single [`Gene`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn remove_gene(&mut self, gene: usize) -> Gene<ReactionType, StateType, InformationType> {
        if self.number_of_genes().get() <= 1 {
            panic!("A genome needs to contain at least one gene, so no gene can be removed.");
        }
        // Remove all inputs pointing to the removed gene.
        self.input.remove_associations_with_gene(gene);
        // Remove all outputs pointing to the removed gene.
        for output_value in &mut self.output {
            if output_value
                .and_then(|output| Some(output.is_gene(gene)))
                .unwrap_or(false)
            {
                *output_value = None;
            }
        }
        // Remove all substrate associations pointing to the removed gene.
        for association in &mut self.associations {
            association.remove_associated_gene(gene);
        }
        self.genes.remove(gene)
    }

    /// Removes the substrate of the specified [`Gene`] and returns it.
    /// All input-, output- and substrate-associations with the specified substrate
    /// are adjusted.
    /// The number of substrates of a [`Gene`] cannot be reduced to zero.
    ///
    /// # Parameters
    ///
    /// * `substrate` - the substrate to remove
    ///
    /// # Panics
    ///
    /// If the indices are out of bounds or the [`Gene`] contains only a single substrate.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn remove_substrate(&mut self, substrate: GeneSubstrate) -> InformationType {
        if self.get_gene(substrate.gene()).number_of_substrates().get() <= 1 {
            panic!(
                "A gene needs to contain at least one substrate, so no substrate can be removed."
            );
        }
        let removed = self.genes[substrate.gene()].remove_substrate(substrate.substrate());
        self.input_mut()
            .adjust_after_gene_substrate_removal(substrate);
        self.adjust_output_after_gene_substrate_removal(substrate);
        self.adjust_associations_after_gene_substrate_removal(substrate);
        self.validate_associations();
        removed
    }

    /// Duplicates the [`Gene`] at the specified index and returns it
    /// for further processing.
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene to duplicate
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn duplicate_gene(&self, gene: usize) -> Gene<ReactionType, StateType, InformationType> {
        self.genes[gene].duplicate()
    }

    /// Duplicates a random [`Gene`] and returns it
    /// for further processing.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn duplicate_random_gene(&self) -> Gene<ReactionType, StateType, InformationType> {
        self.get_gene(self.get_random_gene()).duplicate()
    }

    /// Returns the [`Gene`] at the specified index.
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene to duplicate
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn get_gene(&self, gene: usize) -> &Gene<ReactionType, StateType, InformationType> {
        &self.genes[gene]
    }

    /// Returns the [`Gene`] at the specified index as mutable.
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene to duplicate
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn get_gene_mut(
        &mut self,
        gene: usize,
    ) -> &mut Gene<ReactionType, StateType, InformationType> {
        &mut self.genes[gene]
    }

    /// Checks if the [`Substrate`] is present in the `Genome`.
    ///
    /// # Parameters
    ///
    /// * `genes` - the [`Gene`]'s of the `Genome`
    /// * `substrate` - the substrate to check for
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn has_substrate(
        genes: &Vec<Gene<ReactionType, StateType, InformationType>>,
        substrate: &GeneSubstrate,
    ) -> bool {
        if let Some(gene) = genes.get(substrate.gene) {
            gene.substrates.len() > substrate.substrate
        } else {
            false
        }
    }

    /// Checks if the [`Substrate`] is an output to the `Genome`.
    ///
    /// # Parameters
    /// * `substrate` - the substrate to check for
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn has_output_substrate(&self, substrate: &GeneSubstrate) -> bool {
        self.output
            .iter()
            .filter_map(|potential_output| potential_output.as_ref())
            .any(|output| output == substrate)
    }

    /// Returns all substrates of the specified [`Gene`] combined with their [`GeneSubstrate`] pointer.
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    fn translate_get_gene_substrates(
        &self,
        gene: usize,
    ) -> Vec<(GeneSubstrate, Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>)>
    {
        let gene_reference = &self.genes[gene];
        (0..gene_reference.number_of_substrates().get())
            .map(|substrate| GeneSubstrate::new(gene, substrate))
            .map(|gene_substrate| {
                (
                    gene_substrate,
                    Rc::new(RefCell::new(Substrate::new(
                        gene_reference.substrates[gene_substrate.substrate].clone(),
                    ))),
                )
            })
            .collect()
    }

    /// Returns the index of a random [`Gene`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn get_random_gene(&self) -> usize {
        thread_rng().gen_range(0, self.number_of_genes().get())
    }

    /// Returns the index of a random output [`GeneSubstrate`] if there is any.
    ///
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    pub fn get_random_output(&self) -> Option<usize> {
        if self.number_of_outputs() > 0 {
            Some(thread_rng().gen_range(0, self.number_of_outputs()))
        } else {
            None
        }
    }

    /// Returns the index of a random [`GeneAssociation`] if there is any.
    ///
    /// [`GeneSubstrate`]: ./struct.GeneAssociation.html
    pub fn get_random_association(&self) -> Option<usize> {
        if self.number_of_associations() > 0 {
            Some(thread_rng().gen_range(0, self.number_of_associations()))
        } else {
            None
        }
    }

    /// Duplicates the `Genome` and all its contents.
    pub fn duplicate(&self) -> Self {
        // At the moment this is just a wrapper for cloning.
        self.clone()
    }

    /// Generates a random [`GeneSubstrate`] pointing to a [`Substrate`] of a [`Gene`]
    /// contained within the `Genome`.
    ///
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Gene`]: ./struct.Gene.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn random_gene_substrate(&self) -> GeneSubstrate {
        let random_gene = self.get_random_gene();
        let random_substrate = self.get_gene(random_gene).get_random_substrate();
        GeneSubstrate::new(random_gene, random_substrate)
    }

    /// Adds a [`GeneAssociation`] to the `Genome` if possible and returns the index of the
    /// new gene.
    /// This function will fail if the underlying vector would overflow due to the addition.
    ///
    /// # Parameters
    ///
    /// * `association` - the [`GeneAssociation`] to add to the `Genome`
    ///
    /// [`GeneAssociation`]: ./struct.GeneAssociation.html
    pub fn add_association(
        &mut self,
        association: GeneAssociation<InformationType>,
    ) -> Option<usize> {
        if let Some(new_index) = self.associations.len().checked_add(1) {
            self.associations.push(association);
            Some(new_index)
        } else {
            None
        }
    }

    /// Removes the [`GeneAssociation`] at the specified index and returns it.
    ///
    /// # Parameters
    ///
    /// * `association` - the index of the association to remove
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`GeneAssociation`]: ./struct.GeneAssociation.html
    pub fn remove_association(&mut self, association: usize) -> GeneAssociation<InformationType> {
        self.associations.remove(association)
    }

    /// Adjust the output [`GeneSubstrate`] references after the binary [`Substrate`] of a contained
    /// [`Gene`] was removed.
    ///
    /// # Parameters
    ///
    /// * `removed_substrate` - an index based pointer to the removed [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_output_after_gene_substrate_removal(&mut self, removed_substrate: GeneSubstrate) {
        for output_ref in &mut self.output {
            if output_ref.is_some() {
                let mut current_reference = output_ref.unwrap();
                if current_reference == removed_substrate {
                    // Remove the output association if it points to the removed substrate.
                    *output_ref = None;
                } else if current_reference.is_gene(removed_substrate.gene) {
                    // Adjust the output association if it points to the gene from which the
                    // substrate was removed.
                    current_reference.substrate =
                        adjust_index(current_reference.substrate, removed_substrate.substrate);
                    *output_ref = Some(current_reference);
                }
                // Otherwise, leave the output untouched.
            }
        }
    }

    /// Adjust the [`GeneAssociation`] references after the binary [`Substrate`] of a contained
    /// [`Gene`] was removed.
    ///
    /// # Parameters
    ///
    /// * `removed_substrate` - an index based pointer to the removed [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneAssociation`]: ./struct.GeneAssociation.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_associations_after_gene_substrate_removal(
        &mut self,
        removed_substrate: GeneSubstrate,
    ) {
        for association in &mut self.associations {
            association.adjust_after_gene_substrate_removal(removed_substrate);
        }
    }

    /// Validates the [`GeneSubstrate`] references of input and output after a recombination
    /// event and removes invalid ones.
    ///
    /// # Parameters
    ///
    /// * `genes` - the [`Gene`]s of the [`Genome`]
    /// * `io_substrates` - the input or output [`GeneSubstrate`]s
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    fn validate_input_ouput_associations(
        genes: &Vec<Gene<ReactionType, StateType, InformationType>>,
        io_substrates: &mut Vec<Option<GeneSubstrate>>,
    ) {
        let mut duplicate_checker = HashSet::new();
        for io_substrate in io_substrates {
            if let Some(current_substrate) = io_substrate {
                if !duplicate_checker.insert(current_substrate.clone())
                    || !Self::has_substrate(genes, current_substrate)
                {
                    // Remove the io association if it points to a duplicate or an invalid substrate.
                    *io_substrate = None;
                }
                // Otherwise, leave it untouched.
            }
        }
    }

    /// Validates the [`GeneAssociation`] references
    /// after a recombination event and removes invalid ones.
    ///
    /// [`GeneAssociation`]: ./struct.GeneAssociation.html
    fn validate_gene_substrate_associations(&mut self) {
        let genes = &self.genes;
        for association in &mut self.associations {
            association
                .associations
                .retain(|a| Self::has_substrate(genes, a));
        }
    }

    /// Validates the [`GeneAssociation`] references and the input-output [`GeneSubstrate`]s
    /// after a recombination event and removes invalid ones.
    ///
    /// [`GeneAssociation`]: ./struct.GeneAssociation.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    fn validate_associations(&mut self) {
        self.input.validate(&self.genes);
        Self::validate_input_ouput_associations(&self.genes, &mut self.output);
        self.validate_gene_substrate_associations();
    }

    /// Load a `Genome` from a JSON file if possible.
    /// An error will be returned if parsing the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file from which the `Genome` should be loaded
    pub fn load_from_file<P>(path_to_file: P) -> Result<Self, Box<dyn Error>>
    where
        P: AsRef<Path>,
    {
        let mut file = File::open(&path_to_file)?;
        let mut file_content = Vec::new();
        file.read_to_end(&mut file_content)?;
        let genome = rmp_serde::from_read_ref(&file_content)?;
        Ok(genome)
    }

    /// Write a `Genome` to a JSON file if possible.
    /// An error will be returned if writing to the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file the `Genome` should be written to
    pub fn write_to_file<P>(&self, path_to_file: P) -> Result<(), Box<dyn Error>>
    where
        P: AsRef<Path>,
    {
        let mut file = File::create(path_to_file)?;
        let ser = rmp_serde::to_vec(&self)?;
        file.write_all(&ser)?;
        Ok(file.sync_all()?)
    }

    /// Returns the binary size of this `Genome` in byte.
    ///
    /// # Panics
    ///
    /// If the underlying serialisation fails.
    pub fn binary_size(&self) -> usize {
        rmp_serde::to_vec(&self)
            .expect("Serialisation of the genome failed.")
            .len()
    }

    pub fn translate(
        &self,
    ) -> Organism<ReactionType, StateType, InformationType, InputElementType, InputSensorType> {
        let mut gene_substrate_map: HashMap<
            GeneSubstrate,
            Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
        > = HashMap::new();
        // Insert all genome level substrates.
        for gene_association in &self.associations {
            let genome_level_substrate =
                Rc::new(RefCell::new(Substrate::new(gene_association.substrate.clone())));
            for gene_substrate in &gene_association.associations {
                gene_substrate_map
                    .entry(gene_substrate.clone())
                    .or_insert(genome_level_substrate.clone());
            }
        }
        // Insert all gene level substrates without overwriting genome level ones.
        for gene_index in 0..self.number_of_genes().get() {
            for (gene_substrate, substrate) in
                self.translate_get_gene_substrates(gene_index).into_iter()
            {
                gene_substrate_map
                    .entry(gene_substrate)
                    .or_insert(substrate);
            }
        }
        // Translate receptors and catalytic centres.
        for gene_index in 0..self.number_of_genes().get() {
            for receptor in self.get_gene(gene_index).receptors.iter() {
                receptor.translate(gene_index, &gene_substrate_map);
            }
        }
        let substrates = gene_substrate_map.values().map(|sub| sub.clone()).collect();
        let input = self.input.translate(&gene_substrate_map);
        let output = self
            .output
            .iter()
            .map(|substrate| {
                substrate.and_then(|gene_substrate| {
                    gene_substrate_map
                        .get(&gene_substrate)
                        .and_then(|inner| Some(inner.clone()))
                })
            })
            .collect();
        Organism::new(substrates, input, output)
    }
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information + Default,
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType> + Default,
    > Default
    for Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    fn default() -> Self {
        Genome {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_input_element: PhantomData,
            phantom_input_sensor: PhantomData,
            input: GenomicInputSensor::default(),
            output: Vec::default(),
            genes: vec![Gene::<ReactionType, StateType, InformationType>::default()],
            associations: Vec::default(),
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
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
    > CrossOver
    for Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    fn is_similar(&self, other: &Self) -> bool {
        self.number_of_outputs() == other.number_of_outputs()
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let input = self.input.cross_over(&other.input);
            let output = self.output.cross_over(&other.output);
            let genes = self.genes.cross_over(&other.genes);
            let associations = self.associations.cross_over(&other.associations);
            let mut recombined = Genome {
                phantom_reaction: PhantomData,
                phantom_state: PhantomData,
                phantom_information: PhantomData,
                phantom_input_element: PhantomData,
                phantom_input_sensor: PhantomData,
                input,
                output,
                genes,
                associations,
            };
            // Remove invalid gene-substrate-associations.
            recombined.validate_associations();
            recombined
        } else {
            // If the two genomes are not similar return a random one.
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy, Serialize, Deserialize)]
/// A `GeneSubstrate` is an index based pointer to a [`Substrate`] encoded in a [`Gene`]
/// contained within the respective [`Genome`].
///
/// [`Gene`]: ./struct.Gene.html
/// [`Genome`]: ./struct.Genome.html
/// [`Substrate`]: ../protein/struct.Substrate.html
pub struct GeneSubstrate {
    // index of gene inside the containing genome
    gene: usize,
    // index of substrate inside the respective gene
    substrate: usize,
}

impl GeneSubstrate {
    /// Creates a new 'GeneSubstrate' reference.
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the referenced [`Gene`]
    /// * `substrate` - the index of the referenced [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn new(gene: usize, substrate: usize) -> Self {
        GeneSubstrate { gene, substrate }
    }

    /// Returns the gene index the reference is pointing to.
    pub fn gene(&self) -> usize {
        self.gene
    }

    /// Returns the substrate index the reference is pointing to.
    pub fn substrate(&self) -> usize {
        self.substrate
    }

    /// Checks wether this `GeneSubstrate` references the [`Gene`] at the specified index.
    ///
    /// # Parameters
    ///
    /// * `gene_index` - the index of the gene to check for
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn is_gene(&self, gene_index: usize) -> bool {
        self.gene == gene_index
    }

    // This function was never used but might be usefull at some point.
    /// Checks wether this `GeneSubstrate` references the [`Substrate`] at the specified index.
    ///
    /// # Parameters
    ///
    /// * `substrate_index` - the index of the substrate to check for
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn is_substrate(&self, substrate_index: usize) -> bool {
        self.substrate == substrate_index
    }
}

impl CrossOver for GeneSubstrate {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        let gene = self.gene.cross_over(&other.gene);
        let substrate = self.substrate.cross_over(&other.substrate);
        GeneSubstrate::new(gene, substrate)
    }
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
/// A `GeneAssociation` is a [`Substrate`] defined on [`Genome`] level that replaces [`Gene`] specific
/// [`Substrate`]s upon translation. This process is used to interconnect different [`Gene`]s on a
/// [`Genome`] level.
///
/// [`Gene`]: ./struct.Gene.html
/// [`Genome`]: ./struct.Genome.html
/// [`Substrate`]: ../protein/struct.Substrate.html
pub struct GeneAssociation<T> {
    // substrate value defined in the genome and shared between genes
    substrate: T,
    // gene specific substrates pointing to the shared substrate
    associations: Vec<GeneSubstrate>,
}

impl<T> GeneAssociation<T>
where
    T: Information,
{
    /// Creates a new `GeneAssociation` from the specified [`Substrate`]
    /// with no initial associations.
    ///
    /// # Parameters
    ///
    /// * `substrate` - the common [`Substrate`] for this association
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn new(substrate: T) -> Self {
        GeneAssociation {
            substrate,
            associations: Vec::new(),
        }
    }

    /// Returns the substrate of this association.
    pub fn substrate(&self) -> &T {
        &self.substrate
    }

    /// Returns the substrate of this association.
    pub fn substrate_mut(&mut self) -> &mut T {
        &mut self.substrate
    }

    /// Returns the number of substrates this association references.
    pub fn number_of_associated_substrates(&self) -> usize {
        self.associations.len()
    }

    /// Remove all associations with the specified [`Gene`].
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene to remove associations to
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn remove_associated_gene(&mut self, gene: usize) {
        self.associations.retain(|a| a.gene != gene);
    }

    /// Returns the index of a random [`GeneSubstrate`] if there are any.
    ///
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    pub fn get_random_gene_substrate(&self) -> Option<usize> {
        if self.associations.len() > 0 {
            Some(thread_rng().gen_range(0, self.associations.len()))
        } else {
            None
        }
    }

    /// Adds a [`GeneSubstrate`] to the `GeneAssociation` if possible and returns the index of the
    /// new [`GeneSubstrate`].
    /// This function will fail if the underlying vector would overflow due to the addition.
    ///
    /// # Parameters
    ///
    /// * `association` - the [`GeneSubstrate`] to add to the `GeneAssociation`
    ///
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    pub fn add_association(&mut self, association: GeneSubstrate) -> Option<usize> {
        if let Some(new_index) = self.associations.len().checked_add(1) {
            self.associations.push(association);
            Some(new_index)
        } else {
            None
        }
    }

    /// Removes a [`GeneSubstrate`] from the `GeneAssociation` and returns it.
    ///
    /// # Parameters
    ///
    /// * `association_index` - the [`GeneSubstrate`]s index
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    pub fn remove_association(&mut self, association_index: usize) -> GeneSubstrate {
        self.associations.remove(association_index)
    }

    /// Adjust the [`GeneSubstrate`] references after the binary [`Substrate`] of a contained
    /// [`Gene`] was removed.
    ///
    /// # Parameters
    ///
    /// * `removed_substrate` - an index based pointer to the removed [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_after_gene_substrate_removal(&mut self, removed_substrate: GeneSubstrate) {
        self.associations.retain(|a| a != &removed_substrate);
        for association in &mut self.associations {
            if association.is_gene(removed_substrate.gene) {
                // Adjust the input association if it points to the gene from which the
                // substrate was removed.
                association.substrate =
                    adjust_index(association.substrate, removed_substrate.substrate);
            }
            // Otherwise, leave the input untouched.
        }
    }
}

impl<T> CrossOver for GeneAssociation<T>
where
    T: Information,
{
    fn is_similar(&self, other: &Self) -> bool {
        self.substrate.is_similar(&other.substrate)
            && self.associations.len() == other.associations.len()
    }

    fn cross_over(&self, other: &Self) -> Self {
        let substrate = self.substrate.cross_over(&other.substrate);
        let associations = self.associations.cross_over(&other.associations);
        GeneAssociation {
            substrate,
            associations,
        }
    }
}

/// A `Gene` is an immutable structure encoding a self-contained network, but without
/// explicite function. It can be transcribed into a functional protein network.
/// A `Gene` is required to encode at least 1 [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct Gene<R, S, T> {
    phantom_r: PhantomData<R>,
    phantom_s: PhantomData<S>,
    phantom_t: PhantomData<T>,
    substrates: Vec<T>,
    receptors: Vec<GenomicReceptor<R, S, T>>,
}

impl<R: Reaction<T>, S: State<T>, T: Information> Gene<R, S, T> {
    /// Creates a new `Gene` containing the specified substrates.
    /// A `Gene` needs at least 1 substrate.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the initial substrates the gene contains
    ///
    /// # Panics
    ///
    /// If the vector of `substrates` is empty.
    pub fn new(substrates: Vec<T>) -> Self {
        if substrates.is_empty() {
            panic!("A gene needs at least one substrate.");
        }
        Gene {
            phantom_r: PhantomData,
            phantom_s: PhantomData,
            phantom_t: PhantomData,
            substrates,
            receptors: Vec::new(),
        }
    }

    /// Gets the number of [`Substrate`]s encoded by this `Gene`.
    /// A `Gene` must encode for 1 or more [`Substrate`]s.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn number_of_substrates(&self) -> NonZeroUsize {
        NonZeroUsize::new(self.substrates.len())
            .expect("No substrate is encoded by this gene. This is forbidden by the contract.")
    }

    /// Gets the number of [`GenomicReceptor`]s encoded by this `Gene`.
    ///
    /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
    pub fn number_of_receptors(&self) -> usize {
        self.receptors.len()
    }

    /// Returns the index of a random [`Substrate`] encoded by this `Gene`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_random_substrate(&self) -> usize {
        thread_rng().gen_range(0, self.number_of_substrates().get())
    }

    /// Returns the index of a random [`GenomicReceptor`] encoded by this `Gene` if there is any.
    ///
    /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
    pub fn get_random_receptor(&self) -> Option<usize> {
        if self.number_of_receptors() > 0 {
            Some(thread_rng().gen_range(0, self.number_of_receptors()))
        } else {
            None
        }
    }

    /// Adds a binary [`Substrate`] to the `Gene` if possible and returns the index of the
    /// new [`Substrate`].
    /// This function will fail if the underlying vector would overflow due to the addition.
    ///
    /// # Parameters
    ///
    /// * `substrate` - the [`Substrate`] to add to the `Gene`
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn add_substrate(&mut self, substrate: T) -> Option<usize> {
        if let Some(new_index) = self.substrates.len().checked_add(1) {
            self.substrates.push(substrate);
            Some(new_index)
        } else {
            None
        }
    }

    /// Adds a [`GenomicReceptor`] to the `Gene` if possible and returns the index of the
    /// new [`GenomicReceptor`].
    /// This function will fail if the underlying vector would overflow due to the addition.
    ///
    /// # Parameters
    ///
    /// * `receptor` - the [`Substrate`] to add to the `Gene`
    ///
    /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
    pub fn add_receptor(&mut self, receptor: GenomicReceptor<R, S, T>) -> Option<usize> {
        if let Some(new_index) = self.receptors.len().checked_add(1) {
            self.receptors.push(receptor);
            Some(new_index)
        } else {
            None
        }
    }

    /// Removes the receptor at the specified index from the `Gene` if possible and
    /// returns the removed receptor.
    ///
    /// # Parameters
    ///
    /// * `receptor_index` - the receptor's index to remove from the `Gene`
    ///
    /// # Panics
    ///
    /// If the specified receptor's index is out of bounds.
    pub fn remove_receptor(&mut self, receptor_index: usize) -> GenomicReceptor<R, S, T> {
        self.receptors.remove(receptor_index)
    }

    /// Returns the receptor at the specified index if possible.
    ///
    /// # Parameters
    ///
    /// * `receptor_index` - the receptor's index
    pub fn receptor(&self, receptor_index: usize) -> Option<&GenomicReceptor<R, S, T>> {
        self.receptors.get(receptor_index)
    }

    /// Returns the receptor at the specified index if possible as mutable.
    ///
    /// # Parameters
    ///
    /// * `receptor_index` - the receptor's index
    pub fn receptor_mut(&mut self, receptor_index: usize) -> Option<&mut GenomicReceptor<R, S, T>> {
        self.receptors.get_mut(receptor_index)
    }

    /// Removes the binary [`Substrate`] at the specified index from the `Gene` if possible and
    /// returns the removed [`Substrate`].
    ///
    /// # Parameters
    ///
    /// * `substrate_index` - the [`Substrate`]'s index to remove from the `Gene`
    ///
    /// # Panics
    ///
    /// If the specified [`Substrate`]'s index is out of bounds or there is just a single
    /// [`Substrate`] present, as a `Gene` needs at least 1 [`Substrate`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn remove_substrate(&mut self, substrate_index: usize) -> T {
        if self.number_of_substrates().get() <= 1 {
            panic!(
                "A genome needs to contain at least one substrate, so no substrate can be removed."
            );
        }
        // Remove all receptors and catalytic centres referencing the substrate.
        self.receptors
            .retain(|r| !r.referes_to_substrate(substrate_index));
        // Adjust the indices of the remaining pointers.
        for receptor in &mut self.receptors {
            receptor.adjust_indices(substrate_index);
        }
        // Remove the substrate.
        self.substrates.remove(substrate_index)
    }

    /// Fuse this `Gene` with the specified `Gene` if possible and return the fusion.
    /// This function will fail if the fusion would overflow the internal [`Substrate`] or
    /// [`Receptor`] vector.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Receptor`]: ../protein/struct.Receptor.html
    pub fn fuse(&self, other_gene: &Gene<R, S, T>) -> Option<Gene<R, S, T>> {
        // Prevent overflow of substrates.
        if let None = self
            .substrates
            .len()
            .checked_add(other_gene.substrates.len())
        {
            return None;
        };
        // Prevent overflow of receptors.
        if let None = self.receptors.len().checked_add(other_gene.receptors.len()) {
            return None;
        };
        let mut fusion_gene = self.duplicate();
        fusion_gene
            .substrates
            .append(&mut other_gene.substrates.clone());
        let mut receptors = other_gene.receptors.clone();
        for mut receptor in &mut receptors {
            // Update all substrate indices.
            receptor.substrates = receptor
                .substrates
                .iter()
                .map(|i| i + self.number_of_substrates().get())
                .collect();
            receptor.triggers = receptor
                .triggers
                .iter()
                .map(|i| i + self.number_of_substrates().get())
                .collect();
            receptor.enzyme.educts = receptor
                .enzyme
                .educts
                .iter()
                .map(|i| i + self.number_of_substrates().get())
                .collect();
            receptor.enzyme.products = receptor
                .enzyme
                .products
                .iter()
                .map(|i| i + self.number_of_substrates().get())
                .collect();
        }
        fusion_gene.receptors.append(&mut receptors);
        Some(fusion_gene)
    }

    /// Duplicates the `Gene` and all its contents.
    pub fn duplicate(&self) -> Self {
        // At the moment this is just a wrapper for cloning.
        self.clone()
    }

    /// Creates a random [`GenomicCatalyticCentre`] specific to this `Gene`.
    pub fn random_catalytic_centre(&self) -> GenomicCatalyticCentre<R, S, T> {
        let reaction = R::random();
        let educts = (0..reaction.get_educt_number())
            .map(|_| self.get_random_substrate())
            .collect();
        let products = (0..reaction.get_product_number())
            .map(|_| self.get_random_substrate())
            .collect();
        GenomicCatalyticCentre {
            phantom_s: PhantomData,
            phantom_t: PhantomData,
            educts,
            products,
            reaction,
        }
    }

    /// Creates a random [`GenomicReceptor`] specific to this `Gene`.
    pub fn random_receptor(&self) -> GenomicReceptor<R, S, T> {
        let state = S::random();
        let enzyme = self.random_catalytic_centre();
        let substrates = (0..state.get_substrate_number())
            .map(|_| self.get_random_substrate())
            .collect();
        let triggers = vec![self.get_random_substrate()];
        GenomicReceptor {
            phantom_r: PhantomData,
            phantom_t: PhantomData,
            triggers,
            substrates,
            state,
            enzyme,
        }
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information + Default> Default for Gene<R, S, T> {
    fn default() -> Self {
        Gene {
            phantom_r: PhantomData,
            phantom_s: PhantomData,
            phantom_t: PhantomData,
            substrates: vec![T::default()],
            receptors: Vec::default(),
        }
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> CrossOver for Gene<R, S, T> {
    fn is_similar(&self, other: &Self) -> bool {
        other.number_of_substrates() == self.number_of_substrates()
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let substrates = self.substrates.cross_over(&other.substrates);
            let receptors = self.receptors.cross_over(&other.receptors);
            Gene {
                phantom_r: PhantomData,
                phantom_s: PhantomData,
                phantom_t: PhantomData,
                substrates,
                receptors,
            }
        } else {
            // If the two genes are not similar return a random one.
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

/// Adjust the index of [`Substrate`] pointers after removal of a [`Substrate`] from a [`Gene`].
///
/// # Parameters
///
/// * `current_index` - the index of the [`Substrate`] pointer before removal of [`Substrate`]
/// at `removed_index`
/// * `removed_index` - the index of the removed [`Substrate`]
///
/// # Panics
///
/// If `current_index` equals `removed_index`.
///
/// [`Gene`]: ./struct.Gene.html
/// [`Substrate`]: ../protein/struct.Substrate.html
fn adjust_index(current_index: usize, removed_index: usize) -> usize {
    if current_index < removed_index {
        current_index
    } else if current_index > removed_index {
        current_index - 1
    } else {
        panic!(
            "Index {} should have been removed, but was still passed as current index.",
            current_index
        );
    }
}

/// A `GenomicReceptor` represents the information of an actual [`Receptor`] that
/// is triggered upon specific [`Substrate`] changes and compares [`Substrate`], eventually
/// performing a [`Reaction`]. It is contained within a [`Gene`].
///
/// [`Receptor`]: ../protein/struct.Receptor.html
/// [`Gene`]: ./struct.Gene.html
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenomicReceptor<R, S, T> {
    phantom_r: PhantomData<R>,
    phantom_t: PhantomData<T>,
    // TODO: Triggers should be a set, as it does not make any sense,
    /// that a substrate can trigger a reaction twice.
    triggers: Vec<usize>,
    substrates: Vec<usize>,
    state: S,
    enzyme: GenomicCatalyticCentre<R, S, T>,
}

impl<R: Reaction<T>, S: State<T>, T: Information> GenomicReceptor<R, S, T> {
    /// Creates a `GenomicReceptor` encoding an actual [`Receptor`] detecting the specified [`State`] of its substrates and triggering
    /// the reaction encoded in its [`GenomicCatalyticCentre`].
    ///
    /// # Parameters
    ///
    /// * `triggers` - the [`Substrate`] indices in the respective containing [`Gene`] that trigger the Receptor upon change
    /// * `substrates` - the [`Substrate`] indices in the respective containing [`Gene`] the encoded [`State`] should check
    /// * `state` - the [`State`] to check for
    /// * `enzyme` - the [`GenomicCatalyticCentre`] to trigger if the [`State`] is appropriate
    ///
    /// # Panics
    ///
    /// If the number of [`Substrate`]s is not exactly equal to the one
    /// required by the [`State`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Receptor`]: ../protein/struct.Receptor.html
    /// [`State`]: ../chemistry/struct.State.html
    /// [`Gene`]: ./struct.Gene.html
    /// [`GenomicCatalyticCentre`]: ./struct.GenomicCatalyticCentre.html
    pub fn new(
        triggers: Vec<usize>,
        substrates: Vec<usize>,
        state: S,
        enzyme: GenomicCatalyticCentre<R, S, T>,
    ) -> Self {
        assert_eq!(substrates.len(), state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            state, state.get_substrate_number(), substrates.len());
        GenomicReceptor {
            phantom_r: PhantomData,
            phantom_t: PhantomData,
            triggers,
            substrates,
            state,
            enzyme,
        }
    }

    /// Adds a triggering [`Substrate`] to this `GenomicReceptor` if possible and returns the index of the
    /// new [`GenomicReceptor`].
    /// This function will fail if the underlying vector would overflow due to the addition.
    ///
    /// # Parameters
    ///
    /// * `trigger` - the triggering [`Substrate`] to add to the `GenomicReceptor`
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn add_trigger(&mut self, trigger: usize) -> Option<usize> {
        if let Some(new_index) = self.triggers.len().checked_add(1) {
            self.triggers.push(trigger);
            Some(new_index)
        } else {
            None
        }
    }

    /// Replaces the [`State`] of this `GenomicReceptor`.
    ///
    /// # Parameters
    ///
    /// * `state` - the [`State`] to check for
    /// * `substrates` - the [`Substrate`] indices in the respective containing [`Gene`]
    /// the encoded [`State`] should check
    ///
    /// # Panics
    ///
    /// If the number of [`Substrate`]s is not exactly equal to the one
    /// required by the [`State`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`State`]: ../chemistry/struct.State.html
    /// [`Gene`]: ./struct.Gene.html
    pub fn replace_state(&mut self, state: S, substrates: Vec<usize>) {
        assert_eq!(substrates.len(), state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            state, state.get_substrate_number(), substrates.len());
        self.state = state;
        self.substrates = substrates;
    }

    /// Replaces the [`Substrate`]s of this `GenomicReceptor`.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the [`Substrate`] indices in the respective containing [`Gene`]
    /// the encoded [`State`] should check
    ///
    /// # Panics
    ///
    /// If the number of [`Substrate`]s is not exactly equal to the one
    /// required by the [`State`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`State`]: ../chemistry/struct.State.html
    /// [`Gene`]: ./struct.Gene.html
    pub fn replace_substrates(&mut self, substrates: Vec<usize>) {
        assert_eq!(substrates.len(), self.state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            self.state, self.state.get_substrate_number(), substrates.len());
        self.substrates = substrates;
    }

    /// Replaces the [`Substrate`] at the specified index.
    ///
    /// # Parameters
    ///
    /// * `index` - the index of the [`Substrate`] reference to replace
    /// * `substrate` - the new [`Substrate`] reference value
    ///
    /// # Panics
    ///
    /// If the `index` is out of bounds.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn replace_substrate(&mut self, index: usize, substrate: usize) {
        self.substrates[index] = substrate;
    }

    /// Replaces the [`GenomicCatalyticCentre`].
    ///
    /// # Parameters
    ///
    /// * `enzyme` - the [`GenomicCatalyticCentre`] to trigger if the [`State`] is appropriate
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`GenomicCatalyticCentre`]: ./struct.GenomicCatalyticCentre.html
    pub fn replace_enzyme(&mut self, enzyme: GenomicCatalyticCentre<R, S, T>) {
        self.enzyme = enzyme;
    }

    /// Returns the [`GenomicCatalyticCentre`].
    ///
    /// [`GenomicCatalyticCentre`]: ./struct.GenomicCatalyticCentre.html
    pub fn enzyme(&self) -> &GenomicCatalyticCentre<R, S, T> {
        &self.enzyme
    }

    /// Returns the [`State`] that triggers this receptor.
    ///
    /// [`State`]: ../chemistry/trait.State.htm
    pub fn state(&self) -> &S {
        &self.state
    }

    /// Returns the substrates of this receptor.
    pub fn substrates(&self) -> &Vec<usize> {
        &self.substrates
    }

    /// Returns the [`GenomicCatalyticCentre`] as mutable.
    ///
    /// [`GenomicCatalyticCentre`]: ./struct.GenomicCatalyticCentre.html
    pub fn enzyme_mut(&mut self) -> &mut GenomicCatalyticCentre<R, S, T> {
        &mut self.enzyme
    }

    /// Removes a triggering [`Substrate`] from this `GenomicReceptor` and returns the
    /// removed [`Substrate`].
    ///
    /// # Parameters
    ///
    /// * `trigger` - the index of the triggering [`Substrate`] to remove from the `GenomicReceptor`
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn remove_trigger(&mut self, trigger: usize) -> usize {
        self.triggers.remove(trigger)
    }

    /// Returns the index of a random triggering [`Substrate`] of this `GenomicReceptor`
    /// if there are any triggers.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_random_trigger(&self) -> Option<usize> {
        if self.triggers.len() > 0 {
            Some(thread_rng().gen_range(0, self.triggers.len()))
        } else {
            None
        }
    }

    /// Returns the index of a random [`Substrate`] of this `GenomicReceptor`
    /// if there are any.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_random_substrate(&self) -> Option<usize> {
        if self.substrates.len() > 0 {
            Some(thread_rng().gen_range(0, self.substrates.len()))
        } else {
            None
        }
    }

    /// Checks wether this `GenomicReceptor` contains any reference to the specified [`Substrate`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn referes_to_substrate(&self, substrate_index: usize) -> bool {
        self.enzyme.referes_to_substrate(substrate_index)
            || self
                .triggers
                .iter()
                .chain(self.substrates.iter())
                .any(|element| element == &substrate_index)
    }

    /// Adjust the indices of [`Substrate`] pointers after removal of a [`Substrate`].
    ///
    /// # Parameters
    ///
    /// * `removed_index` - the index of the removed [`Substrate`]
    ///
    /// # Panics
    ///
    /// If `removed_index` is still present as index pointer.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_indices(&mut self, removed_index: usize) {
        for trigger in &mut self.triggers {
            *trigger = adjust_index(*trigger, removed_index);
        }
        for substrate in &mut self.substrates {
            *substrate = adjust_index(*substrate, removed_index);
        }
        self.enzyme.adjust_indices(removed_index);
    }

    /// Translates the `GenomicReceptor` into a [`Receptor`] and directly associates them with their
    /// triggering [`Substrate`]s.
    ///
    /// # Parameters
    /// * `gene_index` - the index of the [`Gene`] containing this `Receptor` inside the
    /// [`Genome`].
    /// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`]s of the containing [`Genome`]
    ///
    /// # Panics
    ///
    /// If the `substrate_lookup` map does not contain one of the requested [`Substrate`]s.
    ///
    /// [`Receptor`]: ../protein/struct.Receptor.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`Gene`]: ./struct.Gene.html
    fn translate(
        &self,
        gene_index: usize,
        substrate_lookup: &HashMap<GeneSubstrate, Rc<RefCell<Substrate<R, S, T>>>>,
    ) {
        let substrates = self
            .substrates
            .iter()
            .map(|substrate_index| GeneSubstrate::new(gene_index, *substrate_index))
            .map(|gene_substrate| {
                substrate_lookup.get(&gene_substrate).expect(&format!(
                    "The substrate lookup map did not contain {:?}.",
                    &gene_substrate
                ))
            })
            .map(|strong| Rc::downgrade(strong))
            .collect();
        let enzyme = self.enzyme.translate(gene_index, substrate_lookup);
        let state = self.state.clone();
        let receptor = Rc::new(Receptor::new(substrates, state, enzyme));
        for trigger_index in self.triggers.iter() {
            let gene_substrate = GeneSubstrate::new(gene_index, *trigger_index);
            substrate_lookup
                .get(&gene_substrate)
                .expect(&format!("The substrate lookup map did not contain {:?}.", &gene_substrate))
                .borrow_mut()
                .add_receptor(receptor.clone());
        }
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> CrossOver for GenomicReceptor<R, S, T> {
    fn is_similar(&self, other: &Self) -> bool {
        self.state.is_similar(&other.state)
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let triggers = self.triggers.cross_over(&other.triggers);
            let substrates = self.substrates.cross_over(&other.substrates);
            let state = self.state.cross_over(&other.state);
            let enzyme = self.enzyme.cross_over(&other.enzyme);
            GenomicReceptor {
                phantom_r: PhantomData,
                phantom_t: PhantomData,
                triggers,
                substrates,
                state,
                enzyme,
            }
        } else {
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

/// A `GenomicCatalyticCentre` represents the information of an actual [`CatalyticCentre`] that
/// produces products from educt [`Substrate`]s
/// by performing a [`Reaction`]. It is contained within a [`Gene`].
///
/// [`CatalyticCentre`]: ../protein/struct.CatalyticCentre.html
/// [`Gene`]: ./struct.Gene.html
/// [`Substrate`]: ../protein/struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenomicCatalyticCentre<R, S, T> {
    phantom_s: PhantomData<S>,
    phantom_t: PhantomData<T>,
    educts: Vec<usize>,
    products: Vec<usize>,
    reaction: R,
}

impl<R: Reaction<T>, S: State<T>, T: Information> GenomicCatalyticCentre<R, S, T> {
    /// Creates a new `GenomicCatalyticCentre` containing the information to produce
    /// an actual [`CatalyticCentre`] via the containing [`Gene`].
    ///
    /// # Parameters
    ///
    /// * `educts` - the educt [`Substrate`]s indices in the respective containing [`Gene`]
    /// * `products` - the product [`Substrate`]s indices in the respective containing [`Gene`]
    /// * `reaction` - the [`Reaction`] to catalyse
    ///
    /// # Panics
    ///
    /// If the number of educt or product [`Substrate`]s is not exactly equal to the one
    /// required by the [`Reaction`].
    ///
    /// [`CatalyticCentre`]: ../protein/struct.CatalyticCentre.html
    /// [`Gene`]: ./struct.Gene.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Reaction`]: ../chemistry/trait.Reaction.html
    pub fn new(educts: Vec<usize>, products: Vec<usize>, reaction: R) -> Self {
        assert_eq!(
            educts.len(),
            reaction.get_educt_number(),
            "The number of required educts for reaction {:?} is {}, but {} educts were supplied.",
            reaction,
            reaction.get_educt_number(),
            educts.len()
        );
        assert_eq!(products.len(), reaction.get_product_number(),
            "The number of required products for reaction {:?} is {}, but {} products were supplied.",
            reaction, reaction.get_product_number(), products.len());
        GenomicCatalyticCentre {
            phantom_s: PhantomData,
            phantom_t: PhantomData,
            educts,
            products,
            reaction,
        }
    }

    /// Returns the [`Reaction`] catalysed by this `GenomicCatalyticCentre`.
    ///
    /// [`Reaction`]: ../chemistry/trait.Reaction.html
    pub fn reaction(&self) -> &R {
        &self.reaction
    }

    /// Returns the reaction educts.
    pub fn educts(&self) -> &Vec<usize> {
        &self.educts
    }

    /// Returns the reaction educts.
    pub fn products(&self) -> &Vec<usize> {
        &self.products
    }

    /// Checks wether this `GenomicCatalyticCentre` contains any reference to the specified [`Substrate`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn referes_to_substrate(&self, substrate_index: usize) -> bool {
        self.educts
            .iter()
            .chain(self.products.iter())
            .any(|element| element == &substrate_index)
    }

    /// Adjust the indices of [`Substrate`] pointers after removal of a [`Substrate`].
    ///
    /// # Parameters
    ///
    /// * `removed_index` - the index of the removed [`Substrate`]
    ///
    /// # Panics
    ///
    /// If `removed_index` is still present as index pointer.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_indices(&mut self, removed_index: usize) {
        for educt in &mut self.educts {
            *educt = adjust_index(*educt, removed_index);
        }
        for product in &mut self.products {
            *product = adjust_index(*product, removed_index);
        }
    }

    /// Returns the index of a random educt [`Substrate`] of this `GenomicCatalyticCentre`
    /// if there are any.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_random_educt(&self) -> Option<usize> {
        if self.educts.len() > 0 {
            Some(thread_rng().gen_range(0, self.educts.len()))
        } else {
            None
        }
    }

    /// Returns the index of a random product [`Substrate`] of this `GenomicCatalyticCentre`
    /// if there are any.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_random_product(&self) -> Option<usize> {
        if self.products.len() > 0 {
            Some(thread_rng().gen_range(0, self.products.len()))
        } else {
            None
        }
    }

    /// Replaces a single educt.
    ///
    /// # Parameters
    ///
    /// * `index` - the index of the educt to replace
    /// * `educt` - the new educt reference
    ///
    /// # Panics
    ///
    /// If the `index` is out of bounds.
    pub fn replace_educt(&mut self, index: usize, educt: usize) {
        self.educts[index] = educt;
    }

    /// Replaces a single product.
    ///
    /// # Parameters
    ///
    /// * `index` - the index of the educt to replace
    /// * `product` - the new product reference
    ///
    /// # Panics
    ///
    /// If the `index` is out of bounds.
    pub fn replace_product(&mut self, index: usize, product: usize) {
        self.products[index] = product;
    }

    /// Translates the `GenomicCatalyticCentre` into a [`CatalyticCentre`].
    ///
    /// # Parameters
    /// * `gene_index` - the index of the [`Gene`] containing this `GenomicCatalyticCentre` inside the
    /// [`Genome`].
    /// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`]s of the containing [`Genome`]
    ///
    /// # Panics
    ///
    /// If the `substrate_lookup` map does not contain one of the requested [`Substrate`]s.
    ///
    /// [`CatalyticCentre`]: ../protein/struct.CatalyticCentre.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`Gene`]: ./struct.Gene.html
    fn translate(
        &self,
        gene_index: usize,
        substrate_lookup: &HashMap<GeneSubstrate, Rc<RefCell<Substrate<R, S, T>>>>,
    ) -> CatalyticCentre<R, S, T> {
        let educts = self
            .educts
            .iter()
            .map(|substrate_index| GeneSubstrate::new(gene_index, *substrate_index))
            .map(|gene_substrate| {
                substrate_lookup.get(&gene_substrate).expect(&format!(
                    "The substrate lookup map did not contain {:?}.",
                    &gene_substrate
                ))
            })
            .map(|strong| Rc::downgrade(strong))
            .collect();
        let products = self
            .products
            .iter()
            .map(|substrate_index| GeneSubstrate::new(gene_index, *substrate_index))
            .map(|gene_substrate| {
                substrate_lookup.get(&gene_substrate).expect(&format!(
                    "The substrate lookup map did not contain {:?}.",
                    &gene_substrate
                ))
            })
            .map(|strong| Rc::downgrade(strong))
            .collect();
        let reaction = self.reaction.clone();
        CatalyticCentre::new(educts, products, reaction)
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> CrossOver for GenomicCatalyticCentre<R, S, T> {
    fn is_similar(&self, other: &Self) -> bool {
        self.reaction.is_similar(&other.reaction)
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let educts = self.educts.cross_over(&other.educts);
            let products = self.products.cross_over(&other.products);
            let reaction = self.reaction.cross_over(&other.reaction);
            GenomicCatalyticCentre {
                phantom_s: PhantomData,
                phantom_t: PhantomData,
                educts,
                products,
                reaction,
            }
        } else {
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

/// A `GenomicInputSensor` represents the information of an actual
/// [`InputSensor`](oben::evolution::protein::InputSensor)
/// that provides input
/// [`Substrate`](oben::evolution::protein::Substrate)s.
/// It is contained within a
/// [`Gene`](oben::evolution::gene::Gene).
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenomicInputSensor<
    ReactionType,
    StateType,
    InformationType,
    InputElementType,
    InputSensorType,
> {
    phantom_reaction: PhantomData<ReactionType>,
    phantom_state: PhantomData<StateType>,
    phantom_information: PhantomData<InformationType>,
    phantom_input_element: PhantomData<InputElementType>,
    input_substrates: Vec<Option<GeneSubstrate>>,
    feedback_substrates: Vec<Option<GeneSubstrate>>,
    input: InputSensorType,
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
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
    >
    GenomicInputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Get the number of associated input [`Substrate`](oben::evolution::protein::Substrate)s.
    pub fn number_of_associated_inputs(&self) -> usize {
        self.input_substrates.iter().filter(|i| i.is_some()).count()
    }

    /// Adjust the input [`GeneSubstrate`] references after the [`Substrate`] of a contained
    /// [`Gene`] was removed.
    ///
    /// # Parameters
    ///
    /// * `removed_substrate` - an index based pointer to the removed [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_after_gene_substrate_removal(&mut self, removed_substrate: GeneSubstrate) {
        Self::adjust_substrates(&mut self.input_substrates, removed_substrate);
        Self::adjust_substrates(&mut self.feedback_substrates, removed_substrate);
    }

    fn adjust_substrates(
        substrates: &mut Vec<Option<GeneSubstrate>>,
        removed_substrate: GeneSubstrate,
    ) {
        for substrate in substrates {
            *substrate = Self::adjust_substrate(*substrate, removed_substrate);
        }
    }

    fn adjust_substrate(
        substrate_option: Option<GeneSubstrate>,
        removed_substrate: GeneSubstrate,
    ) -> Option<GeneSubstrate> {
        substrate_option.and_then(|substrate| {
            if substrate == removed_substrate {
                None
            } else if substrate.is_gene(removed_substrate.gene()) {
                Some(GeneSubstrate::new(
                    substrate.gene(),
                    adjust_index(substrate.substrate(), removed_substrate.substrate()),
                ))
            } else {
                Some(substrate)
            }
        })
    }

    /// Validates the [`GeneSubstrate`] references after a recombination
    /// event and removes invalid ones.
    ///
    /// # Parameters
    ///
    /// * `genes` - the [`Gene`]s of the [`Genome`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    fn validate(&mut self, genes: &Vec<Gene<ReactionType, StateType, InformationType>>) {
        Self::validate_substrates(&mut self.input_substrates, genes);
        Self::validate_substrates(&mut self.feedback_substrates, genes);
    }

    fn validate_substrates(
        substrates: &mut Vec<Option<GeneSubstrate>>,
        genes: &Vec<Gene<ReactionType, StateType, InformationType>>,
    ) {
        let mut duplicate_checker = HashSet::new();
        for substrate_option in substrates {
            if let Some(substrate) = substrate_option {
                // Remove the association if it points to a duplicate or an invalid substrate.
                if !duplicate_checker.insert(*substrate)
                    || !Genome::<
                        ReactionType,
                        StateType,
                        InformationType,
                        InputElementType,
                        InputSensorType,
                    >::has_substrate(genes, substrate)
                {
                    *substrate_option = None;
                }
                // Otherwise, leave it untouched.
            }
        }
    }

    /// Translates the `GenomicInputSensor` into an [`InputSensor`](oben::evolution::protein::InputSensor).
    ///
    /// # Parameters
    /// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`](oben::evolution::protein::Substrate)s
    /// of the containing [`Genome`](oben::evolution::gene::Genome)
    ///
    /// # Panics
    ///
    /// If the `substrate_lookup` map does not contain one of the requested [`Substrate`](oben::evolution::protein::Substrate)s.
    fn translate(
        &self,
        substrate_lookup: &HashMap<
            GeneSubstrate,
            Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
        >,
    ) -> InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    {
        let input_substrates = Self::translate_substrates(&self.input_substrates, substrate_lookup);
        let feedback_substrates =
            Self::translate_substrates(&self.feedback_substrates, substrate_lookup);
        let input = self.input.clone();
        InputSensor::new(input, input_substrates, feedback_substrates)
    }

    pub fn remove_associations_with_gene(&mut self, gene_index: usize) {
        Self::remove_associations_with_gene_from_substrates(&mut self.input_substrates, gene_index);
        Self::remove_associations_with_gene_from_substrates(
            &mut self.feedback_substrates,
            gene_index,
        );
    }

    fn remove_associations_with_gene_from_substrates(
        substrates: &mut Vec<Option<GeneSubstrate>>,
        gene_index: usize,
    ) {
        for substrate_option in substrates {
            if let Some(substrate) = substrate_option {
                if substrate.is_gene(gene_index) {
                    *substrate_option = None;
                }
            }
        }
    }

    /// Translates the [`Substrate`](oben::evolution::protein::Substrate)s of a `GenomicInputSensor`
    ///
    /// # Parameters
    /// * `substrates` - the [`Substrate`](oben::evolution::protein::Substrate)s to transcribe
    /// * `substrate_lookup` - a map for lookup of all the translated [`Substrate`](oben::evolution::protein::Substrate)s
    /// of the containing [`Genome`](oben::evolution::gene::Genome)
    ///
    /// # Panics
    ///
    /// If the `substrate_lookup` map does not contain one of the requested [`Substrate`](oben::evolution::protein::Substrate)s.
    fn translate_substrates(
        substrates: &Vec<Option<GeneSubstrate>>,
        substrate_lookup: &HashMap<
            GeneSubstrate,
            Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
        >,
    ) -> Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        substrates
            .iter()
            .map(Option::as_ref)
            .map(|gene_substrate_option| {
                gene_substrate_option.map(|gene_substrate| {
                    substrate_lookup.get(gene_substrate).expect(&format!(
                        "The substrate lookup map did not contain {:?}.",
                        gene_substrate
                    ))
                })
            })
            .map(|strong_option| strong_option.map(|strong| Rc::downgrade(strong)))
            .collect()
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
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType> + Default,
    > Default
    for GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    >
{
    fn default() -> Self {
        GenomicInputSensor {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_input_element: PhantomData,
            input_substrates: Vec::default(),
            feedback_substrates: Vec::default(),
            input: InputSensorType::default(),
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
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
    > CrossOver
    for GenomicInputSensor<
        ReactionType,
        StateType,
        InformationType,
        InputElementType,
        InputSensorType,
    >
{
    fn is_similar(&self, other: &Self) -> bool {
        self.input.is_similar(&other.input)
    }

    fn cross_over(&self, other: &Self) -> Self {
        if self.is_similar(other) {
            let input_substrates = self.input_substrates.cross_over(&other.input_substrates);
            let feedback_substrates = self
                .feedback_substrates
                .cross_over(&other.feedback_substrates);
            let input = self.input.cross_over(&other.input);
            GenomicInputSensor {
                phantom_reaction: PhantomData,
                phantom_state: PhantomData,
                phantom_information: PhantomData,
                phantom_input_element: PhantomData,
                input_substrates,
                feedback_substrates,
                input,
            }
        } else {
            do_a_or_b(|| self.clone(), || other.clone())
        }
    }
}

pub trait GenomeMutation<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
    InputElementType: Clone + std::fmt::Debug + PartialEq + Send + Sync + CrossOver + Serialize + DeserializeOwned,
    InputSensorType: Input<InputElementType, InformationType>,
>: Sized + Send + Sync
{
    /// Generates a mutated version of the specified [`Genome`] based on the kind of
    /// `GenomeMutation`.
    ///
    /// This will fail if the mutated [`Genome`] would be identical to the input or if the
    /// mutation is impossible for the specified [`Genome`].
    ///
    /// # Parameters
    ///
    /// `genome` - the base [`Genome`] to generate a mutated version of
    ///
    /// [`Genome`]: ./struct.Genome.html
    fn mutate(
        &self,
        genome: &Genome<
            ReactionType,
            StateType,
            InformationType,
            InputElementType,
            InputSensorType,
        >,
    ) -> Option<Genome<ReactionType, StateType, InformationType, InputElementType, InputSensorType>>;

    /// Creates a random `GenomeMutation`.
    fn random() -> Self;
}

/// The `CrossOver` trait allows for genetic elements to be compared for similarity and be
/// recombined.
pub trait CrossOver {
    /// Checks if two genetic components are similar enough for cross-over.
    ///
    /// # Parameters
    ///
    /// * `other` - a matching genetic component of the other individual
    fn is_similar(&self, other: &Self) -> bool;

    /// Creates a randomly recombined genetic component based on two initial components.
    ///
    /// # Parameters
    ///
    /// * `other` - a matching genetic component of the other individual
    fn cross_over(&self, other: &Self) -> Self;
}

// This implementation improves the usability of the `CrossOver` trait for genomic elements.
impl CrossOver for usize {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        a_or_b(*self, *other)
    }
}

// This implementation improves the usability of the `CrossOver` trait for genomic elements.
impl<T> CrossOver for Option<T>
where
    T: CrossOver + Sized + Clone,
{
    fn is_similar(&self, other: &Self) -> bool {
        match (self, other) {
            (None, None) => true,
            (Some(a), Some(b)) => a.is_similar(b),
            _ => false,
        }
    }

    fn cross_over(&self, other: &Self) -> Self {
        match (self, other) {
            (None, None) => None,
            (Some(a), Some(b)) => Some(a.cross_over(b)),
            _ => do_a_or_b(|| self.clone(), || other.clone()),
        }
    }
}

// This implementation improves the usability of the `CrossOver` trait for genomic elements.
impl<T> CrossOver for Vec<T>
where
    T: CrossOver + Clone,
{
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        let mut recombined_vec = Vec::new();
        // Determine shorter and longer vector.
        let short;
        let long;
        if self.len() >= other.len() {
            short = other;
            long = self;
        } else {
            short = self;
            long = other;
        }
        for (index, element_long) in long.iter().enumerate() {
            if let Some(element_short) = short.get(index) {
                // As long as there are corresponding elements in the shorter vector,
                // perform recombination.
                recombined_vec.push(element_long.cross_over(element_short));
            } else {
                // If there are corresponding elements in the shorter vector
                // randomly add the element from the longer vector or skip it.
                do_a_or_b(|| (), || recombined_vec.push(element_long.clone()));
            }
        }
        recombined_vec
    }
}

#[cfg(test)]
mod tests;
