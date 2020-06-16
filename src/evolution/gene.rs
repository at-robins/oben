//! The `gene` module contains information structures to describe,
//! replicate and mutate the evolutionary network.
extern crate bitvec;
extern crate rand;
extern crate serde;

use super::chemistry::{Reaction, State};
use super::population::Organism;
use super::protein::{CatalyticCentre, Receptor, Substrate};
use bitvec::{boxed::BitBox, order::Local, vec::BitVec};
use rand::{distributions::{Distribution, Standard}, thread_rng, Rng};
use std::num::NonZeroUsize;
use std::error::Error;
use std::fs::File;
use std::path::Path;
use serde::{Deserialize, Serialize};
use core::cell::RefCell;
use std::rc::Rc;
use std::collections::HashMap;

/// The minimal length in byte of a randomly created binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
const RANDOM_SUBSTRATE_MIN_LENGTH: usize = 0;
/// The maximal length in byte of a randomly created binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
const RANDOM_SUBSTRATE_MAX_LENGTH: usize = 64;

/// A `Genome` is a collection of individual [`Gene`]s and associations between them.
/// A `Genome` is required to consist of 1 or more genes.
///
/// [`Gene`]: ./struct.Gene.html
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct Genome {
    input: Vec<Option<GeneSubstrate>>,
    output: Vec<Option<GeneSubstrate>>,
    genes: Vec<Gene>,
    associations: Vec<GeneAssociation>,
}

impl Genome {
    /// Get the number of [`Gene`]s in this `Genome`.
    /// A `Genome` must encode 1 or more [`Gene`]s.
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn number_of_genes(&self) -> NonZeroUsize {
        NonZeroUsize::new(self.genes.len()).expect("No gene is encoded by this genome. This is forbidden by the contract.")
    }

    /// Get the number of input [`Substrate`]s in this `Genome`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn number_of_inputs(&self) -> usize {
        self.input.len()
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

    /// Sets the input association at the specified index to the specified value
    /// and returns the previous value.
    ///
    /// # Parameters
    ///
    /// * `input_index` - the index of the input association to changes
    /// * `input_value` - the new value for the specified input association
    ///
    /// # Panics
    ///
    /// If the index is out of bounds.
    fn set_input(&mut self, input_index: usize, input_value: Option<GeneSubstrate>) -> Option<GeneSubstrate> {
        if input_index >= self.input.len() {
            panic!("The input vector of this genome is of length {}, but insertion of element {:#?} at index {} was attempted.", self.number_of_inputs(), input_value, input_index);
        }
        let old_value = self.input[input_index];
        self.input[input_index] = input_value;
        old_value
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
    fn set_output(&mut self, output_index: usize, output_value: Option<GeneSubstrate>) -> Option<GeneSubstrate> {
        if output_index >= self.output.len() {
            panic!("The output vector of this genome is of length {}, but insertion of element {:#?} at index {} was attempted.", self.number_of_outputs(), output_value, output_index);
        }
        let old_value = self.output[output_index];
        self.output[output_index] = output_value;
        old_value
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
    fn add_gene(&mut self, gene: Gene) -> Option<usize>{
        if let Some(new_index) = self.genes.len().checked_add(1) {
            self.genes.push(gene);
            Some(new_index)
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
    fn remove_gene(&mut self, gene: usize) -> Gene {
        if self.number_of_genes().get() <= 1 {
            panic!("A genome needs to contain at least one gene, so no gene can be removed.");
        }
        // Remove all inputs pointing to the removed gene.
        for input_value in &mut self.input {
            if input_value.and_then(|input| Some(input.gene == gene)).unwrap_or(false) {
                *input_value = None;
            }
        }
        // Remove all outputs pointing to the removed gene.
        for output_value in &mut self.output {
            if output_value.and_then(|output| Some(output.gene == gene)).unwrap_or(false) {
                *output_value = None;
            }
        }
        // Remove all substrate associations pointing to the removed gene.
        for association in &mut self.associations {
            association.remove_associated_gene(gene);
        }
        self.genes.remove(gene)
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
    pub fn duplicate_gene(&self, gene: usize) -> Gene {
        self.genes[gene].duplicate()
    }

    /// Duplicates the [`Gene`] at the specified index and adds it directly to the genome
    /// if possible.
    /// This function will fail if the underlying vector would overflow due to the duplication.
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
    fn duplicate_gene_internal(&mut self, gene: usize) -> Option<usize> {
        self.add_gene(self.genes[gene].duplicate())
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
    pub fn get_gene(&self, gene: usize) -> &Gene {
        &self.genes[gene]
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
    fn get_gene_substrates(&self, gene: usize) -> Vec<GeneSubstrate, Rc<RefCell<BitBox<Local, u8>>>> {
        let gene = &self.genes[gene];
        (0..gene.number_of_substrates())
            .map()
            .map()
            .collect()
    }

    /// Returns the index of a random [`Gene`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn get_random_gene(&self) -> usize {
        thread_rng().gen_range(0, self.number_of_genes().get())
    }

    /// Returns the index of a random input [`GeneSubstrate`] if there is any.
    ///
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    pub fn get_random_input(&self) -> Option<usize> {
        if self.number_of_inputs() > 0 {
            Some(thread_rng().gen_range(0, self.number_of_inputs()))
        } else {
            None
        }
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
    fn random_gene_substrate(&self) -> GeneSubstrate {
        let random_gene = self.get_random_gene();
        let random_substrate = self.get_gene(random_gene).get_random_substrate();
        GeneSubstrate{
            gene: random_gene,
            substrate: random_substrate
        }
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
    fn add_association(&mut self, association: GeneAssociation) -> Option<usize>{
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
    fn remove_association(&mut self, association: usize) -> GeneAssociation {
        self.associations.remove(association)
    }

    /// Adjust the input [`GeneSubstrate`] references after the binary [`Substrate`] of a contained
    /// [`Gene`] was removed.
    ///
    /// # Parameters
    ///
    /// * `removed_substrate` - an index based pointer to the removed [`Substrate`]
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_input_after_gene_substrate_removal(&mut self, removed_substrate: GeneSubstrate) {
        for input_ref in &mut self.input {
            if input_ref.is_some() {
                let mut current_reference = input_ref.unwrap();
                if current_reference == removed_substrate {
                    // Remove the input association if it points to the removed substrate.
                    *input_ref = None;
                } else if current_reference.is_gene(removed_substrate.gene) {
                    // Adjust the input association if it points to the gene from which the
                    // substrate was removed.
                    current_reference.substrate = Gene::adjust_index(current_reference.substrate, removed_substrate.substrate);
                    *input_ref = Some(current_reference);
                }
                // Otherwise, leave the input untouched.
            }
        }
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
                    current_reference.substrate = Gene::adjust_index(current_reference.substrate, removed_substrate.substrate);
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
    fn adjust_associations_after_gene_substrate_removal(&mut self, removed_substrate: GeneSubstrate) {
        for association in &mut self.associations {
            association.adjust_after_gene_substrate_removal(removed_substrate);
        }
    }

    /// Load a `Genome` from a JSON file if possible.
    /// An error will be returned if parsing the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file from which the `Genome` should be loaded
    pub fn load_from_file<P>(path_to_file: P) -> Result<Self, Box<dyn Error>> where P: AsRef<Path> {
       let file = File::open(path_to_file)?;
       let genome = serde_json::from_reader(file)?;
       Ok(genome)
    }

    /// Write a `Genome` to a JSON file if possible.
    /// An error will be returned if writing to the file failed.
    ///
    /// # Parameters
    ///
    /// * `path_to_file` - the JSON file the `Genome` should be written to
    pub fn write_to_file<P>(&self, path_to_file: P) -> Result<(), Box<dyn Error>> where P: AsRef<Path> {
       let file = File::create(path_to_file)?;
       let write_success = serde_json::to_writer(file, self)?;
       Ok(write_success)
    }

    pub fn translate(&self) -> Organism {
        let mut gene_substrate_map: HashMap<GeneSubstrate, Rc<RefCell<BitBox<Local, u8>>>> = HashMap::new();
        for gene_association in &self.associations {
            let genome_level_substrate = Rc::new(RefCell::new(gene_association.substrate.clone()));
            for gene_substrate in &gene_association.associations {
                gene_substrate_map.entry(gene_substrate.clone()).or_insert(genome_level_substrate.clone());
            }
        }

        let substrates = vec!();
        let input = vec!();
        let output = vec!();
        Organism{substrates, input, output}
    }
}

impl Default for Genome {
    fn default() -> Self {
        Genome {
            input: Vec::default(),
            output: Vec::default(),
            genes: vec!(Gene::default()),
            associations: Vec::default(),
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
struct GeneSubstrate {
    // index of gene inside the containing genome
    gene: usize,
    // index of substrate inside the respective gene
    substrate: usize,
}

impl GeneSubstrate {

    /// Checks wether this `GeneSubstrate` references the [`Gene`] at the specified index.
    ///
    /// # Parameters
    ///
    /// * `gene_index` - the index of the gene to check for
    ///
    /// [`Gene`]: ./struct.Gene.html
    fn is_gene(&self, gene_index: usize) -> bool {
        self.gene == gene_index
    }

    /// Checks wether this `GeneSubstrate` references the [`Substrate`] at the specified index.
    ///
    /// # Parameters
    ///
    /// * `substrate_index` - the index of the substrate to check for
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn is_substrate(&self, substrate_index: usize) -> bool {
        self.substrate == substrate_index
    }

}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
/// A `GeneAssociation` is a [`Substrate`] defined on [`Genome`] level that replaces [`Gene`] specific
/// [`Substrate`]s upon translation. This process is used to interconnect different [`Gene`]s on a
/// [`Genome`] level.
///
/// [`Gene`]: ./struct.Gene.html
/// [`Genome`]: ./struct.Genome.html
/// [`Substrate`]: ../protein/struct.Substrate.html
struct GeneAssociation {
    // substrate value defined in the genome and shared between genes
    substrate: BitBox<Local, u8>,
    // gene specific substrates pointing to the shared substrate
    associations: Vec<GeneSubstrate>,
}

impl GeneAssociation {
    /// Remove all associations with the specified [`Gene`].
    ///
    /// # Parameters
    ///
    /// * `gene` - the index of the gene to remove associations to
    ///
    /// [`Gene`]: ./struct.Gene.html
    fn remove_associated_gene(&mut self, gene: usize) {
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
    fn add_association(&mut self, association: GeneSubstrate) -> Option<usize>{
        if let Some(new_index) = self.associations.len().checked_add(1) {
            self.associations.push(association);
            Some(new_index)
        } else {
            None
        }
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
                association.substrate = Gene::adjust_index(association.substrate, removed_substrate.substrate);
            }
            // Otherwise, leave the input untouched.
        }
    }
}

/// A `Gene` is an immutable structure encoding a self-contained network, but without
/// explicite function. It can be transcribed into a functional protein network.
/// A `Gene` is required to encode at least 1 [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct Gene {
    substrates: Vec<BitBox<Local, u8>>,
    receptors: Vec<GenomicReceptor>,
}

impl Gene {
    /// Gets the number of [`Substrate`]s encoded by this `Gene`.
    /// A `Gene` must encode for 1 or more [`Substrate`]s.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn number_of_substrates(&self) -> NonZeroUsize {
        NonZeroUsize::new(self.substrates.len()).expect("No substrate is encoded by this gene. This is forbidden by the contract.")
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
    fn add_substrate(&mut self, substrate: BitBox<Local, u8>) -> Option<usize>{
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
    fn add_receptor(&mut self, receptor: GenomicReceptor) -> Option<usize>{
        if let Some(new_index) = self.receptors.len().checked_add(1) {
            self.receptors.push(receptor);
            Some(new_index)
        } else {
            None
        }
    }

    /// Removes the binary [`Substrate`]a the specified index from the `Gene` if possible and
    /// returns the removed [`Substrate`].
    ///
    /// # Parameters
    ///
    /// * `substrate_index` - the [`Substrate`]'s index to remove to the `Gene`
    ///
    /// # Panics
    ///
    /// If the specified [`Substrate`]'s index is out of bounds or there is just a single
    /// [`Substrate`] present, as a `Gene` needs at least 1 [`Substrate`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn remove_substrate(&mut self, substrate_index: usize) -> BitBox<Local, u8>{
        if self.number_of_substrates().get() <= 1 {
            panic!("A genome needs to contain at least one substrate, so no substrate can be removed.");
        }
        // Remove all receptors and catalytic centres referencing the substrate.
        self.receptors.retain(|r| !r.referes_to_substrate(substrate_index));
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
    pub fn fuse(&self, other_gene: &Gene) -> Option<Gene> {
        // Prevent overflow of substrates.
        if let None = self.substrates.len().checked_add(other_gene.substrates.len()) {
            return None;
        };
        // Prevent overflow of receptors.
        if let None = self.receptors.len().checked_add(other_gene.receptors.len()) {
            return None;
        };
        let mut fusion_gene = self.duplicate();
        fusion_gene.substrates.append(&mut other_gene.substrates.clone());
        let mut receptors = other_gene.receptors.clone();
        for mut receptor in &mut receptors {
            // Update all substrate indices.
            receptor.substrates = receptor.substrates.iter().map(|i| i + self.number_of_substrates().get()).collect();
            receptor.triggers = receptor.triggers.iter().map(|i| i + self.number_of_substrates().get()).collect();
            receptor.enzyme.educts = receptor.enzyme.educts.iter().map(|i| i + self.number_of_substrates().get()).collect();
            receptor.enzyme.products = receptor.enzyme.products.iter().map(|i| i + self.number_of_substrates().get()).collect();
        }
        fusion_gene.receptors.append(&mut receptors);
        Some(fusion_gene)
    }

    /// Duplicates the `Gene` and all its contents.
    pub fn duplicate(&self) -> Self {
        // At the moment this is just a wrapper for cloning.
        self.clone()
    }

    /// Adjust the index of [`Substrate`] pointers after removal of a [`Substrate`].
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
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn adjust_index(current_index: usize, removed_index: usize) -> usize {
        if current_index < removed_index {
            current_index
        } else if current_index > removed_index {
            current_index - 1
        } else {
            panic!("Index {} should have been removed, but was still passed as current index.", current_index);
        }
    }

    /// Creates a random [`GenomicCatalyticCentre`] specific to this `Gene`.
    fn random_catalytic_centre(&self) -> GenomicCatalyticCentre {
        let reaction: Reaction = rand::random();
        let educts = (0..reaction.get_educt_number()).map(|_| self.get_random_substrate()).collect();
        let products = (0..reaction.get_product_number()).map(|_| self.get_random_substrate()).collect();
        GenomicCatalyticCentre{educts, products, reaction}
    }

    /// Creates a random [`GenomicReceptor`] specific to this `Gene`.
    fn random_receptor(&self) -> GenomicReceptor {
        let state: State = rand::random();
        let enzyme = self.random_catalytic_centre();
        let substrates = (0..state.get_substrate_number()).map(|_| self.get_random_substrate()).collect();
        let triggers = vec!(self.get_random_substrate());
        GenomicReceptor{triggers, substrates, state, enzyme}
    }
}

impl Default for Gene {
    fn default() -> Self {
        Gene {
            substrates: vec!(BitBox::empty()),
            receptors: Vec::default(),
        }
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
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct GenomicReceptor {
    triggers: Vec<usize>,
    substrates: Vec<usize>,
    state: State,
    enzyme: GenomicCatalyticCentre,
}

impl GenomicReceptor {
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
    pub fn new(triggers: Vec<usize>, substrates: Vec<usize>, state: State, enzyme: GenomicCatalyticCentre) -> Self {
        assert_eq!(substrates.len(), state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            state, state.get_substrate_number(), substrates.len());
        GenomicReceptor{triggers, substrates, state, enzyme}
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
    fn add_trigger(&mut self, trigger: usize) -> Option<usize>{
        if let Some(new_index) = self.triggers.len().checked_add(1) {
            self.triggers.push(trigger);
            Some(new_index)
        } else {
            None
        }
    }

    /// Returns the index of a random triggering [`Substrate`] of this `GenomicReceptor`
    /// if there are any triggers.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn get_random_trigger(&self) -> Option<usize> {
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
    fn get_random_substrate(&self) -> Option<usize> {
        if self.substrates.len() > 0 {
            Some(thread_rng().gen_range(0, self.substrates.len()))
        } else {
            None
        }
    }

    /// Checks wether this `GenomicReceptor` contains any reference to the specified [`Substrate`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn referes_to_substrate(&self, substrate_index: usize) -> bool {
        self.enzyme.referes_to_substrate(substrate_index) || self.triggers.iter().chain(self.substrates.iter()).any(|element| element == &substrate_index)
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
            *trigger = Gene::adjust_index(*trigger, removed_index);
        }
        for substrate in &mut self.substrates {
            *substrate = Gene::adjust_index(*substrate, removed_index);
        }
        self.enzyme.adjust_indices(removed_index);
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
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct GenomicCatalyticCentre {
    educts: Vec<usize>,
    products: Vec<usize>,
    reaction: Reaction,
}

impl GenomicCatalyticCentre {
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
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    pub fn new(educts: Vec<usize>, products: Vec<usize>, reaction: Reaction) -> Self {
        assert_eq!(educts.len(), reaction.get_educt_number(),
            "The number of required educts for reaction {:?} is {}, but {} educts were supplied.",
            reaction, reaction.get_educt_number(), educts.len());
        assert_eq!(products.len(), reaction.get_product_number(),
            "The number of required products for reaction {:?} is {}, but {} products were supplied.",
            reaction, reaction.get_product_number(), products.len());
        GenomicCatalyticCentre{educts, products, reaction}
    }

    /// Checks wether this `GenomicCatalyticCentre` contains any reference to the specified [`Substrate`].
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn referes_to_substrate(&self, substrate_index: usize) -> bool {
        self.educts.iter().chain(self.products.iter()).any(|element| element == &substrate_index)
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
            *educt = Gene::adjust_index(*educt, removed_index);
        }
        for product in &mut self.products {
            *product = Gene::adjust_index(*product, removed_index);
        }
    }

    /// Returns the index of a random educt [`Substrate`] of this `GenomicCatalyticCentre`
    /// if there are any.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    fn get_random_educt(&self) -> Option<usize> {
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
    fn get_random_product(&self) -> Option<usize> {
        if self.products.len() > 0 {
            Some(thread_rng().gen_range(0, self.products.len()))
        } else {
            None
        }
    }

}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
/// All mutations that can be used to randomly modify a  [`Genome`].
///
/// [`Genome`]: ./struct.Genome.html
pub enum GenomeMutation {
    /// Random association of an input.
    InputAssociation,
    /// Random removal of an input.
    InputDissociation,
    /// Random association of an output.
    OutputAssociation,
    /// Random removal an output.
    OutputDissociation,
    /// Random fusion of two genes.
    GeneFusion,
    /// Random removal of a gene.
    GeneDeletion,
    /// Random duplication of a gene.
    GeneDuplication,
    /// Random addition of an genome level substrate association.
    AssociationInsertion,
    /// Random removal of an genome level substrate association.
    AssociationDeletion,
    /// Random mutation of the value of a random genome level substrate association.
    AssociationMutationSubstrate,
    /// Random mutation of a single bit of a random genome level substrate association.
    AssociationMutationSubstrateMutationFlip,
    /// Random addition of a single bit of a random genome level substrate association.
    AssociationMutationSubstrateMutationInsertion,
    /// Random removal of a single bit of a random genome level substrate association.
    AssociationMutationSubstrateMutationDeletion,
    /// Random addition of a gene's substrate to be referenced by a genome level substrate association.
    AssociationMutationGeneInsertion,
    /// Random removal of a gene's substrate to be referenced by a genome level substrate association.
    AssociationMutationGeneDeletion,
    /// Random addition of a substrate to a gene.
    GeneMutationSubstrateInsertion,
    /// Random removal of a substrate to a gene.
    GeneMutationSubstrateDeletion,
    /// Random alteration of a substrate of a gene.
    GeneMutationSubstrateMutation,
    /// Random flipping of a bit of a substrate of a gene.
    GeneMutationSubstrateMutationFlip,
    /// Random addition of a bit of a substrate of a gene.
    GeneMutationSubstrateMutationInsertion,
    /// Random removal of a bit of a substrate of a gene.
    GeneMutationSubstrateMutationDeletion,
    /// Random addition of a receptor to a gene.
    GeneMutationReceptorInsertion,
    /// Random removal of a receptor from a gene.
    GeneMutationReceptorDeletion,
    /// Random addition of a receptor trigger to a gene.
    GeneMutationReceptorMutationTriggerInsertion,
    /// Random removal of a receptor trigger from a gene.
    GeneMutationReceptorMutationTriggerDeletion,
    /// Random mutation of a receptor state of a gene.
    GeneMutationReceptorMutationStateMutation,
    /// Random mutation of a receptor's substrates of a gene.
    GeneMutationReceptorMutationSubstratesMutation,
    /// Random mutation of a receptor's enzyme of a gene.
    GeneMutationReceptorMutationEnzymeMutation,
    /// Random mutation of a catalytic centres's educts.
    GeneMutationCatalyticCentreMutationEductMutation,
    /// Random mutation of a catalytic centres's products.
    GeneMutationCatalyticCentreMutationProductMutation,
    // TODO: Lateral gene transfer
}

impl GenomeMutation {
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
    pub fn mutate(&self, genome: &Genome) -> Option<Genome> {
        match self {
            GenomeMutation::InputAssociation => GenomeMutation::mutate_input_association(genome),
            GenomeMutation::InputDissociation => GenomeMutation::mutate_input_dissociation(genome),
            GenomeMutation::OutputAssociation => GenomeMutation::mutate_output_association(genome),
            GenomeMutation::OutputDissociation => GenomeMutation::mutate_output_dissociation(genome),
            GenomeMutation::GeneFusion => GenomeMutation::mutate_gene_fusion(genome),
            GenomeMutation::GeneDeletion => GenomeMutation::mutate_gene_deletion(genome),
            GenomeMutation::GeneDuplication => GenomeMutation::mutate_gene_duplication(genome),
            GenomeMutation::AssociationInsertion => GenomeMutation::mutate_association_insertion(genome),
            GenomeMutation::AssociationDeletion => GenomeMutation::mutate_association_deletion(genome),
            GenomeMutation::AssociationMutationSubstrate => GenomeMutation::mutate_association_substrate(genome),
            GenomeMutation::AssociationMutationSubstrateMutationFlip => GenomeMutation::mutate_association_substrate_mutation_flip(genome),
            GenomeMutation::AssociationMutationSubstrateMutationInsertion => GenomeMutation::mutate_association_substrate_mutation_insertion(genome),
            GenomeMutation::AssociationMutationSubstrateMutationDeletion => GenomeMutation::mutate_association_substrate_mutation_deletion(genome),
            GenomeMutation::AssociationMutationGeneInsertion => GenomeMutation::mutate_association_gene_insertion(genome),
            GenomeMutation::AssociationMutationGeneDeletion => GenomeMutation::mutate_association_gene_deletion(genome),
            GenomeMutation::GeneMutationSubstrateInsertion => GenomeMutation::mutate_gene_substrate_insertion(genome),
            GenomeMutation::GeneMutationSubstrateDeletion => GenomeMutation::mutate_gene_substrate_deletion(genome),
            GenomeMutation::GeneMutationSubstrateMutation => GenomeMutation::mutate_gene_substrate_mutation(genome),
            GenomeMutation::GeneMutationSubstrateMutationFlip => GenomeMutation::mutate_gene_substrate_mutation_flip(genome),
            GenomeMutation::GeneMutationSubstrateMutationInsertion => GenomeMutation::mutate_gene_substrate_mutation_insertion(genome),
            GenomeMutation::GeneMutationSubstrateMutationDeletion => GenomeMutation::mutate_gene_substrate_mutation_deletion(genome),
            GenomeMutation::GeneMutationReceptorInsertion => GenomeMutation::mutate_gene_receptor_insertion(genome),
            GenomeMutation::GeneMutationReceptorDeletion => GenomeMutation::mutate_gene_receptor_deletion(genome),
            GenomeMutation::GeneMutationReceptorMutationTriggerInsertion => GenomeMutation::mutate_gene_receptor_trigger_insertion(genome),
            GenomeMutation::GeneMutationReceptorMutationTriggerDeletion => GenomeMutation::mutate_gene_receptor_trigger_deletion(genome),
            GenomeMutation::GeneMutationReceptorMutationStateMutation => GenomeMutation::mutate_gene_receptor_state_mutation(genome),
            GenomeMutation::GeneMutationReceptorMutationSubstratesMutation => GenomeMutation::mutate_gene_receptor_substrate_mutation(genome),
            GenomeMutation::GeneMutationReceptorMutationEnzymeMutation => GenomeMutation::mutate_gene_receptor_enzyme_mutation(genome),
            GenomeMutation::GeneMutationCatalyticCentreMutationEductMutation => GenomeMutation::mutate_gene_catalytic_centre_educt_mutation(genome),
            GenomeMutation::GeneMutationCatalyticCentreMutationProductMutation => GenomeMutation::mutate_gene_catalytic_centre_product_mutation(genome),
            // GenomeMutation::LateralGeneTransfer => None, //TODO: implement a global gene pool
        }
    }

    /// Duplicates the [`Genome`] and then duplicates a random [`Gene`] inside of the [`Genome`].
    /// Returns the altered [`Genome`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Genome`]: ./struct.Genome.html
    fn mutate_gene_duplication(genome: &Genome) -> Option<Genome> {
        let mut mutated_genome = genome.duplicate();
        mutated_genome.duplicate_gene_internal(mutated_genome.get_random_gene());
        Some(mutated_genome)
    }

    /// Duplicates the [`Genome`] and then associates a random input [`GeneSubstrate`] of the [`Genome`]
    /// with a random [`Substrate`] of a random [`Gene`].
    /// Returns the altered [`Genome`] if there are any input [`GeneSubstrate`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_input_association(genome: &Genome) -> Option<Genome> {
       if genome.number_of_inputs() > 0 {
           let mut mutated_genome = genome.duplicate();
           // It was checked before if there are any inputs, so the unwrap must work.
           let random_input = mutated_genome.get_random_input().unwrap();
           let random_gene = mutated_genome.get_random_gene();
           let random_gene_substrate = mutated_genome.get_gene(random_gene).get_random_substrate();
           mutated_genome.set_input(random_input, Some(GeneSubstrate{
               gene: random_gene,
               substrate: random_gene_substrate
           }));
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and then dissociates a random input [`GeneSubstrate`] of the [`Genome`].
   /// Returns the altered [`Genome`] if there are any input [`GeneSubstrate`] and the randomly selected
   /// input [`GeneSubstrate`] was previously associated.
   ///
   /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_input_dissociation(genome: &Genome) -> Option<Genome> {
       if genome.number_of_inputs() > 0 {
           let mut mutated_genome = genome.duplicate();
           // It was checked before if there are any inputs, so the unwrap must work.
           let random_input = mutated_genome.get_random_input().unwrap();
           if let None = mutated_genome.input[random_input] {
               // If there is no association, return no change.
               None
           } else {
               // If there is an association, remove it.
               mutated_genome.set_input(random_input, None);
               Some(mutated_genome)
           }
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and then associates a random output [`GeneSubstrate`] of the [`Genome`]
   /// with a random [`Substrate`] of a random [`Gene`].
   /// Returns the altered [`Genome`] if there are any output [`GeneSubstrate`].
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_output_association(genome: &Genome) -> Option<Genome> {
       if genome.number_of_outputs() > 0 {
           let mut mutated_genome = genome.duplicate();
           // It was checked before if there are any outputs, so the unwrap must work.
           let random_output = mutated_genome.get_random_output().unwrap();
           let random_gene = mutated_genome.get_random_gene();
           let random_gene_substrate = mutated_genome.get_gene(random_gene).get_random_substrate();
           mutated_genome.set_output(random_output, Some(GeneSubstrate{
               gene: random_gene,
               substrate: random_gene_substrate
           }));
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and then dissociates a random output [`GeneSubstrate`] of the [`Genome`].
   /// Returns the altered [`Genome`] if there are any output [`GeneSubstrate`] and the randomly selected
   /// output [`GeneSubstrate`] was previously associated.
   ///
   /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_output_dissociation(genome: &Genome) -> Option<Genome> {
       if genome.number_of_outputs() > 0 {
           let mut mutated_genome = genome.duplicate();
           // It was checked before if there are any outputs, so the unwrap must work.
           let random_output = mutated_genome.get_random_output().unwrap();
           if let None = mutated_genome.input[random_output] {
               // If there is no association, return no change.
               None
           } else {
               // If there is an association, remove it.
               mutated_genome.set_input(random_output, None);
               Some(mutated_genome)
           }
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and then fuses two random [`Gene`]s of the [`Genome`].
   /// The fusion product is added as a new [`Gene`].
   /// Returns the altered [`Genome`] if there are 2 or more different [`Gene`]s in the [`Genome`]
   /// and if the fusion process was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_gene_fusion(genome: &Genome) -> Option<Genome> {
       // Only allow this mutation if there is more than one gene in the genome, since
       // 2 different genes are needed for fusion.
       if genome.number_of_genes().get() > 1 {
           let mut mutated_genome = genome.duplicate();
           // Select two different genes from the genome.
           let mut genes: Vec<usize> = (0..mutated_genome.number_of_genes().get()).collect();
           let gene_a = genes.remove(thread_rng().gen_range(0, genes.len()));
           let gene_b = genes.remove(thread_rng().gen_range(0, genes.len()));
           if let Some(fusion_gene) = mutated_genome.duplicate_gene(gene_a).fuse(mutated_genome.get_gene(gene_b)) {
               if let Some(_) = mutated_genome.add_gene(fusion_gene) {
                   return Some(mutated_genome)
               }
           }
       }
       // Return nothing if any of the prequisites is not met.
       None
   }

   /// Duplicates the [`Genome`] and then removes a random [`Gene`] from it.
   /// Returns the altered [`Genome`] if there was more than 1 [`Gene`] in the [`Genome`].
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_gene_deletion(genome: &Genome) -> Option<Genome> {
       // Only allow this mutation if there is more than one gene in the genome,
       // otherwise the genome would end up in an invalid state.
       if genome.number_of_genes().get() > 1 {
           let mut mutated_genome = genome.duplicate();
           let random_gene = mutated_genome.get_random_gene();
           mutated_genome.remove_gene(random_gene);
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and adds a random [`GeneAssociation`] to it.
   /// Returns the altered [`Genome`] if the [`GeneAssociation`] could be established.
   ///
   /// [`Gene`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_association_insertion(genome: &Genome) -> Option<Genome> {
       let mut mutated_genome = genome.duplicate();
       let random_association = GeneAssociation {
           substrate: random_substrate(),
           associations: vec!(mutated_genome.random_gene_substrate())
       };
       mutated_genome.add_association(random_association).and(Some(mutated_genome))
   }

   /// Duplicates the [`Genome`] and removes a random [`GeneAssociation`] from it.
   /// Returns the altered [`Genome`] if 1 or more [`GeneAssociation`]s were present.
   ///
   /// [`Gene`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_association_deletion(genome: &Genome) -> Option<Genome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           mutated_genome.remove_association(random_association_index);
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly sets an existing [`GeneAssociation`]'s [`Substrate`] value.
   /// Returns the altered [`Genome`] if 1 or more [`GeneAssociation`]s were present.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_association_substrate(genome: &Genome) -> Option<Genome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           let random_substrate = mutate_substrate_based_on(mutated_genome.associations[random_association_index].substrate.len() / 8);
           // If the substrates are the same, report no change.
           if random_substrate != mutated_genome.associations[random_association_index].substrate {
               mutated_genome.associations[random_association_index].substrate = random_substrate;
               Some(mutated_genome)
           } else {
               None
           }
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly modifies a bit of a [`Substrate`] of a [`GeneAssociation`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_association_substrate_mutation_flip(genome: &Genome) -> Option<Genome> {
       genome.get_random_association().and_then(|random_association_index|{
           if genome.associations[random_association_index].substrate.len() > 0 {
               let mut mutated_genome = genome.duplicate();
               mutate_substrate_single_bit(&mut mutated_genome.associations[random_association_index].substrate);
               Some(mutated_genome)
           } else {
               None
           }
       })
   }

   /// Duplicates the [`Genome`] and randomly adds a bit to a [`Substrate`] of a [`GeneAssociation`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_association_substrate_mutation_insertion(genome: &Genome) -> Option<Genome> {
       genome.get_random_association().and_then(|random_association_index|{
               let mut mutated_genome = genome.duplicate();
               let mut mutated_substrate = mutated_genome.associations[random_association_index].substrate.clone();
               mutated_substrate = mutate_substrate_single_bit_insertion(mutated_substrate);
               mutated_genome.associations[random_association_index].substrate = mutated_substrate;
               Some(mutated_genome)
       })
   }

   /// Duplicates the [`Genome`] and randomly removes a bit to a [`Substrate`] of a [`GeneAssociation`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_association_substrate_mutation_deletion(genome: &Genome) -> Option<Genome> {
       genome.get_random_association().and_then(|random_association_index|{
           if genome.associations[random_association_index].substrate.len() > 0 {
               let mut mutated_genome = genome.duplicate();
               let mut mutated_substrate = mutated_genome.associations[random_association_index].substrate.clone();
               mutated_substrate = mutate_substrate_single_bit_deletion(mutated_substrate);
               mutated_genome.associations[random_association_index].substrate = mutated_substrate;
               Some(mutated_genome)
           } else {
               None
           }
       })
   }

   /// Duplicates the [`Genome`] and randomly adds a [`GeneSubstrate`] to an existing [`GeneAssociation`].
   /// Returns the altered [`Genome`] if 1 or more [`GeneAssociation`]s were present.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GeneAssociation`]: ./struct.GeneSubstrate.html
   fn mutate_association_gene_insertion(genome: &Genome) -> Option<Genome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           let random_gene_substrate = mutated_genome.random_gene_substrate();
            mutated_genome.associations[random_association_index].add_association(random_gene_substrate).and(Some(mutated_genome))
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly removes a [`GeneSubstrate`] from an existing [`GeneAssociation`].
   /// Returns the altered [`Genome`] if 1 or more [`GeneAssociation`]s were present and a [`GeneSubstrate`]
   /// could be removed.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GeneAssociation`]: ./struct.GeneSubstrate.html
   fn mutate_association_gene_deletion(genome: &Genome) -> Option<Genome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           if mutated_genome.associations[random_association_index].associations.len() > 0 {
               if let Some(random_deletion_index) = mutated_genome.associations[random_association_index].get_random_gene_substrate() {
                   mutated_genome.associations[random_association_index].associations.remove(random_deletion_index);
                   Some(mutated_genome)
               } else {
                   None
               }
           } else {
               None
           }
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly adds a [`Substrate`] to an existing [`Gene`].
   /// Returns the altered [`Genome`] if the addition was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_insertion(genome: &Genome) -> Option<Genome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       mutated_genome.genes[random_gene_index].add_substrate(random_substrate()).and(Some(mutated_genome))
   }

   /// Duplicates the [`Genome`] and randomly removes a [`Substrate`] from an existing [`Gene`].
   /// Returns the altered [`Genome`] if more than 1 [`Substrate`] was present for the selected [`Gene`].
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_deletion(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       if genome.get_gene(random_gene_index).number_of_substrates().get() > 1 {
           let mut mutated_genome = genome.duplicate();
           let removed_substrate = GeneSubstrate{
               gene: random_gene_index,
               substrate: mutated_genome.genes[random_gene_index].get_random_substrate(),
           };
           mutated_genome.genes[removed_substrate.gene].remove_substrate(removed_substrate.substrate);
           mutated_genome.adjust_input_after_gene_substrate_removal(removed_substrate);
           mutated_genome.adjust_output_after_gene_substrate_removal(removed_substrate);
           mutated_genome.adjust_associations_after_gene_substrate_removal(removed_substrate);
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly modifies a [`Substrate`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_mutation(genome: &Genome) -> Option<Genome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       let random_substrate_index = mutated_genome.get_gene(random_gene_index).get_random_substrate();
       let random_substrate_value = mutate_substrate_based_on(mutated_genome.get_gene(random_gene_index).substrates[random_substrate_index].len() / 8);
       if random_substrate_value != mutated_genome.get_gene(random_gene_index).substrates[random_substrate_index] {
           mutated_genome.genes[random_gene_index].substrates[random_substrate_index] = random_substrate_value;
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly modifies a bit of a [`Substrate`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_mutation_flip(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       let random_substrate_index = genome.get_gene(random_gene_index).get_random_substrate();
       if genome.get_gene(random_gene_index).substrates[random_substrate_index].len() > 0 {
           let mut mutated_genome = genome.duplicate();
           mutate_substrate_single_bit(&mut mutated_genome.genes[random_gene_index].substrates[random_substrate_index]);
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly adds a bit to a [`Substrate`] of an existing [`Gene`].
   /// Returns the altered [`Genome`].
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_mutation_insertion(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       let random_substrate_index = genome.get_gene(random_gene_index).get_random_substrate();
           let mut mutated_genome = genome.duplicate();
           let mut mutated_substrate = mutated_genome.genes[random_gene_index].substrates[random_substrate_index].clone();
           mutated_substrate = mutate_substrate_single_bit_insertion(mutated_substrate);
           mutated_genome.genes[random_gene_index].substrates[random_substrate_index] = mutated_substrate;
           Some(mutated_genome)
   }

   /// Duplicates the [`Genome`] and randomly removes a bit of a [`Substrate`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_mutation_deletion(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       let random_substrate_index = genome.get_gene(random_gene_index).get_random_substrate();
       if genome.get_gene(random_gene_index).substrates[random_substrate_index].len() > 0 {
           let mut mutated_genome = genome.duplicate();
           let mut mutated_substrate = mutated_genome.genes[random_gene_index].substrates[random_substrate_index].clone();
           mutated_substrate = mutate_substrate_single_bit_deletion(mutated_substrate);
           mutated_genome.genes[random_gene_index].substrates[random_substrate_index] = mutated_substrate;
           Some(mutated_genome)
       } else {
           None
       }
   }

   /// Duplicates the [`Genome`] and randomly adds a [`GenomicReceptor`] to an existing [`Gene`].
   /// Returns the altered [`Genome`] if the addition was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   fn mutate_gene_receptor_insertion(genome: &Genome) -> Option<Genome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       let random_receptor = mutated_genome.get_gene(random_gene_index).random_receptor();
       mutated_genome.genes[random_gene_index].add_receptor(random_receptor).and(Some(mutated_genome))
   }

   /// Duplicates the [`Genome`] and randomly removes a [`GenomicReceptor`] from an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   fn mutate_gene_receptor_deletion(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       let mut mutated_genome = genome.duplicate();
        mutated_genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            mutated_genome.genes[random_gene_index].receptors.remove(random_receptor);
            Some(mutated_genome)
       })
   }

   /// Duplicates the [`Genome`] and randomly adds a trigger to a [`GenomicReceptor`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if the addition was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   fn mutate_gene_receptor_trigger_insertion(genome: &Genome) -> Option<Genome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       mutated_genome.genes[random_gene_index].get_random_receptor().and_then(|random_receptor| {
           let random_substrate_index = mutated_genome.get_gene(random_gene_index).get_random_substrate();
           mutated_genome.genes[random_gene_index].receptors[random_receptor].add_trigger(random_substrate_index).and(Some(mutated_genome))
       })
   }

   /// Duplicates the [`Genome`] and randomly removes a trigger from a [`GenomicReceptor`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] and triggers were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   fn mutate_gene_receptor_trigger_deletion(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
           genome.get_gene(random_gene_index).receptors[random_receptor].get_random_trigger().and_then(|random_trigger| {
               let mut mutated_genome = genome.duplicate();
               mutated_genome.genes[random_gene_index].receptors[random_receptor].triggers.remove(random_trigger);
               Some(mutated_genome)
           })
       })
   }

   /// Duplicates the [`Genome`] and randomly mutates a [`State`] of a [`GenomicReceptor`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   /// [`State`]: ../chemistry/struct.State.html
   fn mutate_gene_receptor_state_mutation(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
           let mut mutated_genome = genome.duplicate();
           let state: State = rand::random();
           let substrates = (0..state.get_substrate_number()).map(|_| mutated_genome.get_gene(random_gene_index).get_random_substrate()).collect();
           mutated_genome.genes[random_gene_index].receptors[random_receptor].state = state;
           mutated_genome.genes[random_gene_index].receptors[random_receptor].substrates = substrates;
           Some(mutated_genome)
       })
   }

   /// Duplicates the [`Genome`] and randomly mutates a [`Substrate`] of a [`GenomicReceptor`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] and [`Substrate`]s were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_receptor_substrate_mutation(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
           genome.get_gene(random_gene_index).receptors[random_receptor].get_random_substrate().and_then(|random_substrate| {
               let mut mutated_genome = genome.duplicate();
               mutated_genome.genes[random_gene_index].receptors[random_receptor].substrates[random_substrate] = mutated_genome.get_gene(random_gene_index).get_random_substrate();
               Some(mutated_genome)
           })
       })
   }

   /// Duplicates the [`Genome`] and randomly mutates a [`Reaction`] of a [`GenomicReceptor`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`]s were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   /// [`Reaction`]: ../chemistry/struct.Reaction.html
   fn mutate_gene_receptor_enzyme_mutation(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
           let mut mutated_genome = genome.duplicate();
           let enzyme = mutated_genome.get_gene(random_gene_index).random_catalytic_centre();
           mutated_genome.genes[random_gene_index].receptors[random_receptor].enzyme = enzyme;
           Some(mutated_genome)
       })
   }

   /// Duplicates the [`Genome`] and randomly mutates an educt [`Substrate`] of a [`GenomicCatalyticCentre`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] and educt [`Substrate`]s were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   /// [`GenomicCatalyticCentre`]: ./struct.GenomicCatalyticCentre.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_catalytic_centre_educt_mutation(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
           genome.get_gene(random_gene_index).receptors[random_receptor].enzyme.get_random_educt().and_then(|random_educt| {
               let mut mutated_genome = genome.duplicate();
               mutated_genome.genes[random_gene_index].receptors[random_receptor].enzyme.educts[random_educt] = mutated_genome.get_gene(random_gene_index).get_random_substrate();
               Some(mutated_genome)
           })
       })
   }

   /// Duplicates the [`Genome`] and randomly mutates an product [`Substrate`] of a [`GenomicCatalyticCentre`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] and product [`Substrate`]s were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   /// [`GenomicCatalyticCentre`]: ./struct.GenomicCatalyticCentre.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_catalytic_centre_product_mutation(genome: &Genome) -> Option<Genome> {
       let random_gene_index = genome.get_random_gene();
       genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
           genome.get_gene(random_gene_index).receptors[random_receptor].enzyme.get_random_product().and_then(|random_product| {
               let mut mutated_genome = genome.duplicate();
               mutated_genome.genes[random_gene_index].receptors[random_receptor].enzyme.products[random_product] = mutated_genome.get_gene(random_gene_index).get_random_substrate();
               Some(mutated_genome)
           })
       })
   }

}

/// Generates a random binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn random_substrate() -> BitBox<Local, u8> {
    let length: usize = thread_rng().gen_range(RANDOM_SUBSTRATE_MIN_LENGTH, RANDOM_SUBSTRATE_MAX_LENGTH + 1);
    let random_bytes: Vec<u8> = (0..length).map(|_| thread_rng().gen::<u8>()).collect();
    BitBox::from_slice(&random_bytes)
}

/// Generates a random binary [`Substrate`] based on the specified length in byte.
/// The resulting [`Substrate`] will be of length 1 B if the base length was 0.
/// Otherwise the resulting length will be in range `base length / 2` to `base_length * 2`.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn mutate_substrate_based_on(base_length: usize) -> BitBox<Local, u8> {
    let length: usize;
    if base_length == 0 {
        // If the substrate was empty, create the smallest possible substrate to begin with.
        length = 1;
    } else {
        let min_length: usize = base_length / 2;
        let max_length: usize = if let Some(u) = base_length.checked_mul(2) {
            u
        } else {
            usize::max_value()
        };
        // This limits the generated number to `usize::max_value() - 1`, which should be enough
        // and not pose a problem.
        length = thread_rng().gen_range(min_length, max_length);
    }
    let random_bytes: Vec<u8> = (0..length).map(|_| thread_rng().gen::<u8>()).collect();
    BitBox::from_slice(&random_bytes)
}

/// Flips a random bit of the binary [`Substrate`].
///
/// # Parameters
///
/// * `base_substrate` - the [`Substrate`] to mutate
///
/// Panics
///
/// If `base_substrate` is empty.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn mutate_substrate_single_bit(base_substrate: &mut BitBox<Local, u8>) {
    if base_substrate.is_empty() {
        panic!("No byte can be flipped in an empty substrate.");
    }
    let random_bit_index = thread_rng().gen_range(0, base_substrate.len());
    // The unwrap should always work, since the index being in range was checked for.
    // The clone should be cheap as a bool primitive is cloned.
    let random_bit = base_substrate.get(random_bit_index).unwrap().clone();
    // Flip the bit.
    base_substrate.set(random_bit_index, random_bit ^ true);
}

/// Adds a random bit to the binary [`Substrate`].
///
/// # Parameters
///
/// * `base_substrate` - the [`Substrate`] to mutate
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn mutate_substrate_single_bit_insertion(base_substrate: BitBox<Local, u8>) -> BitBox<Local, u8> {
    let mut mutated_substrate = BitVec::from_boxed_bitslice(base_substrate);
    if mutated_substrate.is_empty() {
        mutated_substrate.push(thread_rng().gen());
    } else {
        let random_bit_index = thread_rng().gen_range(0, mutated_substrate.len());
        mutated_substrate.insert(random_bit_index, thread_rng().gen());
    }
    mutated_substrate.into_boxed_bitslice()
}

/// Removes a random bit from the binary [`Substrate`].
///
/// # Parameters
///
/// * `base_substrate` - the [`Substrate`] to mutate
///
/// Panics
///
/// If `base_substrate` is empty.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn mutate_substrate_single_bit_deletion(base_substrate: BitBox<Local, u8>) -> BitBox<Local, u8> {
    if base_substrate.is_empty() {
        panic!("No byte can be removed from an empty substrate.");
    }
    let mut mutated_substrate = BitVec::from_boxed_bitslice(base_substrate);
    let random_bit_index = thread_rng().gen_range(0, mutated_substrate.len());
    mutated_substrate.remove(random_bit_index);
    mutated_substrate.into_boxed_bitslice()
}

impl Distribution<GenomeMutation> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> GenomeMutation {
        match rng.gen_range(0u8, 30) {
            0 => GenomeMutation::InputAssociation,
            1 => GenomeMutation::InputDissociation,
            2 => GenomeMutation::OutputAssociation,
            3 => GenomeMutation::OutputDissociation,
            4 => GenomeMutation::GeneFusion,
            5 => GenomeMutation::GeneDeletion,
            6 => GenomeMutation::GeneDuplication,
            7 => GenomeMutation::AssociationInsertion,
            8 => GenomeMutation::AssociationDeletion,
            9 => GenomeMutation::AssociationMutationSubstrate,
            10 => GenomeMutation::AssociationMutationGeneInsertion,
            11 => GenomeMutation::AssociationMutationGeneDeletion,
            12 => GenomeMutation::GeneMutationSubstrateInsertion,
            13 => GenomeMutation::GeneMutationSubstrateDeletion,
            14 => GenomeMutation::GeneMutationSubstrateMutation,
            15 => GenomeMutation::GeneMutationReceptorInsertion,
            16 => GenomeMutation::GeneMutationReceptorDeletion,
            17 => GenomeMutation::GeneMutationReceptorMutationTriggerInsertion,
            18 => GenomeMutation::GeneMutationReceptorMutationTriggerDeletion,
            19 => GenomeMutation::GeneMutationReceptorMutationStateMutation,
            20 => GenomeMutation::GeneMutationReceptorMutationSubstratesMutation,
            21 => GenomeMutation::GeneMutationReceptorMutationEnzymeMutation,
            22 => GenomeMutation::GeneMutationCatalyticCentreMutationEductMutation,
            23 => GenomeMutation::GeneMutationCatalyticCentreMutationProductMutation,
            24 => GenomeMutation::GeneMutationSubstrateMutationFlip,
            25 => GenomeMutation::GeneMutationSubstrateMutationInsertion,
            26 => GenomeMutation::GeneMutationSubstrateMutationDeletion,
            27 => GenomeMutation::AssociationMutationSubstrateMutationFlip,
            28 => GenomeMutation::AssociationMutationSubstrateMutationInsertion,
            29 => GenomeMutation::AssociationMutationSubstrateMutationDeletion,
            _ => panic!("A random number with no matching genomic mutation was created."),
        }
    }
}
