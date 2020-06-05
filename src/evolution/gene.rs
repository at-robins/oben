//! The `gene` module contains information structures to describe,
//! replicate and mutate the evolutionary network.
extern crate bitvec;
extern crate rand;

use super::chemistry::{Reaction, State};
use bitvec::{boxed::BitBox, order::Local};
use rand::{distributions::{Distribution, Standard}, thread_rng, Rng};
use std::num::NonZeroUsize;

/// The minimal length in byte of a randomly created binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
const RANDOM_SUBSTRATE_MIN_LENGTH: usize = 0;
/// The maximal length in byte of a randomly created binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
const RANDOM_SUBSTRATE_MAX_LENGTH: usize = 128;

/// A `Genome` is a collection of individual [`Gene`]s and associations between them.
/// A `Genome` is required to consist of 1 or more genes.
///
/// [`Gene`]: ./struct.Gene.html
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
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
    /// [`Gene`]: ./struct.GeneAssociation.html
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
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
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

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
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
        // TODO: Replace with drain_filter once stabilised.
        let mut i = 0;
        while i != self.associations.len() {
            if self.associations[i].gene == gene {
                self.associations.remove(i);
            } else {
                i += 1;
            }
        }
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
}

/// A `Gene` is an immutable structure encoding a self-contained network, but without
/// explicite function. It can be transcribed into a functional protein network.
/// A `Gene` is required to encode at least 1 [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
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


    /// Returns the index of a random [`Substrate`] encoded by this `Gene`.
    ///
    /// [`Substrate`]: ../protein/struct.Substrate.html
    pub fn get_random_substrate(&self) -> usize {
        thread_rng().gen_range(0, self.number_of_substrates().get())
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
}

/// A `GenomicReceptor` represents the information of an actual [`Receptor`] that
/// is triggered upon specific [`Substrate`] changes and compares [`Substrate`], eventually
/// performing a [`Reaction`]. It is contained within a [`Gene`].
///
/// [`Receptor`]: ../protein/struct.Receptor.html
/// [`Gene`]: ./struct.Gene.html
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
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
}

/// A `GenomicCatalyticCentre` represents the information of an actual [`CatalyticCentre`] that
/// produces products from educt [`Substrate`]s
/// by performing a [`Reaction`]. It is contained within a [`Gene`].
///
/// [`CatalyticCentre`]: ../protein/struct.CatalyticCentre.html
/// [`Gene`]: ./struct.Gene.html
/// [`Substrate`]: ../protein/struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
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
}

pub struct MutagenicEnvironment {
    // mutation rate per size in byte
    mutation_rate: f64,
    // chance for mutation on the genome level, otherwise a mutation specific to a gene
    genomic_mutation_rate: f64,

}

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
    // TODO: Lateral gene transfer
    // TODO: Gene mutations
    // TODO: Genomic receptor
    // TODO: Genomic catalytic centre
}

impl GenomeMutation {
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
            GenomeMutation::AssociationMutationGeneInsertion => GenomeMutation::mutate_association_gene_insertion(genome),
            GenomeMutation::AssociationMutationGeneDeletion => GenomeMutation::mutate_association_gene_deletion(genome),
            GenomeMutation::GeneMutationSubstrateInsertion => GenomeMutation::mutate_gene_substrate_insertion(genome),
            GenomeMutation::GeneMutationSubstrateDeletion => GenomeMutation::mutate_gene_substrate_deletion(genome),
            GenomeMutation::GeneMutationSubstrateMutation => GenomeMutation::mutate_gene_substrate_mutation(genome),
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
           let random_substrate_index = mutated_genome.genes[random_gene_index].get_random_substrate();
           mutated_genome.genes[random_gene_index].substrates.remove(random_substrate_index);
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

impl Distribution<GenomeMutation> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> GenomeMutation {
        match rng.gen_range(0u8, 15) {
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
            _ => panic!("A random number with no matching genomic mutation was created.")
        }
    }
}
