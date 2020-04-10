//! The `gene` module contains information structures to describe,
//! replicate and mutate the evolutionary network.
extern crate bitvec;
extern crate rand;

use super::chemistry::{Reaction, State};
use bitvec::{boxed::BitBox};
use rand::{thread_rng, Rng};
use std::num::NonZeroUsize;

/// A `Genome` is a collection of individual [`Gene`]s and associations between them.
/// A `Genome` is required to consist of 1 or more genes.
///
/// [`Gene`]: ./struct.Gene.html
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct Genome {
    input: Vec<GeneSubstrate>,
    output: Vec<GeneSubstrate>,
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

    /// Returns the index of a random [`Gene`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn get_random_gene(&self) -> usize {
        thread_rng().gen_range(0, self.number_of_genes().get())
    }

    /// Duplicates the `Genome` and all its contents.
    pub fn duplicate(&self) -> Self {
        // At the moment this is just a wrapper for cloning.
        self.clone()
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
    substrate: BitBox,
    // gene specific substrates pointing to the shared substrate
    associations: Vec<GeneSubstrate>,
}

/// A `Gene` is an immutable structure encoding a self-contained network, but without
/// explicite function. It defines interfaces with other genes by the means of input/output
/// substrates. It can be transcribed into a functional protein network.
/// A `Gene` is required to encode at least 1 [`Substrate`]s.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct Gene {
    substrates: Vec<BitBox>,
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

    pub fn fuse(&self, other_gene: &Gene) -> Option<Gene> {
        if let None = self.substrates.len().checked_add(other_gene.substrates.len()) {
            return None;
        };
        let mut fusion_gene = self.duplicate();
        fusion_gene.substrates.append(&mut other_gene.substrates.clone());
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
/// [`Substrate`]: ./struct.Substrate.html
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

enum GenomeMutation {
    ///
    InputAssociation,
    OutputAssociation,
    SubstrateInsertion,
    SubstrateDeletion,
    SubstrateMutation,
    GeneFusion,
    GeneDeletion,
    GeneDuplication,
    LateralGeneTransfer,
    AssociationInsertion,
    AssociationDeletion,
    AssociationMuation,
}

impl GenomeMutation {
    fn mutate(&self, genome: &Genome) -> Option<Genome> {
        let mut mutated_genome = genome.duplicate();
        match self {
            GenomeMutation::InputAssociation if mutated_genome.number_of_inputs() > 0 => {
                let random_input = thread_rng().gen_range(0, mutated_genome.number_of_inputs());
                let random_gene = mutated_genome.get_random_gene();
                let random_gene_substrate = thread_rng().gen_range(0, mutated_genome.genes[random_gene].number_of_substrates().get());
                mutated_genome.input[random_input] = GeneSubstrate{
                    gene: random_gene,
                    substrate: random_gene_substrate
                };
            },
            GenomeMutation::OutputAssociation if mutated_genome.number_of_outputs() > 0 => {
                let random_output = thread_rng().gen_range(0, mutated_genome.number_of_outputs());
                let random_gene = mutated_genome.get_random_gene();
                let random_gene_substrate = thread_rng().gen_range(0, mutated_genome.genes[random_gene].number_of_substrates().get());
                mutated_genome.output[random_output] = GeneSubstrate{
                    gene: random_gene,
                    substrate: random_gene_substrate
                };
            },
            GenomeMutation::SubstrateInsertion => {},
            GenomeMutation::SubstrateDeletion => {},
            GenomeMutation::SubstrateMutation => {},
            // Only allow this mutation if there is more than one gene in the genome.
            GenomeMutation::GeneFusion if mutated_genome.number_of_genes().get() > 1 => {
                // Select two different genes from the genome.
                let mut genes: Vec<usize> = (0..mutated_genome.number_of_genes().get()).collect();
                let gene_a = genes.remove(thread_rng().gen_range(0, genes.len()));
                let gene_b = genes.remove(thread_rng().gen_range(0, genes.len()));
                let fusion_gene = mutated_genome.duplicate_gene(gene_a);

                // TODO: finish
            },
            GenomeMutation::GeneDeletion => {
                let random_gene = mutated_genome.get_random_gene();
                // Delete the gene from the genome.
                mutated_genome.remove_gene(random_gene);
                // TODO: reassociate input and output
                // TODO: delete all other substrate associations
                // mutated_genome.associations = mutated_genome.associations.iter().filter(|association| association.gene == random_gene).collect();
            },
            GenomeMutation::GeneDuplication => {
                mutated_genome.duplicate_gene_internal(mutated_genome.get_random_gene());
            },
            GenomeMutation::AssociationInsertion => {},
            GenomeMutation::AssociationDeletion => {},
            GenomeMutation::AssociationMuation => {},
            GenomeMutation::LateralGeneTransfer => return None, //TODO: implement a global gene pool
            _ => return None,
        };
        Some(mutated_genome)
    }
}
