//! The `gene` module contains information structures to describe,
//! replicate and mutate the evolutionary network.
extern crate bitvec;
extern crate rand;

use super::chemistry::{Reaction, State};
use bitvec::{boxed::BitBox};

/// A `Genome` is a collection of individual [`Gene`]s and associations between them.
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
    ///
    /// [`Gene`]: ./struct.Gene.html
    pub fn number_of_genes(&self) -> usize {
        self.genes.len()
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
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct Gene {
    substrates: Vec<BitBox>,
    io: Vec<usize>, // indices of the input/output interface substrates
    receptors: Vec<GenomicReceptor>,
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
    GeneInsertion,
    GeneDeletion,
    AssociationInsertion,
    AssociationDeletion,
    AssociationMuation,
}

impl GenomeMutation {
    fn mutate(&self, genome: &Genome) -> Genome {
        let mutated_genome = genome.clone();
        match self {
            GenomeMutation::InputAssociation => {

            },
            GenomeMutation::OutputAssociation => {},
            GenomeMutation::SubstrateInsertion => {},
            GenomeMutation::SubstrateDeletion => {},
            GenomeMutation::SubstrateMutation => {},
            GenomeMutation::GeneInsertion => {},
            GenomeMutation::GeneDeletion => {},
            GenomeMutation::AssociationInsertion => {},
            GenomeMutation::AssociationDeletion => {},
            GenomeMutation::AssociationMuation => {},
        };
        mutated_genome
    }
}
