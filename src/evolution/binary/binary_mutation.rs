//! The `gene` module contains information structures to describe,
//! replicate and mutate the evolutionary network.
extern crate bitvec;
extern crate rand;

use bitvec::vec::BitVec;
use rand::{distributions::{Distribution, Standard}, thread_rng, Rng};
use super::{BinaryGenome, BinaryReaction, BinarySubstrate, BinaryState};
use super::super::gene::{Gene, GeneAssociation, GenomeMutation};
use super::super::chemistry::State;

/// The minimal length in byte of a randomly created binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
const RANDOM_SUBSTRATE_MIN_LENGTH: usize = 0;
/// The maximal length in byte of a randomly created binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
const RANDOM_SUBSTRATE_MAX_LENGTH: usize = 1;

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
/// All mutations that can be used to randomly modify a binary [`Genome`].
///
/// [`Genome`]: ./struct.Genome.html
pub enum BinaryMutation {
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

impl BinaryMutation {
    /// Duplicates the [`Genome`] and then duplicates a random [`Gene`] inside of the [`Genome`].
    /// Returns the altered [`Genome`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Genome`]: ./struct.Genome.html
    fn mutate_gene_duplication(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let mut mutated_genome = genome.duplicate();
        mutated_genome.add_gene(mutated_genome.duplicate_random_gene()).and_then(|_| Some(mutated_genome))
    }

    /// Duplicates the [`Genome`] and then associates a random input [`GeneSubstrate`] of the [`Genome`]
    /// with a random [`Substrate`] of a random [`Gene`].
    /// Returns the altered [`Genome`] if there are any input [`GeneSubstrate`].
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_input_association(genome: &BinaryGenome) -> Option<BinaryGenome> {
       if genome.number_of_inputs() > 0 {
           // It was checked before if there are any inputs, so the unwrap must work.
           let random_input = genome.get_random_input().unwrap();
           let random_gene_substrate = genome.random_gene_substrate();
           if !genome.has_input_substrate(&random_gene_substrate) {
               let mut mutated_genome = genome.duplicate();
               mutated_genome.set_input(random_input, Some(random_gene_substrate));
               Some(mutated_genome)
           } else {
               // Do not allow two input substrates to be the same. This would violate the
               // specification of an I/O-substrate.
               None
           }
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
   fn mutate_input_dissociation(genome: &BinaryGenome) -> Option<BinaryGenome> {
       genome.get_random_input()
            .and_then(|random_in| genome.input(random_in)
                .and_then(|inner| inner)
                .and_then(|_| {
                    let mut mutated_genome = genome.duplicate();
                    // If there is an association, remove it.
                    mutated_genome.set_input(random_in, None);
                    Some(mutated_genome)
                })
            )
   }

   /// Duplicates the [`Genome`] and then associates a random output [`GeneSubstrate`] of the [`Genome`]
   /// with a random [`Substrate`] of a random [`Gene`].
   /// Returns the altered [`Genome`] if there are any output [`GeneSubstrate`].
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`GeneSubstrate`]: ./struct.GeneSubstrate.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_output_association(genome: &BinaryGenome) -> Option<BinaryGenome> {
       if genome.number_of_outputs() > 0 {
           // It was checked before if there are any outputs, so the unwrap must work.
           let random_output = genome.get_random_output().unwrap();
           let random_gene_substrate = genome.random_gene_substrate();
           if !genome.has_output_substrate(&random_gene_substrate) {
               let mut mutated_genome = genome.duplicate();
               mutated_genome.set_output(random_output, Some(random_gene_substrate));
               Some(mutated_genome)
           } else {
               // Do not allow two output substrates to be the same. This would violate the
               // specification of an I/O-substrate.
               None
           }
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
   fn mutate_output_dissociation(genome: &BinaryGenome) -> Option<BinaryGenome> {
       genome.get_random_output()
            .and_then(|random_out| genome.input(random_out)
                .and_then(|inner| inner)
                .and_then(|_| {
                    let mut mutated_genome = genome.duplicate();
                    // If there is an association, remove it.
                    mutated_genome.set_output(random_out, None);
                    Some(mutated_genome)
                })
            )
   }

   /// Duplicates the [`Genome`] and then fuses two random [`Gene`]s of the [`Genome`].
   /// The fusion product is added as a new [`Gene`].
   /// Returns the altered [`Genome`] if there are 2 or more different [`Gene`]s in the [`Genome`]
   /// and if the fusion process was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_gene_fusion(genome: &BinaryGenome) -> Option<BinaryGenome> {
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
   fn mutate_gene_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
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
   fn mutate_association_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let mut mutated_genome = genome.duplicate();
       let mut random_association = GeneAssociation::new(random_substrate());
       random_association.add_association(mutated_genome.random_gene_substrate());
       mutated_genome.add_association(random_association).and(Some(mutated_genome))
   }

   /// Duplicates the [`Genome`] and removes a random [`GeneAssociation`] from it.
   /// Returns the altered [`Genome`] if 1 or more [`GeneAssociation`]s were present.
   ///
   /// [`Gene`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   fn mutate_association_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
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
   fn mutate_association_substrate(genome: &BinaryGenome) -> Option<BinaryGenome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           let random_substrate = mutate_substrate_based_on(
               mutated_genome.association(random_association_index)
               // Unwrap must work since the index was just created.
               .unwrap()
               .substrate().len() / 8
           );
           // If the substrates are the same, report no change.
           if &random_substrate != mutated_genome.association(random_association_index).unwrap().substrate() {
               *mutated_genome.association_mut(random_association_index).unwrap().substrate_mut() = random_substrate;
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
   fn mutate_association_substrate_mutation_flip(genome: &BinaryGenome) -> Option<BinaryGenome> {
       genome.get_random_association().and_then(|random_association_index|{
           if genome.association(random_association_index).unwrap().substrate().len() > 0 {
               let mut mutated_genome = genome.duplicate();
               mutate_substrate_single_bit(mutated_genome.association_mut(random_association_index).unwrap().substrate_mut());
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
   fn mutate_association_substrate_mutation_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       genome.get_random_association().and_then(|random_association_index|{
               let mut mutated_genome = genome.duplicate();
               let mut mutated_substrate = mutated_genome.association(random_association_index).unwrap().substrate().clone();
               mutated_substrate = mutate_substrate_single_bit_insertion(mutated_substrate);
               *mutated_genome.association_mut(random_association_index).unwrap().substrate_mut() = mutated_substrate;
               Some(mutated_genome)
       })
   }

   /// Duplicates the [`Genome`] and randomly removes a bit to a [`Substrate`] of a [`GeneAssociation`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`GeneAssociation`]: ./struct.GeneAssociation.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_association_substrate_mutation_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       genome.get_random_association().and_then(|random_association_index|{
           if genome.association(random_association_index).unwrap().substrate().len() > 0 {
               let mut mutated_genome = genome.duplicate();
               let mut mutated_substrate = mutated_genome.association(random_association_index).unwrap().substrate().clone();
               mutated_substrate = mutate_substrate_single_bit_deletion(mutated_substrate);
               *mutated_genome.association_mut(random_association_index).unwrap().substrate_mut() = mutated_substrate;
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
   fn mutate_association_gene_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           let random_gene_substrate = mutated_genome.random_gene_substrate();
            mutated_genome.association_mut(random_association_index).unwrap().add_association(random_gene_substrate).and(Some(mutated_genome))
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
   fn mutate_association_gene_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       if let Some(random_association_index) = genome.get_random_association() {
           let mut mutated_genome = genome.duplicate();
           if mutated_genome.association(random_association_index).unwrap().number_of_associated_substrates() > 0 {
               if let Some(random_deletion_index) = mutated_genome.association(random_association_index).unwrap().get_random_gene_substrate() {
                   mutated_genome.association_mut(random_association_index).unwrap().remove_association(random_deletion_index);
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
   fn mutate_gene_substrate_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       mutated_genome.add_substrate_to_gene(random_gene_index, random_substrate()).and(Some(mutated_genome))
   }

   /// Duplicates the [`Genome`] and randomly removes a [`Substrate`] from an existing [`Gene`].
   /// Returns the altered [`Genome`] if more than 1 [`Substrate`] was present for the selected [`Gene`].
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let random_substrate = genome.random_gene_substrate();
       if genome.get_gene(random_substrate.gene()).number_of_substrates().get() > 1 {
           let mut mutated_genome = genome.duplicate();
           mutated_genome.remove_substrate(random_substrate);
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
   fn mutate_gene_substrate_mutation(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let mut mutated_genome = genome.duplicate();
       let random_substrate = mutated_genome.random_gene_substrate();
       let random_substrate_value = mutate_substrate_based_on(mutated_genome.get_substrate(random_substrate).unwrap().len() / 8);
       if &random_substrate_value != mutated_genome.get_substrate(random_substrate).unwrap() {
           *mutated_genome.get_substrate_mut(random_substrate).unwrap() = random_substrate_value;
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
   fn mutate_gene_substrate_mutation_flip(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let random_substrate = genome.random_gene_substrate();
       if genome.get_substrate(random_substrate).unwrap().len() > 0 {
           let mut mutated_genome = genome.duplicate();
           mutate_substrate_single_bit( mutated_genome.get_substrate_mut(random_substrate).unwrap());
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
   fn mutate_gene_substrate_mutation_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
           let mut mutated_genome = genome.duplicate();
           let random_substrate = mutated_genome.random_gene_substrate();
           let mut mutated_substrate = mutated_genome.get_substrate(random_substrate).unwrap().clone();
           mutated_substrate = mutate_substrate_single_bit_insertion(mutated_substrate);
           *mutated_genome.get_substrate_mut(random_substrate).unwrap() = mutated_substrate;
           Some(mutated_genome)
   }

   /// Duplicates the [`Genome`] and randomly removes a bit of a [`Substrate`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if the modification was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`Substrate`]: ../protein/struct.Substrate.html
   fn mutate_gene_substrate_mutation_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let random_substrate = genome.random_gene_substrate();
       if genome.get_substrate(random_substrate).unwrap().len() > 0 {
           let mut mutated_genome = genome.duplicate();
           let mut mutated_substrate = mutated_genome.get_substrate(random_substrate).unwrap().clone();
           mutated_substrate = mutate_substrate_single_bit_deletion(mutated_substrate);
           *mutated_genome.get_substrate_mut(random_substrate).unwrap() = mutated_substrate;
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
   fn mutate_gene_receptor_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       let random_receptor = mutated_genome.get_gene(random_gene_index).random_receptor();
       mutated_genome.get_gene_mut(random_gene_index).add_receptor(random_receptor).and(Some(mutated_genome))
   }

   /// Duplicates the [`Genome`] and randomly removes a [`GenomicReceptor`] from an existing [`Gene`].
   /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] were present.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   fn mutate_gene_receptor_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let random_gene_index = genome.get_random_gene();
       let mut mutated_genome = genome.duplicate();
        mutated_genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            mutated_genome.get_gene_mut(random_gene_index).remove_receptor(random_receptor);
            Some(mutated_genome)
       })
   }

   /// Duplicates the [`Genome`] and randomly adds a trigger to a [`GenomicReceptor`] of an existing [`Gene`].
   /// Returns the altered [`Genome`] if the addition was successful.
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
   fn mutate_gene_receptor_trigger_insertion(genome: &BinaryGenome) -> Option<BinaryGenome> {
       let mut mutated_genome = genome.duplicate();
       let random_gene_index = mutated_genome.get_random_gene();
       mutated_genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
       let random_substrate_index = mutated_genome.get_gene(random_gene_index).get_random_substrate();
       mutated_genome.get_gene_mut(random_gene_index)
            .receptor_mut(random_receptor)
            .unwrap()
            .add_trigger(random_substrate_index)
            .and(Some(mutated_genome))
       })
   }

    /// Duplicates the [`Genome`] and randomly removes a trigger from a [`GenomicReceptor`] of an existing [`Gene`].
    /// Returns the altered [`Genome`] if 1 or more [`GenomicReceptor`] and triggers were present.
    ///
    /// [`Gene`]: ./struct.Gene.html
    /// [`Genome`]: ./struct.Genome.html
    /// [`GenomicReceptor`]: ./struct.GenomicReceptor.html
    fn mutate_gene_receptor_trigger_deletion(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let random_gene_index = genome.get_random_gene();
        genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            genome.get_gene(random_gene_index).receptor(random_receptor).unwrap().get_random_trigger().and_then(|random_trigger| {
                let mut mutated_genome = genome.duplicate();
                mutated_genome.get_gene_mut(random_gene_index)
                    .receptor_mut(random_receptor)
                    .unwrap()
                    .remove_trigger(random_trigger);
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
    fn mutate_gene_receptor_state_mutation(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let random_gene_index = genome.get_random_gene();
        genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            let mut mutated_genome = genome.duplicate();
            let state: BinaryState = BinaryState::random();
            let substrates = (0..state.get_substrate_number()).map(|_| mutated_genome.get_gene(random_gene_index).get_random_substrate()).collect();
            mutated_genome.get_gene_mut(random_gene_index)
                .receptor_mut(random_receptor)
                .unwrap()
                .replace_state(state, substrates);
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
    fn mutate_gene_receptor_substrate_mutation(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let random_gene_index = genome.get_random_gene();
        genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            genome.get_gene(random_gene_index).receptor(random_receptor).unwrap().get_random_substrate().and_then(|random_substrate| {
                let mut mutated_genome = genome.duplicate();
                let new_substrate = mutated_genome.get_gene(random_gene_index).get_random_substrate();
                mutated_genome.get_gene_mut(random_gene_index)
                    .receptor_mut(random_receptor)
                    .unwrap()
                    .replace_substrate(random_substrate, new_substrate);
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
    fn mutate_gene_receptor_enzyme_mutation(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let random_gene_index = genome.get_random_gene();
        genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            let mut mutated_genome = genome.duplicate();
            let enzyme = mutated_genome.get_gene(random_gene_index).random_catalytic_centre();
            mutated_genome.get_gene_mut(random_gene_index)
                .receptor_mut(random_receptor)
                .unwrap()
                .replace_enzyme(enzyme);
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
    fn mutate_gene_catalytic_centre_educt_mutation(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let random_gene_index = genome.get_random_gene();
        genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            genome.get_gene(random_gene_index).receptor(random_receptor).unwrap().enzyme().get_random_educt().and_then(|random_educt| {
                let mut mutated_genome = genome.duplicate();
                let new_educt = mutated_genome.get_gene(random_gene_index).get_random_substrate();
                mutated_genome.get_gene_mut(random_gene_index)
                    .receptor_mut(random_receptor)
                    .unwrap()
                    .enzyme_mut()
                    .replace_educt(random_educt, new_educt);
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
    fn mutate_gene_catalytic_centre_product_mutation(genome: &BinaryGenome) -> Option<BinaryGenome> {
        let random_gene_index = genome.get_random_gene();
        genome.get_gene(random_gene_index).get_random_receptor().and_then(|random_receptor| {
            genome.get_gene(random_gene_index).receptor(random_receptor).unwrap().enzyme().get_random_product().and_then(|random_product| {
                let mut mutated_genome = genome.duplicate();
                let new_product = mutated_genome.get_gene(random_gene_index).get_random_substrate();
                mutated_genome.get_gene_mut(random_gene_index)
                    .receptor_mut(random_receptor)
                    .unwrap()
                    .enzyme_mut()
                    .replace_product(random_product, new_product);
                Some(mutated_genome)
            })
        })
    }

   /// Transfers the specified [`Gene`] into the specified acceptor [`Genome`] if possible and
   /// Returns the altered [`Genome`].
   ///
   /// # Parameters
   ///
   /// * `acceptor` - the [`Genome`] accepting a lateral gene transfer
   /// * `donor_gene` - the gene to be inserted into the acceptor [`Genome`]
   ///
   /// [`Gene`]: ./struct.Gene.html
   /// [`Genome`]: ./struct.Genome.html
   pub fn lateral_gene_transfer(acceptor: &BinaryGenome, donor_gene: Gene<BinaryReaction, BinaryState, BinarySubstrate>) -> Option<BinaryGenome> {
       let mut mutated_genome = acceptor.duplicate();
       mutated_genome.add_gene(donor_gene).and(Some(mutated_genome))
   }

   /// Mutates the specified [`Genome`] the specified amount of times. If any of the mutations
   /// was not successful `None` is returned.
   ///
   /// # Parameters
   ///
   /// * `number_of_mutations` - the number of times the genome should be mutated
   /// * `genome` - the [`Genome`] to mutate
   ///
   /// [`Genome`]: ./struct.Genome.html
    pub fn mutate_n_times(number_of_mutations: usize, genome: &BinaryGenome) -> Option<BinaryGenome> {
       if number_of_mutations == 0 {
           Some(genome.duplicate())
       } else {
           let mut mutated_genome = None;
           for i in 0..number_of_mutations {
               if i == 0 {
                   // If this is the first mutation, mutate the original genome.
                   mutated_genome = rand::random::<BinaryMutation>().mutate(genome);
               } else if mutated_genome.is_some() {
                   // If this is not the first mutation and the previous mutation was
                   // successful, continue mutating.
                   // Unwrapping must succeed as it was checked beforehand.
                   mutated_genome = rand::random::<BinaryMutation>().mutate(&mutated_genome.unwrap());
               } else {
                   // If this is not the first mutation and the previous mutation was not
                   // successful, stop mutating.
                   break;
               }
           }
           mutated_genome
       }
    }
}

/// Generates a random binary [`Substrate`].
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn random_substrate() -> BinarySubstrate {
    let length: usize = thread_rng().gen_range(RANDOM_SUBSTRATE_MIN_LENGTH, RANDOM_SUBSTRATE_MAX_LENGTH + 1);
    let random_bytes: Vec<u8> = (0..length).map(|_| thread_rng().gen::<u8>()).collect();
    BitVec::from_vec(random_bytes).into_boxed_bitslice()
}

/// Generates a random binary [`Substrate`] based on the specified length in byte.
/// The resulting [`Substrate`] will be of length 1 B if the base length was 0.
/// Otherwise the resulting length will be in range `base length / 2` to `base_length * 2`.
///
/// [`Substrate`]: ../protein/struct.Substrate.html
fn mutate_substrate_based_on(base_length: usize) -> BinarySubstrate {
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
    BitVec::from_vec(random_bytes).into_boxed_bitslice()
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
fn mutate_substrate_single_bit(base_substrate: &mut BinarySubstrate) {
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
fn mutate_substrate_single_bit_insertion(base_substrate: BinarySubstrate) -> BinarySubstrate {
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
fn mutate_substrate_single_bit_deletion(base_substrate: BinarySubstrate) -> BinarySubstrate {
    if base_substrate.is_empty() {
        panic!("No byte can be removed from an empty substrate.");
    }
    let mut mutated_substrate = BitVec::from_boxed_bitslice(base_substrate);
    let random_bit_index = thread_rng().gen_range(0, mutated_substrate.len());
    mutated_substrate.remove(random_bit_index);
    mutated_substrate.into_boxed_bitslice()
}

impl Distribution<BinaryMutation> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BinaryMutation {
        match rng.gen_range(0u8, 30) {
            0 => BinaryMutation::InputAssociation,
            1 => BinaryMutation::InputDissociation,
            2 => BinaryMutation::OutputAssociation,
            3 => BinaryMutation::OutputDissociation,
            4 => BinaryMutation::GeneFusion,
            5 => BinaryMutation::GeneDeletion,
            6 => BinaryMutation::GeneDuplication,
            7 => BinaryMutation::AssociationInsertion,
            8 => BinaryMutation::AssociationDeletion,
            9 => BinaryMutation::AssociationMutationSubstrate,
            10 => BinaryMutation::AssociationMutationGeneInsertion,
            11 => BinaryMutation::AssociationMutationGeneDeletion,
            12 => BinaryMutation::GeneMutationSubstrateInsertion,
            13 => BinaryMutation::GeneMutationSubstrateDeletion,
            14 => BinaryMutation::GeneMutationSubstrateMutation,
            15 => BinaryMutation::GeneMutationReceptorInsertion,
            16 => BinaryMutation::GeneMutationReceptorDeletion,
            17 => BinaryMutation::GeneMutationReceptorMutationTriggerInsertion,
            18 => BinaryMutation::GeneMutationReceptorMutationTriggerDeletion,
            19 => BinaryMutation::GeneMutationReceptorMutationStateMutation,
            20 => BinaryMutation::GeneMutationReceptorMutationSubstratesMutation,
            21 => BinaryMutation::GeneMutationReceptorMutationEnzymeMutation,
            22 => BinaryMutation::GeneMutationCatalyticCentreMutationEductMutation,
            23 => BinaryMutation::GeneMutationCatalyticCentreMutationProductMutation,
            24 => BinaryMutation::GeneMutationSubstrateMutationFlip,
            25 => BinaryMutation::GeneMutationSubstrateMutationInsertion,
            26 => BinaryMutation::GeneMutationSubstrateMutationDeletion,
            27 => BinaryMutation::AssociationMutationSubstrateMutationFlip,
            28 => BinaryMutation::AssociationMutationSubstrateMutationInsertion,
            29 => BinaryMutation::AssociationMutationSubstrateMutationDeletion,
            _ => panic!("A random number with no matching genomic mutation was created."),
        }
    }
}

impl GenomeMutation<BinaryReaction, BinaryState, BinarySubstrate> for BinaryMutation {
    /// Generates a mutated version of the specified [`Genome`] based on the kind of
    /// `BinaryMutation`.
    ///
    /// This will fail if the mutated [`Genome`] would be identical to the input or if the
    /// mutation is impossible for the specified [`Genome`].
    ///
    /// # Parameters
    ///
    /// `genome` - the base [`Genome`] to generate a mutated version of
    ///
    /// [`Genome`]: ./struct.Genome.html
    fn mutate(&self, genome: &BinaryGenome) -> Option<BinaryGenome> {
        match self {
            BinaryMutation::InputAssociation => BinaryMutation::mutate_input_association(genome),
            BinaryMutation::InputDissociation => BinaryMutation::mutate_input_dissociation(genome),
            BinaryMutation::OutputAssociation => BinaryMutation::mutate_output_association(genome),
            BinaryMutation::OutputDissociation => BinaryMutation::mutate_output_dissociation(genome),
            BinaryMutation::GeneFusion => BinaryMutation::mutate_gene_fusion(genome),
            BinaryMutation::GeneDeletion => BinaryMutation::mutate_gene_deletion(genome),
            BinaryMutation::GeneDuplication => BinaryMutation::mutate_gene_duplication(genome),
            BinaryMutation::AssociationInsertion => BinaryMutation::mutate_association_insertion(genome),
            BinaryMutation::AssociationDeletion => BinaryMutation::mutate_association_deletion(genome),
            BinaryMutation::AssociationMutationSubstrate => BinaryMutation::mutate_association_substrate(genome),
            BinaryMutation::AssociationMutationSubstrateMutationFlip => BinaryMutation::mutate_association_substrate_mutation_flip(genome),
            BinaryMutation::AssociationMutationSubstrateMutationInsertion => BinaryMutation::mutate_association_substrate_mutation_insertion(genome),
            BinaryMutation::AssociationMutationSubstrateMutationDeletion => BinaryMutation::mutate_association_substrate_mutation_deletion(genome),
            BinaryMutation::AssociationMutationGeneInsertion => BinaryMutation::mutate_association_gene_insertion(genome),
            BinaryMutation::AssociationMutationGeneDeletion => BinaryMutation::mutate_association_gene_deletion(genome),
            BinaryMutation::GeneMutationSubstrateInsertion => BinaryMutation::mutate_gene_substrate_insertion(genome),
            BinaryMutation::GeneMutationSubstrateDeletion => BinaryMutation::mutate_gene_substrate_deletion(genome),
            BinaryMutation::GeneMutationSubstrateMutation => BinaryMutation::mutate_gene_substrate_mutation(genome),
            BinaryMutation::GeneMutationSubstrateMutationFlip => BinaryMutation::mutate_gene_substrate_mutation_flip(genome),
            BinaryMutation::GeneMutationSubstrateMutationInsertion => BinaryMutation::mutate_gene_substrate_mutation_insertion(genome),
            BinaryMutation::GeneMutationSubstrateMutationDeletion => BinaryMutation::mutate_gene_substrate_mutation_deletion(genome),
            BinaryMutation::GeneMutationReceptorInsertion => BinaryMutation::mutate_gene_receptor_insertion(genome),
            BinaryMutation::GeneMutationReceptorDeletion => BinaryMutation::mutate_gene_receptor_deletion(genome),
            BinaryMutation::GeneMutationReceptorMutationTriggerInsertion => BinaryMutation::mutate_gene_receptor_trigger_insertion(genome),
            BinaryMutation::GeneMutationReceptorMutationTriggerDeletion => BinaryMutation::mutate_gene_receptor_trigger_deletion(genome),
            BinaryMutation::GeneMutationReceptorMutationStateMutation => BinaryMutation::mutate_gene_receptor_state_mutation(genome),
            BinaryMutation::GeneMutationReceptorMutationSubstratesMutation => BinaryMutation::mutate_gene_receptor_substrate_mutation(genome),
            BinaryMutation::GeneMutationReceptorMutationEnzymeMutation => BinaryMutation::mutate_gene_receptor_enzyme_mutation(genome),
            BinaryMutation::GeneMutationCatalyticCentreMutationEductMutation => BinaryMutation::mutate_gene_catalytic_centre_educt_mutation(genome),
            BinaryMutation::GeneMutationCatalyticCentreMutationProductMutation => BinaryMutation::mutate_gene_catalytic_centre_product_mutation(genome),
            // BinaryMutation::LateralGeneTransfer => None, //TODO: implement a global gene pool
        }
    }

    fn random() -> Self {
        thread_rng().gen()
    }
}
