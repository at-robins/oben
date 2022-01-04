use crate::evolution::{
    gene::{Gene, GeneSubstrate, GenomicInputSensor},
    helper::testing::{TestGenome, TestInformation},
};

use super::*;

#[test]
/// Tests if the function `new` of the [`Mutation`] struct correctly creates a [`Mutation`].
fn test_mutation_new() {
    let chance = 0.17;
    let mutation = Mutation::new(chance, |genome: &TestGenome| Some(genome.duplicate()));
    assert_ulps_eq!(chance, mutation.mutation_chance().value());
    let unmutated_genome: TestGenome = Genome::new(
        GenomicInputSensor::default(),
        Vec::new(),
        vec![Gene::new(vec![TestInformation { value: 0 }])],
    );
    assert_eq!(mutation.mutate(&unmutated_genome), Some(unmutated_genome));
}

#[test]
/// Tests if the function `eq` of the [`Mutation`] struct correctly compares [`Mutation`]s.
fn test_mutation_eq() {
    let chance = 0.17;
    let mutation_a = Mutation::new(chance, |genome: &TestGenome| Some(genome.duplicate()));
    let mutation_b = Mutation::new(chance, |genome: &TestGenome| Some(genome.duplicate()));
    assert_eq!(mutation_a, mutation_a);
    assert_eq!(mutation_b, mutation_b);
    assert_ne!(mutation_a, mutation_b);
}

#[test]
/// Tests if the function `mutate` of the [`Mutation`] struct correctly mutates a [`Genome`].
fn test_mutation_mutate() {
    let chance = 0.17;
    let mutation = Mutation::new(chance, |genome: &TestGenome| {
        let mut mutated_genome = genome.duplicate();
        *mutated_genome
            .get_substrate_mut(GeneSubstrate::new(0, 0))
            .unwrap() = TestInformation { value: 5 };
        Some(mutated_genome)
    });
    let unmutated_genome: TestGenome = Genome::new(
        GenomicInputSensor::default(),
        Vec::new(),
        vec![Gene::new(vec![TestInformation { value: 0 }])],
    );
    let expected_mutated_genome: TestGenome = Genome::new(
        GenomicInputSensor::default(),
        Vec::new(),
        vec![Gene::new(vec![TestInformation { value: 5 }])],
    );
    let mutated_genome = mutation.mutate(&unmutated_genome);
    assert_eq!(Some(expected_mutated_genome), mutated_genome);
    assert_ne!(Some(unmutated_genome), mutated_genome);
}
