use crate::evolution::{
    gene::{Gene, GeneSubstrate, GenomicInputSensor, GenomicOutputSensor},
    helper::{
        noop::{NoOpInputElement, NoOpOutputElement},
        testing::{TestGenome, TestInformation, TestInput, TestReaction, TestState, TestOutput},
    },
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
        GenomicOutputSensor::default(),
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
        GenomicOutputSensor::default(),
        vec![Gene::new(vec![TestInformation { value: 0 }])],
    );
    let expected_mutated_genome: TestGenome = Genome::new(
        GenomicInputSensor::default(),
        GenomicOutputSensor::default(),
        vec![Gene::new(vec![TestInformation { value: 5 }])],
    );
    let mutated_genome = mutation.mutate(&unmutated_genome);
    assert_eq!(Some(expected_mutated_genome), mutated_genome);
    assert_ne!(Some(unmutated_genome), mutated_genome);
}

#[test]
/// Tests if the function `test_number_of_mutation_events` of the [`Mutation`] struct correctly simulates the mutation process.
fn test_number_of_mutation_events() {
    {
        let chance = 0.0;
        let mutation = Mutation::new(chance, |genome: &TestGenome| Some(genome.duplicate()));
        assert_eq!(mutation.number_of_mutation_events(), 0);
    }
    {
        let chance = 1.0;
        let mutation = Mutation::new(chance, |genome: &TestGenome| Some(genome.duplicate()));
        assert_eq!(mutation.number_of_mutation_events(), MAX_MUTATION_EVENTS);
    }
    {
        let chance = 0.05;
        let mutation = Mutation::new(chance, |genome: &TestGenome| Some(genome.duplicate()));
        let n = mutation.number_of_mutation_events();
        assert!(n <= MAX_MUTATION_EVENTS);
    }
}

#[test]
/// Tests if the function `new` of the [`MutationCompendium`] struct correctly creates an empty [`MutationCompendium`].
fn test_mutation_compendium_new() {
    let compendium: MutationCompendium<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
        NoOpOutputElement,
        TestOutput,
    > = MutationCompendium::new();
    let unmutated_genome: TestGenome = Genome::new(
        GenomicInputSensor::default(),
        GenomicOutputSensor::default(),
        vec![Gene::new(vec![TestInformation { value: 0 }])],
    );
    assert_eq!(compendium.size(), 0);
    assert!(compendium.mutate(&unmutated_genome).is_none());
}

#[test]
/// Tests if the function `add` of the [`MutationCompendium`] struct correctly adds [`Mutation`]s.
fn test_mutation_compendium_add() {
    let mut compendium: MutationCompendium<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
        NoOpOutputElement,
        TestOutput,
    > = MutationCompendium::new();
    assert_eq!(compendium.size(), 0);
    compendium.add(Mutation::new(0.2, |genome: &TestGenome| Some(genome.duplicate())));
    assert_eq!(compendium.size(), 1);
    compendium.add(Mutation::new(0.3, |genome: &TestGenome| Some(genome.duplicate())));
    assert_eq!(compendium.size(), 2);
}

#[test]
/// Tests if the function `mutate` of the [`MutationCompendium`] struct correctly mutates a [`Genome`].
fn test_mutation_compendium_mutate() {
    {
        // Zero percent chance mutations are never applied.
        let compendium: MutationCompendium<
            TestReaction,
            TestState,
            TestInformation,
            NoOpInputElement,
            TestInput,
            NoOpOutputElement,
            TestOutput,
        > = vec![Mutation::new(0.0, |genome: &TestGenome| {
            Some(genome.duplicate())
        })]
        .into();
        let unmutated_genome: TestGenome = Genome::new(
            GenomicInputSensor::default(),
            GenomicOutputSensor::default(),
            vec![Gene::new(vec![TestInformation { value: 0 }])],
        );
        assert_eq!(compendium.size(), 1);
        assert!(compendium.mutate(&unmutated_genome).is_none());
    }
    {
        // 100 percent chance mutations are always applied.
        let mutation = Mutation::new(1.0, |genome: &TestGenome| {
            let mut mutated_genome = genome.duplicate();
            let substrate = mutated_genome
                .get_substrate_mut(GeneSubstrate::new(0, 0))
                .unwrap();
            *substrate = TestInformation {
                value: substrate.value + 1,
            };
            Some(mutated_genome)
        });
        let compendium: MutationCompendium<
            TestReaction,
            TestState,
            TestInformation,
            NoOpInputElement,
            TestInput,
            NoOpOutputElement,
            TestOutput,
        > = vec![mutation].into();
        let unmutated_genome: TestGenome = Genome::new(
            GenomicInputSensor::default(),
            GenomicOutputSensor::default(),
            vec![Gene::new(vec![TestInformation { value: 0 }])],
        );
        let expected_mutated_genome: TestGenome = Genome::new(
            GenomicInputSensor::default(),
            GenomicOutputSensor::default(),
            vec![Gene::new(vec![TestInformation {
                value: MAX_MUTATION_EVENTS,
            }])],
        );
        assert_eq!(compendium.size(), 1);
        assert_eq!(Some(expected_mutated_genome), compendium.mutate(&unmutated_genome));
    }
    {
        // Other mutations are always applied on a statistical basis.
        let mutation_a = Mutation::new(0.1, |genome: &TestGenome| {
            let mut mutated_genome = genome.duplicate();
            let substrate = mutated_genome
                .get_substrate_mut(GeneSubstrate::new(0, 0))
                .unwrap();
            if substrate.value < MAX_MUTATION_EVENTS {
                *substrate = TestInformation {
                    value: substrate.value + 1,
                };
            }
            Some(mutated_genome)
        });
        let mutation_b = Mutation::new(0.2, |genome: &TestGenome| {
            let mut mutated_genome = genome.duplicate();
            let substrate = mutated_genome
                .get_substrate_mut(GeneSubstrate::new(0, 0))
                .unwrap();
            if substrate.value < MAX_MUTATION_EVENTS {
                *substrate = TestInformation {
                    value: substrate.value + 1,
                };
            }
            Some(mutated_genome)
        });
        let mutation_c = Mutation::new(0.3, |genome: &TestGenome| {
            let mut mutated_genome = genome.duplicate();
            let substrate = mutated_genome
                .get_substrate_mut(GeneSubstrate::new(0, 0))
                .unwrap();
            if substrate.value < MAX_MUTATION_EVENTS {
                *substrate = TestInformation {
                    value: substrate.value + 1,
                };
            }
            Some(mutated_genome)
        });
        let compendium: MutationCompendium<
            TestReaction,
            TestState,
            TestInformation,
            NoOpInputElement,
            TestInput,
            NoOpOutputElement,
            TestOutput,
        > = vec![mutation_a, mutation_b, mutation_c].into();
        let unmutated_genome: TestGenome = Genome::new(
            GenomicInputSensor::default(),
            GenomicOutputSensor::default(),
            vec![Gene::new(vec![TestInformation { value: 0 }])],
        );
        assert_eq!(compendium.size(), 3);
        if let Some(mutated_genome) = compendium.mutate(&unmutated_genome) {
            assert!(
                mutated_genome
                    .get_substrate(GeneSubstrate::new(0, 0))
                    .unwrap()
                    .value
                    <= MAX_MUTATION_EVENTS
            );
        }
    }
}
