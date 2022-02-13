use crate::evolution::helper::noop::{NoOpGene, NoOpGenome};

use super::*;

#[test]
/// Tests if the function `validate_associations` of the `Genome` struct.
fn test_validate_associations() {
    // Test invalid substrates.
    {
        let gene: NoOpGene = Gene::new(vec![()]);
        let invalid_gene_substrate = GeneSubstrate {
            gene: 1,
            substrate: 1,
        };
        let input_sensor =
            GenomicInputSensor::new(vec![Some(invalid_gene_substrate)], HashMap::new(), ());
        let output_sensor =
            GenomicOutputSensor::new(vec![Some(invalid_gene_substrate)], Some(invalid_gene_substrate), ());
        let invalid_association = GeneAssociation {
            substrate: (),
            associations: vec![invalid_gene_substrate.clone()],
        };
        let mut genome: NoOpGenome = Genome {
            phantom_reaction: PhantomData,
            phantom_state: PhantomData,
            phantom_information: PhantomData,
            phantom_input_element: PhantomData,
            phantom_input_sensor: PhantomData,
            phantom_output_element: PhantomData,
            phantom_output_sensor: PhantomData,
            input: input_sensor,
            output: output_sensor,
            genes: vec![gene],
            associations: vec![invalid_association],
        };
        genome.validate_associations();
        assert_eq!(genome.input.input_substrates(), &vec!(None));
        assert_eq!(genome.output.output_substrates(), &vec!(None));
        assert_eq!(genome.output.finish_substrate(), &None);
        assert_eq!(
            genome.associations,
            vec!(GeneAssociation {
                substrate: (),
                associations: Vec::new()
            })
        );
    }
    // Test duplicate I/O-substrates.
    {
        let gene: NoOpGene = Gene::new(vec![()]);
        let gene_substrate = GeneSubstrate {
            gene: 0,
            substrate: 0,
        };
        let input_sensor = GenomicInputSensor::new(
            vec![Some(gene_substrate.clone()), Some(gene_substrate.clone())],
            HashMap::new(),
            (),
        );
        let output_sensor = GenomicOutputSensor::new(
            vec![Some(gene_substrate.clone()), Some(gene_substrate.clone())],
            None,
            (),
        );
        let mut genome: NoOpGenome = Genome::new(
            input_sensor,
            output_sensor,
            vec![gene],
        );
        genome.validate_associations();
        assert_eq!(genome.input.input_substrates(), &vec!(Some(gene_substrate.clone()), None));
        assert_eq!(genome.output.output_substrates(), &vec!(Some(gene_substrate.clone()), None));
    }
}

#[test]
/// Tests if the function `has_substrate` of the `Genome` struct.
fn test_has_substrate() {
    {
        let gene: NoOpGene = Gene::new(vec![()]);
        let input_sensor = GenomicInputSensor::new(Vec::new(), HashMap::new(), ());
        let output_sensor = GenomicOutputSensor::new(Vec::new(), None, ());
        let genome: NoOpGenome = Genome::new(input_sensor, output_sensor, vec![gene]);
        let positive = GeneSubstrate {
            gene: 0,
            substrate: 0,
        };
        let negative_a = GeneSubstrate {
            gene: 0,
            substrate: 1,
        };
        let negative_b = GeneSubstrate {
            gene: 1,
            substrate: 1,
        };
        assert!(has_substrate(&genome.genes, &positive));
        assert!(!has_substrate(&genome.genes, &negative_a));
        assert!(!has_substrate(&genome.genes, &negative_b));
    }
}
