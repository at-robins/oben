use super::*;

#[test]
/// Tests if the function `validate_associations` of the `Genome` struct.
fn test_validate_associations() {
    // Test invalid substrates.
    {
        let mut gene = Gene{substrates: Vec::new(), receptors: Vec::new()};
            gene.add_substrate(BitBox::empty());
        let invalid_gene_substrate = GeneSubstrate{gene: 1, substrate: 1};
        let invalid_association = GeneAssociation{
            substrate: BitBox::empty(),
            associations: vec!(invalid_gene_substrate.clone())
        };
        let mut genome = Genome {
            input: vec!(Some(invalid_gene_substrate.clone())),
            output: vec!(Some(invalid_gene_substrate.clone())),
            genes: vec!(gene),
            associations: vec!(invalid_association),
        };
        genome.validate_associations();
        assert_eq!(genome.input, vec!(None));
        assert_eq!(genome.output, vec!(None));
        assert_eq!(genome.associations, vec!(GeneAssociation{substrate: BitBox::empty(), associations: Vec::new()}));
    }
    // Test duplicate I/O-substrates.
    {
        let mut gene = Gene{substrates: Vec::new(), receptors: Vec::new()};
            gene.add_substrate(BitBox::empty());
        let gene_substrate = GeneSubstrate{gene: 0, substrate: 0};
        let mut genome = Genome {
            input: vec!(Some(gene_substrate.clone()), Some(gene_substrate.clone())),
            output: vec!(Some(gene_substrate.clone()), Some(gene_substrate.clone())),
            genes: vec!(gene),
            associations: Vec::new(),
        };
        genome.validate_associations();
        assert_eq!(genome.input, vec!(Some(gene_substrate.clone()), None));
        assert_eq!(genome.output, vec!(Some(gene_substrate.clone()), None));
    }
}

#[test]
/// Tests if the function `has_substrate` of the `Genome` struct.
fn test_has_substrate() {
    {
        let mut gene = Gene{substrates: Vec::new(), receptors: Vec::new()};
            gene.add_substrate(BitBox::empty());
        let genome = Genome {
            input: Vec::new(),
            output: Vec::new(),
            genes: vec!(gene),
            associations: Vec::new(),
        };
        let positive = GeneSubstrate{gene: 0, substrate: 0};
        let negative_a = GeneSubstrate{gene: 0, substrate: 1};
        let negative_b = GeneSubstrate{gene: 1, substrate: 1};
        assert!(Genome::has_substrate(&genome.genes, &positive));
        assert!(!Genome::has_substrate(&genome.genes, &negative_a));
        assert!(!Genome::has_substrate(&genome.genes, &negative_b));
    }
}
