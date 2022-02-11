use super::*;

// #[test]
// /// Tests if the function `remove_associations_with_gene_from_substrates` correctly removes associations.
// fn test_remove_associations_with_gene_from_substrates() {
//     let mut substrates = vec![
//         Some(GeneSubstrate::new(0, 0)),
//         Some(GeneSubstrate::new(1, 0)),
//         Some(GeneSubstrate::new(1, 1)),
//         None,
//         Some(GeneSubstrate::new(1, 2)),
//         Some(GeneSubstrate::new(2, 0)),
//         None,
//     ];
//     let expected_substrates = vec![
//         Some(GeneSubstrate::new(0, 0)),
//         None,
//         None,
//         None,
//         None,
//         Some(GeneSubstrate::new(2, 0)),
//         None,
//     ];
//     remove_associations_with_gene_from_substrates(&mut substrates, 1);
//     assert_eq!(substrates, expected_substrates);
// }
