use crate::evolution::helper::noop::{
    NoOpInputElement, NoOpInputSensor, NoOpReaction, NoOpState, NoOpSubstrate,
};

use super::*;

#[test]
/// Tests if the function `new` correctly creates a [`GenomicInputSensor`].
fn test_new() {
    let input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        None,
    ];
    let feedback_substrates = vec![
        Some(GeneSubstrate::new(3, 0)),
        None,
        Some(GeneSubstrate::new(1, 4)),
    ];
    let input: NoOpInputSensor = ();
    let input_sensor: GenomicInputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = GenomicInputSensor::new(
        input_substrates.clone(),
        feedback_substrates.clone(),
        input.clone(),
    );
    assert_eq!(input_sensor.input_substrates(), &input_substrates);
    assert_eq!(input_sensor.feedback_substrates(), &feedback_substrates);
    assert_eq!(input_sensor.input(), &input);
}

#[test]
/// Tests if the function `number_of_associated_inputs` correctly returns the number of associated inputs.
fn test_number_of_associated_inputs() {
    let input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        None,
    ];
    let feedback_substrates = Vec::new();
    let input: NoOpInputSensor = ();
    let input_sensor: GenomicInputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = GenomicInputSensor::new(input_substrates, feedback_substrates, input);
    assert_eq!(input_sensor.number_of_associated_inputs(), 2);
}

#[test]
/// Tests if the function `adjust_after_gene_removal` correctly adjusts the gene indices.
fn test_adjust_after_gene_removal() {
    let input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        Some(GeneSubstrate::new(1, 7)),
        Some(GeneSubstrate::new(2, 1)),
        Some(GeneSubstrate::new(3, 2)),
    ];
    let feedback_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        Some(GeneSubstrate::new(3, 3)),
        Some(GeneSubstrate::new(4, 4)),
    ];
    let expected_input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        None,
        None,
        None,
        Some(GeneSubstrate::new(1, 1)),
        Some(GeneSubstrate::new(2, 2)),
    ];
    let expected_feedback_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        Some(GeneSubstrate::new(2, 3)),
        Some(GeneSubstrate::new(3, 4)),
    ];
    let input: NoOpInputSensor = ();
    let mut input_sensor: GenomicInputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = GenomicInputSensor::new(input_substrates, feedback_substrates, input);
    input_sensor.adjust_after_gene_removal(1);
    assert_eq!(input_sensor.input_substrates(), &expected_input_substrates);
    assert_eq!(input_sensor.feedback_substrates(), &expected_feedback_substrates);
}

#[test]
/// Tests if the function `adjust_after_gene_substrate_removal` correctly adjusts the gene indices.
fn test_adjust_after_gene_substrate_removal() {
    let input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        Some(GeneSubstrate::new(1, 7)),
        Some(GeneSubstrate::new(2, 1)),
        Some(GeneSubstrate::new(3, 2)),
    ];
    let feedback_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        Some(GeneSubstrate::new(3, 3)),
        Some(GeneSubstrate::new(4, 4)),
    ];
    let expected_input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        None,
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 6)),
        Some(GeneSubstrate::new(2, 1)),
        Some(GeneSubstrate::new(3, 2)),
    ];
    let expected_feedback_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        Some(GeneSubstrate::new(3, 3)),
        Some(GeneSubstrate::new(4, 4)),
    ];
    let input: NoOpInputSensor = ();
    let mut input_sensor: GenomicInputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = GenomicInputSensor::new(input_substrates, feedback_substrates, input);
    input_sensor.adjust_after_gene_substrate_removal(GeneSubstrate::new(1, 3));
    assert_eq!(input_sensor.input_substrates(), &expected_input_substrates);
    assert_eq!(input_sensor.feedback_substrates(), &expected_feedback_substrates);
}

#[test]
/// Tests if the function `validate` correctly validates the internal state.
fn test_validate() {
    let input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        Some(GeneSubstrate::new(1, 7)),
        Some(GeneSubstrate::new(2, 1)),
    ];
    let feedback_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        Some(GeneSubstrate::new(3, 3)),
        Some(GeneSubstrate::new(4, 4)),
    ];
    let expected_input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        None,
        Some(GeneSubstrate::new(2, 1)),
    ];
    let expected_feedback_substrates = vec![Some(GeneSubstrate::new(0, 0)), None, None, None];
    let input: NoOpInputSensor = ();
    let mut input_sensor: GenomicInputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = GenomicInputSensor::new(input_substrates, feedback_substrates, input);
    let genes = vec![
        Gene::new(vec![()]),
        Gene::new(vec![(), (), (), (), (), ()]),
        Gene::new(vec![(), ()]),
    ];
    input_sensor.validate(&genes);
    assert_eq!(input_sensor.input_substrates(), &expected_input_substrates);
    assert_eq!(input_sensor.feedback_substrates(), &expected_feedback_substrates);
}

#[test]
/// Tests if the function `translate` correctly translates the genomic information to the protein level.
fn test_translate() {
    let input_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(2, 1)),
    ];
    let feedback_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        Some(GeneSubstrate::new(3, 3)),
        Some(GeneSubstrate::new(4, 4)),
    ];
    let input: NoOpInputSensor = ();
    let input_sensor: GenomicInputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = GenomicInputSensor::new(input_substrates, feedback_substrates, input);
    let substrate_lookup = HashMap::from([
        (GeneSubstrate::new(0, 0), new_noop_substrate()),
        (GeneSubstrate::new(1, 0), new_noop_substrate()),
        (GeneSubstrate::new(1, 3), new_noop_substrate()),
        (GeneSubstrate::new(2, 1), new_noop_substrate()),
        (GeneSubstrate::new(3, 3), new_noop_substrate()),
        (GeneSubstrate::new(4, 4), new_noop_substrate()),
    ]);
    let protein = input_sensor.translate(&substrate_lookup);
    assert_eq!(protein.input_substrates().len(), input_sensor.input_substrates().len());
    for (i, substrate_option) in input_sensor.input_substrates().iter().enumerate() {
        if let Some(substrate) = substrate_option {
            assert!(Rc::downgrade(substrate_lookup.get(substrate).unwrap()).ptr_eq(&protein.input_substrates()[i].clone().unwrap()));
        } else {
            assert!(protein.input_substrates()[i].is_none());
        }
    }
    assert_eq!(protein.feedback_substrates().len(), input_sensor.feedback_substrates().len());
    for (i, substrate_option) in input_sensor.feedback_substrates().iter().enumerate() {
        if let Some(substrate) = substrate_option {
            assert!(Rc::downgrade(substrate_lookup.get(substrate).unwrap()).ptr_eq(&protein.feedback_substrates()[i].clone().unwrap()));
        } else {
            assert!(protein.feedback_substrates()[i].is_none());
        }
    }
}

fn new_noop_substrate() -> Rc<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>> {
    Rc::new(RefCell::new(Substrate::new(())))
}
