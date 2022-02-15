use crate::evolution::helper::noop::{
    NoOpOutputElement, NoOpOutputSensor, NoOpReaction, NoOpState, NoOpSubstrate,
};

use super::*;

#[test]
/// Tests if the function `new` correctly creates a [`GenomicOutputSensor`].
fn test_new() {
    let output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        None,
    ];
    let finish_substrate = Some(GeneSubstrate::new(3, 0));
    let output: NoOpOutputSensor = ();
    let output_sensor: GenomicOutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = GenomicOutputSensor::new(output_substrates.clone(), finish_substrate, output.clone());
    assert_eq!(output_sensor.output_substrates(), &output_substrates);
    assert_eq!(output_sensor.finish_substrate(), &finish_substrate);
    assert_eq!(output_sensor.output(), &output);
}

#[test]
/// Tests if the function `number_of_associated_outputs` correctly returns the number of associated outputs.
fn test_number_of_associated_outputs() {
    let output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        None,
    ];
    let finish_substrate = None;
    let output: NoOpOutputSensor = ();
    let output_sensor: GenomicOutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = GenomicOutputSensor::new(output_substrates, finish_substrate, output);
    assert_eq!(output_sensor.number_of_associated_outputs(), 2);
}

#[test]
/// Tests if the function `adjust_after_gene_removal` correctly adjusts the gene indices.
fn test_adjust_after_gene_removal() {
    let output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        Some(GeneSubstrate::new(1, 7)),
        Some(GeneSubstrate::new(2, 1)),
        Some(GeneSubstrate::new(3, 2)),
    ];
    let finish_substrate = Some(GeneSubstrate::new(3, 3));
    let expected_output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        None,
        None,
        None,
        None,
        Some(GeneSubstrate::new(1, 1)),
        Some(GeneSubstrate::new(2, 2)),
    ];
    let expected_finish_substrate = Some(GeneSubstrate::new(2, 3));
    let output: NoOpOutputSensor = ();
    let mut output_sensor: GenomicOutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = GenomicOutputSensor::new(output_substrates, finish_substrate, output);
    output_sensor.adjust_after_gene_removal(1);
    assert_eq!(output_sensor.output_substrates(), &expected_output_substrates);
    assert_eq!(output_sensor.finish_substrate(), &expected_finish_substrate);
}

#[test]
/// Tests if the function `adjust_after_gene_substrate_removal` correctly adjusts the gene indices.
fn test_adjust_after_gene_substrate_removal() {
    let output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        Some(GeneSubstrate::new(1, 7)),
        Some(GeneSubstrate::new(2, 1)),
        Some(GeneSubstrate::new(3, 2)),
    ];
    let finish_substrate = Some(GeneSubstrate::new(1, 5));
    let expected_output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        None,
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 6)),
        Some(GeneSubstrate::new(2, 1)),
        Some(GeneSubstrate::new(3, 2)),
    ];
    let expected_finish_substrate = Some(GeneSubstrate::new(1, 4));
    let output: NoOpOutputSensor = ();
    let mut output_sensor: GenomicOutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = GenomicOutputSensor::new(output_substrates, finish_substrate, output);
    output_sensor.adjust_after_gene_substrate_removal(GeneSubstrate::new(1, 3));
    assert_eq!(output_sensor.output_substrates(), &expected_output_substrates);
    assert_eq!(output_sensor.finish_substrate(), &expected_finish_substrate);
}

#[test]
/// Tests if the function `validate` correctly validates the internal state.
fn test_validate() {
    let output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        Some(GeneSubstrate::new(1, 7)),
        Some(GeneSubstrate::new(2, 1)),
    ];
    let finish_substrate = Some(GeneSubstrate::new(3, 3));
    let expected_output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(1, 4)),
        None,
        Some(GeneSubstrate::new(2, 1)),
    ];
    let expected_finish_substrate = None;
    let output: NoOpOutputSensor = ();
    let mut output_sensor: GenomicOutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = GenomicOutputSensor::new(output_substrates, finish_substrate, output);
    let genes = vec![
        Gene::new(vec![()]),
        Gene::new(vec![(), (), (), (), (), ()]),
        Gene::new(vec![(), ()]),
    ];
    output_sensor.validate(&genes);
    assert_eq!(output_sensor.output_substrates(), &expected_output_substrates);
    assert_eq!(output_sensor.finish_substrate(), &expected_finish_substrate);
}

#[test]
/// Tests if the function `translate` correctly translates the genomic information to the protein level.
fn test_translate() {
    let output_substrates = vec![
        Some(GeneSubstrate::new(0, 0)),
        Some(GeneSubstrate::new(1, 0)),
        Some(GeneSubstrate::new(1, 3)),
        Some(GeneSubstrate::new(2, 1)),
    ];
    let finish_substrate = Some(GeneSubstrate::new(3, 3));
    let output: NoOpOutputSensor = ();
    let output_sensor: GenomicOutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = GenomicOutputSensor::new(output_substrates, finish_substrate, output);
    let substrate_lookup = HashMap::from([
        (GeneSubstrate::new(0, 0), new_noop_substrate()),
        (GeneSubstrate::new(1, 0), new_noop_substrate()),
        (GeneSubstrate::new(1, 3), new_noop_substrate()),
        (GeneSubstrate::new(2, 1), new_noop_substrate()),
        (GeneSubstrate::new(3, 3), new_noop_substrate()),
    ]);
    let protein = output_sensor.translate(&substrate_lookup);
    assert_eq!(protein.output_substrates().len(), output_sensor.output_substrates().len());
    for (i, substrate_option) in output_sensor.output_substrates().iter().enumerate() {
        if let Some(substrate) = substrate_option {
            assert!(Rc::downgrade(substrate_lookup.get(substrate).unwrap())
                .ptr_eq(&protein.output_substrates()[i].clone().unwrap()));
        } else {
            assert!(protein.output_substrates()[i].is_none());
        }
    }
    assert_eq!(
        SubstrateType::ConventionalSubstrate,
        get_substrate_reference_option_type(substrate_lookup.get(&GeneSubstrate::new(0, 0)))
    );
    assert_eq!(
        SubstrateType::OutputFinishSubstrate,
        get_substrate_reference_option_type(substrate_lookup.get(&GeneSubstrate::new(3, 3)))
    );
    assert_eq!(
        SubstrateType::ConventionalSubstrate,
        get_substrate_reference_option_type(substrate_lookup.get(&GeneSubstrate::new(1, 0)))
    );
    assert_eq!(
        SubstrateType::ConventionalSubstrate,
        get_substrate_reference_option_type(substrate_lookup.get(&GeneSubstrate::new(1, 3)))
    );
    assert_eq!(
        SubstrateType::ConventionalSubstrate,
        get_substrate_reference_option_type(substrate_lookup.get(&GeneSubstrate::new(2, 1)))
    );
}

fn new_noop_substrate() -> Rc<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>> {
    Rc::new(RefCell::new(Substrate::new((), SubstrateType::ConventionalSubstrate)))
}

fn get_substrate_reference_option_type(
    substrate_reference: Option<&Rc<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
) -> SubstrateType {
    substrate_reference
        .unwrap()
        .borrow()
        .substrate_type()
        .clone()
}
