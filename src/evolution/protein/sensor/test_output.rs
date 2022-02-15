use crate::evolution::{
    helper::noop::{NoOpOutputElement, NoOpOutputSensor, NoOpReaction, NoOpState, NoOpSubstrate},
    protein::SubstrateType,
};

use super::*;

#[test]
/// Tests if the function `new` correctly creates an [`OutputSensor`].
fn test_new() {
    let output_substrates = vec![Some(new_noop_substrate()), Some(new_noop_substrate()), None];
    let output_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
    > = output_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let output: NoOpOutputSensor = ();
    let output_sensor: OutputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpOutputElement,
        NoOpOutputSensor,
    > = OutputSensor::new(output.clone(), output_substrates_weak.clone());
    assert!(vec_pointer_equality_optional_weak(
        output_sensor.output_substrates(),
        &output_substrates_weak
    ));
    assert_eq!(output_sensor.output(), &output);
}

fn new_noop_substrate() -> Rc<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>> {
    Rc::new(RefCell::new(Substrate::new((), SubstrateType::ConventionalSubstrate)))
}

fn vec_pointer_equality_optional_weak(
    a: &Vec<Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>>,
    b: &Vec<Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>>,
) -> bool {
    if a.len() != b.len() {
        false
    } else {
        a.iter().zip(b.iter()).all(|(aa, bb)| match (aa, bb) {
            (None, None) => true,
            (Some(ap), Some(bp)) => ap.ptr_eq(bp),
            _ => false,
        })
    }
}
