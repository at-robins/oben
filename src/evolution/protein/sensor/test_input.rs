use crate::evolution::{
    helper::{
        noop::{NoOpInputElement, NoOpInputSensor, NoOpReaction, NoOpState, NoOpSubstrate},
        testing::{TestInformation, TestInput, TestReaction, TestState},
    },
    protein::CatalyticCentre,
};

use super::*;

#[test]
/// Tests if the function `new` correctly creates an [`InputSensor`].
fn test_new() {
    let input_substrates = vec![Some(new_noop_substrate()), Some(new_noop_substrate()), None];
    let feedback_substrates = vec![
        Some(new_noop_substrate()),
        None,
        input_substrates[0].clone(),
    ];
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let feedback_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
    > = feedback_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: NoOpInputSensor = ();
    let mut input_sensor: InputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = InputSensor::new(
        input.clone(),
        input_substrates_weak.clone(),
        feedback_substrates_weak.clone(),
    );
    assert!(vec_pointer_equality_optional_weak(
        input_sensor.input_substrates(),
        &input_substrates_weak
    ));
    assert!(vec_pointer_equality_optional_weak(
        input_sensor.feedback_substrates(),
        &feedback_substrates_weak
    ));
    assert!(!input_sensor.was_changed());
    assert_eq!(input_sensor.input(), &input);
}

#[test]
/// Tests if the function `was_changed` correctly registeres change.
fn test_was_changed() {
    let input_substrates = vec![Some(new_noop_substrate()), Some(new_noop_substrate()), None];
    let feedback_substrates = vec![
        Some(new_noop_substrate()),
        None,
        input_substrates[0].clone(),
    ];
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let feedback_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
    > = feedback_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: NoOpInputSensor = ();
    let mut input_sensor: InputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = InputSensor::new(
        input.clone(),
        input_substrates_weak.clone(),
        feedback_substrates_weak.clone(),
    );
    assert!(!input_sensor.was_changed());
    input_sensor.set_input(());
    assert!(input_sensor.was_changed());
    assert!(!input_sensor.was_changed());
}

#[test]
/// Tests if the function `feedback_update` correctly propagates feedback substrate changes.
fn test_feedback_update() {
    let input_substrates = vec![
        Some(new_test_substrate(0)),
        Some(new_test_substrate(1)),
        None,
    ];
    let mut feedback_substrates = vec![
        Some(new_test_substrate(2)),
        None,
        input_substrates[0].clone(),
    ];
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let feedback_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = feedback_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: TestInput = TestInput {
        last_feedback_update: Vec::new(),
        last_set_input: (),
    };
    let mut input_sensor: InputSensor<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
    > = InputSensor::new(
        input.clone(),
        input_substrates_weak.clone(),
        feedback_substrates_weak.clone(),
    );
    input_sensor.feedback_update();
    assert_eq!(&input_sensor.input().last_feedback_update, &Vec::new());
    assert!(!input_sensor.was_changed());
    feedback_substrates[0]
        .as_mut()
        .unwrap()
        .borrow_mut()
        .set_value(TestInformation { value: 3 });
    input_sensor.feedback_update();
    let expected_changes = vec![Some(TestInformation { value: 3 }), None, None];
    assert_eq!(&input_sensor.input().last_feedback_update, &expected_changes);
    assert!(input_sensor.was_changed());
    assert!(!input_sensor.was_changed());
}

#[test]
/// Tests if the function `cascading_receptors` correctly returns all receptors belonging to the input.
fn test_cascading_receptors() {
    let receptors = vec![
        Rc::new(Receptor::new(
            Vec::new(),
            TestState {},
            CatalyticCentre::new(Vec::new(), Vec::new(), TestReaction {}),
        )),
        Rc::new(Receptor::new(
            Vec::new(),
            TestState {},
            CatalyticCentre::new(Vec::new(), Vec::new(), TestReaction {}),
        )),
        Rc::new(Receptor::new(
            Vec::new(),
            TestState {},
            CatalyticCentre::new(Vec::new(), Vec::new(), TestReaction {}),
        )),
    ];
    let input_substrates = vec![
        Some(new_test_substrate(0)),
        Some(new_test_substrate(1)),
        None,
    ];
    input_substrates[0]
        .as_ref()
        .unwrap()
        .borrow_mut()
        .add_receptor(receptors[0].clone());
    input_substrates[1]
        .as_ref()
        .unwrap()
        .borrow_mut()
        .add_receptor(receptors[1].clone());
    let feedback_substrates = vec![
        Some(new_test_substrate(2)),
        None,
        input_substrates[0].clone(),
    ];
    feedback_substrates[0]
        .as_ref()
        .unwrap()
        .borrow_mut()
        .add_receptor(receptors[2].clone());
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let feedback_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = feedback_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: TestInput = TestInput {
        last_feedback_update: Vec::new(),
        last_set_input: (),
    };
    let input_sensor: InputSensor<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
    > = InputSensor::new(
        input.clone(),
        input_substrates_weak.clone(),
        feedback_substrates_weak.clone(),
    );
    assert_eq!(
        input_sensor.cascading_receptors(),
        vec![receptors[0].clone(), receptors[1].clone()]
    );
}

fn new_noop_substrate() -> Rc<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>> {
    Rc::new(RefCell::new(Substrate::new(())))
}

fn new_test_substrate(
    value: usize,
) -> Rc<RefCell<Substrate<TestReaction, TestState, TestInformation>>> {
    Rc::new(RefCell::new(Substrate::new(TestInformation { value })))
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
