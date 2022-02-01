use crate::evolution::{
    helper::{
        noop::{NoOpInputElement, NoOpInputSensor, NoOpReaction, NoOpState, NoOpSubstrate},
        testing::{TestInformation, TestInput, TestReaction, TestState},
    },
    protein::{CatalyticCentre, SubstrateType},
};

use super::*;

#[test]
/// Tests if the function `new` correctly creates an [`InputSensor`].
fn test_new() {
    let input_substrates = vec![Some(new_noop_substrate()), Some(new_noop_substrate()), None];
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: NoOpInputSensor = ();
    let input_sensor: InputSensor<
        NoOpReaction,
        NoOpState,
        NoOpSubstrate,
        NoOpInputElement,
        NoOpInputSensor,
    > = InputSensor::new(input.clone(), input_substrates_weak.clone());
    assert!(vec_pointer_equality_optional_weak(
        input_sensor.input_substrates(),
        &input_substrates_weak
    ));
    assert_eq!(input_sensor.input(), &input);
}

#[test]
/// Tests if the function `set_input` correctly sets the input.
fn test_set_input() {
    let input_substrates = vec![
        Some(new_test_substrate(10)),
        Some(new_test_substrate(20)),
        None,
    ];
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: TestInput = TestInput {
        last_feedback_update: HashMap::new(),
        last_set_input: (),
        mocked_output_vector: vec![
            TestInformation::default(),
            TestInformation::default(),
            TestInformation::default(),
        ],
    };
    let mut input_sensor: InputSensor<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
    > = InputSensor::new(input.clone(), input_substrates_weak.clone());
    let expected_input_substrate_values: Vec<Option<TestInformation>> = vec![
        Some(TestInformation::default()),
        Some(TestInformation::default()),
        None,
    ];
    input_sensor.set_input(());
    let input_substrate_values: Vec<Option<TestInformation>> =
        substrates_as_owned(input_sensor.input_substrates())
            .into_iter()
            .map(|o| o.map(|s| s.value))
            .collect();
    assert_eq!(expected_input_substrate_values, input_substrate_values);
}

#[test]
#[should_panic]
/// Tests if the function `set_input` correctly panics if the underlying implementation is wrong.
fn test_set_input_length_mismatch() {
    let input_substrates = vec![
        Some(new_test_substrate(1)),
        Some(new_test_substrate(2)),
        None,
    ];
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: TestInput = TestInput {
        last_feedback_update: HashMap::new(),
        last_set_input: (),
        mocked_output_vector: Vec::new(),
    };
    let mut input_sensor: InputSensor<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
    > = InputSensor::new(input.clone(), input_substrates_weak.clone());
    input_sensor.set_input(());
}

#[test]
/// Tests if the function `feedback_update` correctly propagates feedback substrate changes.
fn test_feedback_update() {
    let receptor = Rc::new(Receptor::new(
        Vec::new(),
        TestState {},
        CatalyticCentre::new(Vec::new(), Vec::new(), TestReaction {}),
    ));
    let input_substrates = vec![
        Some(new_test_substrate(0)),
        Some(new_test_substrate(1)),
        None,
    ];
    input_substrates[0]
        .as_ref()
        .unwrap()
        .borrow_mut()
        .add_receptor(receptor.clone());
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: TestInput = TestInput {
        last_feedback_update: HashMap::new(),
        last_set_input: (),
        mocked_output_vector: vec![
            TestInformation::default(),
            TestInformation::default(),
            TestInformation::default(),
        ],
    };
    let mut input_sensor: InputSensor<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
    > = InputSensor::new(input.clone(), input_substrates_weak.clone());
    assert!(input_sensor.feedback_update(HashMap::new()).is_empty());
    assert_eq!(&input_sensor.input().last_feedback_update, &HashMap::new());
    let mut changes: HashMap<usize, TestInformation> = HashMap::new();
    changes.insert(0, TestInformation::default());
    assert_eq!(input_sensor.cascading_receptors(), input_sensor.feedback_update(changes));
    let expected_input_substrate_values: Vec<Option<TestInformation>> = vec![
        Some(TestInformation::default()),
        Some(TestInformation::default()),
        None,
    ];
    let input_substrate_values: Vec<Option<TestInformation>> =
        substrates_as_owned(input_sensor.input_substrates())
            .into_iter()
            .map(|o| o.map(|s| s.value))
            .collect();
    assert_eq!(input_substrate_values, expected_input_substrate_values);
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
    let input_substrates_weak: Vec<
        Option<Weak<RefCell<Substrate<TestReaction, TestState, TestInformation>>>>,
    > = input_substrates
        .iter()
        .map(Option::as_ref)
        .map(|s| s.map(Rc::downgrade))
        .collect();
    let input: TestInput = TestInput {
        last_feedback_update: HashMap::new(),
        last_set_input: (),
        mocked_output_vector: Vec::new(),
    };
    let input_sensor: InputSensor<
        TestReaction,
        TestState,
        TestInformation,
        NoOpInputElement,
        TestInput,
    > = InputSensor::new(input.clone(), input_substrates_weak.clone());
    assert_eq!(
        input_sensor.cascading_receptors(),
        vec![receptors[0].clone(), receptors[1].clone()]
    );
}

fn new_noop_substrate() -> Rc<RefCell<Substrate<NoOpReaction, NoOpState, NoOpSubstrate>>> {
    Rc::new(RefCell::new(Substrate::new((), SubstrateType::ConventionalSubstrate)))
}

fn new_test_substrate(
    value: usize,
) -> Rc<RefCell<Substrate<TestReaction, TestState, TestInformation>>> {
    Rc::new(RefCell::new(Substrate::new(
        TestInformation { value },
        SubstrateType::ConventionalSubstrate,
    )))
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

/// Converts a vector of references to [`Substrate`]s.
///
/// # Parameters
///
/// * `substrates` - the [`Substrate`] references
fn substrates_as_owned<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrates: &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
) -> Vec<Option<Substrate<ReactionType, StateType, InformationType>>> {
    substrates
        .iter()
        .map(|substrate| substrate.as_ref())
        .map(|substrate| substrate.map(substrate_reference_to_owned))
        .collect()
}

/// Converts a reference to a [`Substrate`].
///
/// # Parameters
///
/// * `substrate_reference` - the [`Substrate`] reference
fn substrate_reference_to_owned<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrate_reference: &Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
) -> Substrate<ReactionType, StateType, InformationType> {
    (*(substrate_reference.upgrade().unwrap())).borrow().clone()
}