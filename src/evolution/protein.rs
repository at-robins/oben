//! The `protein` module contains the executive part of the evolutionary network.
extern crate bitvec;
extern crate rmp_serde;

use serde::de::DeserializeOwned;
use serde::Serialize;

use super::chemistry::{Information, Input, Reaction, State};
use super::gene::CrossOver;
use super::helper::Iteration;
use std::borrow::Borrow;
use std::cell::{RefCell, RefMut};
use std::marker::PhantomData;
use std::rc::{Rc, Weak};

/// A `Substrate` represents a chemical entity of a specific value. Additionally
/// a `Substrate` is aware of all [`Receptor`]s detecting its changes.
///
/// [`Receptor`]: ./struct.Receptor.html
#[derive(Debug, Clone)]
pub struct Substrate<ReactionType, StateType, InformationType> {
    value: InformationType,
    receptors: Vec<Rc<Receptor<ReactionType, StateType, InformationType>>>,
    last_updated: Iteration,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
    > Substrate<ReactionType, StateType, InformationType>
{
    /// Creates a new `Substrate` with the specified binary value.
    pub fn new(value: InformationType) -> Self {
        Substrate {
            value,
            receptors: vec![],
            last_updated: Iteration::new(),
        }
    }

    /// Set the value of this substrate.
    /// This method will not notify the corresponding [`Receptor`]s.
    ///
    /// # Parameters
    ///
    /// * `value` - the new value of the substrate
    ///
    /// [`Receptor`]: ./struct.Receptor.html
    pub fn set_value(&mut self, value: InformationType) {
        self.value = value;
    }

    /// Returns the value of this substrate at the specified timepoint.
    ///
    /// # Parameters
    ///
    /// * time - the timepoint the substrate is accessed
    pub fn value(&mut self, time: Iteration) -> &InformationType {
        if time != self.last_updated {
            let time_passed = time - self.last_updated;
            self.value.update_value(time_passed);
            self.last_updated = time;
        }
        &self.value
    }

    /// Returns the number of bits encoded by this `Substrate`.
    pub fn binary_size(&self) -> usize {
        // TODO: Revise this implementation as it is potentially very expensive.
        rmp_serde::to_vec(&self.value)
            .expect("Serialisation of the genome failed.")
            .len()
            * 8
    }

    /// Returns all receptors detecting this substrate.
    pub fn receptors(&self) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        self.receptors.iter().map(Rc::clone).collect()
    }

    /// Adds a [`Receptor`] to the substrate that should be notified upon change.
    ///
    /// [`Receptor`]: ./struct.Receptor.html
    pub fn add_receptor(
        &mut self,
        receptor: Rc<Receptor<ReactionType, StateType, InformationType>>,
    ) {
        self.receptors.push(receptor);
    }
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, Clone)]
pub struct Receptor<ReactionType, StateType, InformationType> {
    substrates: Vec<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
    state: StateType,
    enzyme: CatalyticCentre<ReactionType, StateType, InformationType>,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
    > Receptor<ReactionType, StateType, InformationType>
{
    /// Creates a `Receptor` detecting the specified [`State`] of its substrates and triggering
    /// the [`CatalyticCentre`]'s reaction if appropriate.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the [`Substrate`]s the [`State`] should check
    /// * `state` - the [`State`] to check for
    /// * `enzyme` - the [`CatalyticCentre`] to trigger if the [`State`] is appropriate
    ///
    /// # Panics
    ///
    /// If the number of [`Substrate`]s is not exactly equal to the one
    /// required by the [`State`].
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn new(
        substrates: Vec<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        state: StateType,
        enzyme: CatalyticCentre<ReactionType, StateType, InformationType>,
    ) -> Self {
        assert_eq!(substrates.len(), state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            state, state.get_substrate_number(), substrates.len());
        Receptor {
            substrates,
            state,
            enzyme,
        }
    }

    /// Detects the [`State`] of its substrates and determines if triggering the
    /// [`CatalyticCentre`]'s reaction is appropriate.
    ///
    /// # Parameters
    ///
    /// * `time_of_detection` - the timepoint at which the catalysis happens
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn detect(&self, time_of_detection: Iteration) -> bool {
        // TODO: refactor this ugly code
        let strong: Vec<Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>> = self.substrates.iter()
            // This unwrap must succeed as the containing structure will always be dropped first
            // and no substrate references are leaked.
            .map(|weak| weak.upgrade().unwrap())
            .collect();
        let mut substrates: Vec<RefMut<Substrate<ReactionType, StateType, InformationType>>> =
            strong.iter().map(|sub| sub.borrow_mut()).collect();
        let substrate_values: Vec<&InformationType> = substrates
            .iter_mut()
            .map(|sub| sub.value(time_of_detection))
            .collect();
        self.state.detect(&substrate_values, time_of_detection)
    }

    /// Triggers a [`CatalyticCentre`]'s reaction. Subsequent ("cascading")
    /// receptors, which are supposed to be checked after the reaction was triggered,
    /// are returned.
    ///
    /// # Parameters
    ///
    /// * `time_of_catalysis` - the timepoint at which the catalysis happens
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    ///
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn catalyse(
        &self,
        time_of_catalysis: Iteration,
    ) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        self.enzyme.catalyse(time_of_catalysis);
        self.enzyme.cascading_receptors()
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> PartialEq for Receptor<R, S, T> {
    fn eq(&self, other: &Self) -> bool {
        if self.state != other.state {
            return false;
        }
        if self.enzyme != other.enzyme {
            return false;
        }
        if self.substrates.len() != other.substrates.len() {
            return false;
        }
        // Compare the raw pointers.
        self.substrates
            .iter()
            .zip(other.substrates.iter())
            .all(|(a, b)| a.ptr_eq(b))
    }
}

/// A `CatalyticCentre` produces products from educt [`Substrate`]s
/// by performing a [`Reaction`].
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, Clone)]
pub struct CatalyticCentre<ReactionType, StateType, InformationType> {
    educts: Vec<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
    products: Vec<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
    reaction: ReactionType,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
    > CatalyticCentre<ReactionType, StateType, InformationType>
{
    /// Creates a new `CatalyticCentre` producing products from educt [`Substrate`]s
    /// by performing a [`Reaction`].
    ///
    /// # Parameters
    ///
    /// * `educts` - the educt [`Substrate`]s for the [`Reaction`]
    /// * `products` - the product [`Substrate`]s for the [`Reaction`]
    /// * `reaction` - the [`Reaction`] to catalyse
    ///
    /// # Panics
    ///
    /// If the number of educt or product [`Substrate`]s is not exactly equal to the one
    /// required by the [`Reaction`].
    ///
    /// [`Substrate`]: ./struct.Substrate.html
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    pub fn new(
        educts: Vec<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        products: Vec<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        reaction: ReactionType,
    ) -> Self {
        assert_eq!(
            educts.len(),
            reaction.get_educt_number(),
            "The number of required educts for reaction {:?} is {}, but {} educts were supplied.",
            reaction,
            reaction.get_educt_number(),
            educts.len()
        );
        assert_eq!(products.len(), reaction.get_product_number(),
            "The number of required products for reaction {:?} is {}, but {} products were supplied.",
            reaction, reaction.get_product_number(), products.len());
        CatalyticCentre {
            educts,
            products,
            reaction,
        }
    }

    /// Calculates the values of the products after performing the
    /// [`Reaction`] specific for this catalytic centre.
    ///
    /// # Parameters
    ///
    /// * `time_of_catalysis` - the timepoint at which the catalysis happens
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    ///
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    fn calculate_product_values(&self, time_of_catalysis: Iteration) -> Vec<InformationType> {
        // TODO: refactor this ugly code
        let strong: Vec<Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>>> = self.educts.iter()
            // This unwrap must succeed as the containing structure will always be dropped first
            // and no substrate references are leaked.
            .map(|weak| weak.upgrade().unwrap())
            .collect();
        let mut educts: Vec<RefMut<Substrate<ReactionType, StateType, InformationType>>> =
            strong.iter().map(|sub| sub.borrow_mut()).collect();
        let educts: Vec<&InformationType> = educts
            .iter_mut()
            .map(|sub| sub.value(time_of_catalysis))
            .collect();
        self.reaction.react(&educts, time_of_catalysis)
    }

    /// Catalyses the [`Reaction`] specific for this catalytic centre.
    ///
    /// # Parameters
    ///
    /// * `time_of_catalysis` - the timepoint at which the catalysis happens
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    ///
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    pub fn catalyse(&self, time_of_catalysis: Iteration) {
        let mut product_values = self.calculate_product_values(time_of_catalysis);
        for product in &self.products {
            // TODO: maybe switch to VecDeque and use pop_first()
            // This unwrap must succeed as the containing structure will always be dropped first
            // and no substrate references are leaked.
            product
                .upgrade()
                .unwrap()
                .borrow_mut()
                .set_value(product_values.remove(0));
        }
    }

    /// Returns all receptors that detect the product [`Substrate`]s
    /// of the catalytic centre.
    ///
    /// [`Substrate`]: ./struct.Substrate.html
    pub fn cascading_receptors(
        &self,
    ) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        // This unwrap must succeed as the containing structure will always be dropped first
        // and no substrate references are leaked.
        self.products
            .iter()
            .flat_map(substrate_reference_to_receptors)
            .collect()
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> PartialEq for CatalyticCentre<R, S, T> {
    fn eq(&self, other: &Self) -> bool {
        if self.reaction != other.reaction {
            return false;
        }
        if self.educts.len() != other.educts.len() {
            return false;
        }
        if self.products.len() != other.products.len() {
            return false;
        }
        // Compare the raw pointers.
        let educt_pairs_not_matching = self
            .educts
            .iter()
            .zip(other.educts.iter())
            .any(|(a, b)| !a.ptr_eq(b));
        if educt_pairs_not_matching {
            return false;
        }
        self.products
            .iter()
            .zip(other.products.iter())
            .all(|(a, b)| a.ptr_eq(b))
    }
}

/// An `InputSensor` registers and transforms inputs into internal [`Substrate`] information 
/// for further processing. It can also operate on the input via feedback [`Substrate`]s. 
#[derive(Debug, Clone)]
pub struct InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    phantom_i: PhantomData<InputElementType>,
    input: InputSensorType,
    information_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
    feedback_substrates:
        Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>>,
    old_feedback_substrate_values: Vec<Option<Substrate<ReactionType, StateType, InformationType>>>,
    changed: bool,
}

impl<
        InputElementType: Clone
            + std::fmt::Debug
            + PartialEq
            + Send
            + Sync
            + CrossOver
            + Serialize
            + DeserializeOwned,
        InputSensorType: Input<InputElementType, InformationType>,
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
    > InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
{
    /// Creates a new `InputSensor` to process input and transform it into internal
    /// [`Substrate`] information.
    /// 
    /// # Parameters
    /// 
    /// * `input` - the genetic input sensor definition
    /// * `information_substrates` - the information processing [`Substrate`]s
    /// * `feedback_substrates` - the [`Substrate`]s used to react to input information flow
    pub fn new(
        input: InputSensorType,
        information_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
        feedback_substrates: Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> InputSensor<ReactionType, StateType, InformationType, InputElementType, InputSensorType>
    {
        let old_feedback_substrate_values: Vec<
            Option<Substrate<ReactionType, StateType, InformationType>>,
        > = Self::substrates_as_owned(&feedback_substrates);
        InputSensor {
            phantom_i: PhantomData,
            input,
            information_substrates,
            feedback_substrates,
            old_feedback_substrate_values,
            changed: false,
        }
    }

    /// Returns the input substrates.
    pub fn input_substrates(&self) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.information_substrates
    }

    /// Returns the feedback substrates.
    pub fn feedback_substrates(&self) -> &Vec<Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>> {
        &self.feedback_substrates
    }

    /// Set the input currently registered by the `InputSensor`.
    /// 
    /// # Parameters
    /// 
    /// * `input` - the current input
    pub fn set_input(&mut self, input: InputElementType) {
        self.changed = true;
        self.input.set_input(input);
    }

    /// Checks if the input was changed since the last query.
    /// The change flag is reset upon calling this method.
    pub fn was_changed(&mut self) -> bool {
        if self.changed {
            self.changed = false;
            true
        } else {
            false
        }
    }

    /// Checks all feedback [`Substrate`]s for changens and passes the changes on to the underlying input sensor representation.
    pub fn feedback_update(&mut self) {
        let current_feedback_values = self.feedback_substrates_as_owned();
        let changes: Vec<Option<InformationType>> = self
            .old_feedback_substrate_values
            .iter()
            .map(Option::as_ref)
            .zip(current_feedback_values.iter().map(Option::as_ref))
            .map(|(old_option, current_option)| match (old_option, current_option) {
                (Some(old), Some(current)) if old.value != current.value => {
                    Some(current.value.clone())
                }
                (None, Some(current)) => Some(current.value.clone()),
                _ => None,
            })
            .collect();

        if changes.iter().any(|o| o.is_some()) {
            self.changed = self.input.handle_feedback_substrate_changes(changes);
            self.old_feedback_substrate_values = current_feedback_values;
        }
    }

    /// Returns all [`Receptor`]s affected by input changes.
    pub fn cascading_receptors(
        &self,
    ) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
        // This unwrap must succeed as the containing structure will always be dropped first
        // and no substrate references are leaked.
        self.information_substrates
            .iter()
            .filter(|i| i.is_some())
            .map(|i| i.as_ref())
            .map(|i| i.unwrap())
            .flat_map(substrate_reference_to_receptors)
            .collect()
    }

    /// Returns the feedback [`Substrate`]s as owned structures.
    fn feedback_substrates_as_owned(
        &self,
    ) -> Vec<Option<Substrate<ReactionType, StateType, InformationType>>> {
        InputSensor::<ReactionType, StateType, InformationType, InputElementType, InputSensorType>::substrates_as_owned(&self.feedback_substrates)
    }

    /// Converts a vector of references to [`Substrate`]s.
    /// 
    /// # Parameters
    /// 
    /// * `substrates` - the [`Substrate`] references
    fn substrates_as_owned(
        substrates: &Vec<
            Option<Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>>,
        >,
    ) -> Vec<Option<Substrate<ReactionType, StateType, InformationType>>> {
        substrates
            .iter()
            .map(|substrate| substrate.as_ref())
            .map(|substrate| substrate.map(Self::substrate_reference_to_owned))
            .collect()
    }

    /// Converts a reference to a [`Substrate`].
    /// 
    /// # Parameters
    /// 
    /// * `substrate_reference` - the [`Substrate`] reference
    fn substrate_reference_to_owned(
        substrate_reference: &Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
    ) -> Substrate<ReactionType, StateType, InformationType> {
        (*(substrate_reference.upgrade().unwrap())).borrow().clone()
    }
}

/// Returns the receptors of a weak substrate reference.
///
/// # Parameters
///
/// * `substrate` - the substrate reference
fn substrate_reference_to_receptors<
    ReactionType: Reaction<InformationType>,
    StateType: State<InformationType>,
    InformationType: Information,
>(
    substrate: &Weak<RefCell<Substrate<ReactionType, StateType, InformationType>>>,
) -> Vec<Rc<Receptor<ReactionType, StateType, InformationType>>> {
    let upgrade: Rc<RefCell<Substrate<ReactionType, StateType, InformationType>>> =
        substrate.upgrade().unwrap();
    let s: &RefCell<Substrate<ReactionType, StateType, InformationType>> = upgrade.borrow();
    return s.borrow().receptors();
}
