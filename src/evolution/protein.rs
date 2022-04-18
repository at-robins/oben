//! The `protein` module contains the executive part of the evolutionary network.
extern crate bitvec;
extern crate rmp_serde;

pub use sensor::{InputSensor, OutputSensor};

use super::chemistry::{Information, Reaction, State};
use super::helper::Iteration;
use std::borrow::Borrow;
use std::cell::{RefCell, RefMut};
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
    substrate_type: SubstrateType,
}

impl<
        ReactionType: Reaction<InformationType>,
        StateType: State<InformationType>,
        InformationType: Information,
    > Substrate<ReactionType, StateType, InformationType>
{
    /// Creates a new`Substrate` with the specified value and type
    ///
    /// # Parameters
    ///
    /// * `value` - the value of the `Substrate`
    /// * `substrate_type` - the type of the substrate
    pub fn new(value: InformationType, substrate_type: SubstrateType) -> Self {
        Substrate {
            value,
            receptors: vec![],
            last_updated: Iteration::new(),
            substrate_type,
        }
    }

    /// Returns the type of this `Substrate`.
    pub fn substrate_type(&self) -> &SubstrateType {
        &self.substrate_type
    }

    /// Changes the [`SubstrateType`] as specified.
    ///
    /// # Parameters
    ///
    /// * `substrate_type` - the new [`SubstrateType`]
    pub fn set_substrate_type(&mut self, substrate_type: SubstrateType) {
        self.substrate_type = substrate_type;
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

    /// Returns the input feedback associations an the current value if of type [`SubstrateType::InputFeedbackSubstrate`].
    ///
    /// # Parameters
    ///
    /// * time - the timepoint the substrate is accessed
    pub fn get_input_feedback_associations(
        &mut self,
        time: Iteration,
    ) -> Option<(Vec<usize>, InformationType)> {
        match &self.substrate_type {
            SubstrateType::InputFeedbackSubstrate(feedback_indices) => {
                Some((feedback_indices.clone(), self.value(time).clone()))
            },
            _ => None,
        }
    }

    /// Returns the the current value if of type [`SubstrateType::OutputFinishSubstrate`].
    ///
    /// # Parameters
    ///
    /// * time - the timepoint the substrate is accessed
    pub fn get_output_finish_value(&mut self, time: Iteration) -> Option<InformationType> {
        match &self.substrate_type {
            SubstrateType::OutputFinishSubstrate => Some(self.value(time).clone()),
            _ => None,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
/// The type of [`Substrate`].
/// Special types of [`Substrate`]s provide additional information
/// on unique associations.
pub enum SubstrateType {
    /// A convecntional [`Substrate`] without any special associations.
    ConventionalSubstrate,
    /// A [`Substrate`] that is associated with feedback to an [`InputSensor`].
    InputFeedbackSubstrate(Vec<usize>),
    /// A [`Substrate`] signalling that the current process is finsihed to the [`OutputSensor`].
    OutputFinishSubstrate,
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
    /// and returned together with other associations.
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
    ) -> CatalysisResult<ReactionType, StateType, InformationType> {
        self.enzyme.catalyse(time_of_catalysis);
        CatalysisResult {
            cascading_receptors: self.enzyme.cascading_receptors(),
            input_feedback_associations: self.enzyme.input_feedback_associations(time_of_catalysis),
            output_finish_substrate: self.enzyme.output_finish_substrate(time_of_catalysis),
        }
    }
}

#[derive(Debug, Clone)]
/// The result of a catalysis.
pub struct CatalysisResult<ReactionType, StateType, InformationType> {
    /// The [`Receptor`]s affected by this catalysis.
    pub cascading_receptors: Vec<Rc<Receptor<ReactionType, StateType, InformationType>>>,
    /// The input feedback associations affected by this catalysis.
    pub input_feedback_associations: Vec<(Vec<usize>, InformationType)>,
    pub output_finish_substrate: Option<InformationType>,
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

    /// Returns all input feedback associations and values.
    ///
    /// # Parameters
    ///
    /// * `time` - the timepoint at which the [`Substrate`] values are requested
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    pub fn input_feedback_associations(
        &self,
        time: Iteration,
    ) -> Vec<(Vec<usize>, InformationType)> {
        self.products
            .iter()
            .map(|product| {
                (*product.upgrade().unwrap())
                    .borrow_mut()
                    .get_input_feedback_associations(time)
            })
            .filter(Option::is_some)
            .map(Option::unwrap)
            .collect()
    }

    /// Returns the finish substrate value if the finish output substrate is
    /// altered by this catalytic centre.
    ///
    /// # Parameters
    ///
    /// * `time` - the timepoint at which the [`Substrate`] values are requested
    /// as [`Iteration`](crate::evolution::helper::Iteration)
    pub fn output_finish_substrate(&self, time: Iteration) -> Option<InformationType> {
        for product in &self.products {
            let finish_value = (*product.upgrade().unwrap())
                .borrow_mut()
                .get_output_finish_value(time);
            if finish_value.is_some() {
                return finish_value;
            }
        }
        None
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

mod sensor;
