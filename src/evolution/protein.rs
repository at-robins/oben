//! The `protein` module contains the executive part of the evolutionary network.
extern crate bitvec;
extern crate rmp_serde;

use std::cell::{RefCell, RefMut};
use std::rc::{Rc, Weak};
use super::chemistry::{Information, Reaction, State};
use super::helper::Iteration;

/// A `Substrate` represents a chemical entity of a specific value. Additionally
/// a `Substrate` is aware of all [`Receptor`]s detecting its changes.
///
/// [`Receptor`]: ./struct.Receptor.html
#[derive(Debug, Clone)]
pub struct Substrate<R, S, T> {
    value: T,
    receptors: Vec<Rc<Receptor<R, S, T>>>,
    last_updated: Iteration,
}

impl<R: Reaction<T>, S: State<T>, T: Information> Substrate<R, S, T> {
    /// Creates a new `Substrate` with the specified binary value.
    pub fn new(value: T) -> Self {
        Substrate{value, receptors: vec!(), last_updated: Iteration::new()}
    }

    /// Set the value of this substrate.
    /// This method will not notify the corresponding [`Receptor`]s.
    ///
    /// # Parameters
    ///
    /// * `value` - the new value of the substrate
    ///
    /// [`Receptor`]: ./struct.Receptor.html
    pub fn set_value(&mut self, value: T) {
        self.value = value;
    }

    /// Returns the value of this substrate at the specified timepoint.
    ///
    /// # Parameters
    ///
    /// * time - the timepoint the substrate is accessed
    pub fn value(&mut self, time: Iteration) -> &T {
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
        rmp_serde::to_vec(&self.value).expect("Serialisation of the genome failed.").len() * 8
    }

    /// Returns all receptors detecting this substrate.
    pub fn receptors(&self) -> Vec<Rc<Receptor<R, S, T>>> {
        self.receptors.iter().map(Rc::clone).collect()
    }

    /// Adds a [`Receptor`] to the substrate that should be notified upon change.
    ///
    /// [`Receptor`]: ./struct.Receptor.html
    pub fn add_receptor(&mut self, receptor: Rc<Receptor<R, S, T>>) {
        self.receptors.push(receptor);
    }
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, Clone)]
pub struct Receptor<R, S, T> {
    substrates: Vec<Weak<RefCell<Substrate<R, S, T>>>>,
    state: S,
    enzyme: CatalyticCentre<R, S, T>,
}

impl<R: Reaction<T>, S: State<T>, T: Information> Receptor<R, S, T> {
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
    pub fn new(substrates: Vec<Weak<RefCell<Substrate<R, S, T>>>>, state: S, enzyme: CatalyticCentre<R, S, T>) -> Self {
        assert_eq!(substrates.len(), state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            state, state.get_substrate_number(), substrates.len());
        Receptor{substrates, state, enzyme}
    }

    /// Detects the [`State`] of its substrates and determines if triggering the
    /// [`CatalyticCentre`]'s reaction is appropriate.
    ///
    /// # Parameters
    ///
    /// * `time_of_detection` - the timepoint at which the catalysis happens
    /// as [`Iteration`](oben::evolution::helper::Iteration)
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn detect(&self, time_of_detection: Iteration) ->  bool {
        // TODO: refactor this ugly code
        let strong: Vec<Rc<RefCell<Substrate<R, S, T>>>> = self.substrates.iter()
            // This unwrap must succeed as the containing structure will always be dropped first
            // and no substrate references are leaked.
            .map(|weak| weak.upgrade().unwrap())
            .collect();
        let mut substrates: Vec<RefMut<Substrate<R, S, T>>> = strong.iter()
            .map(|sub| sub.borrow_mut())
            .collect();
        let substrate_values: Vec<&T> = substrates.iter_mut()
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
    /// as [`Iteration`](oben::evolution::helper::Iteration)
    ///
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn catalyse(&self, time_of_catalysis: Iteration) -> Vec<Rc<Receptor<R, S, T>>> {
        self.enzyme.catalyse(time_of_catalysis);
        self.enzyme.cascading_receptors()
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> PartialEq for Receptor<R, S, T>  {
    fn eq(&self, other: &Self) -> bool {
        if self.state != other.state { return false; }
        if self.enzyme != other.enzyme { return false; }
        if self.substrates.len() != other.substrates.len() { return false; }
        // Compare the raw pointers.
        self.substrates.iter()
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
pub struct CatalyticCentre<R, S, T> {
    educts: Vec<Weak<RefCell<Substrate<R, S, T>>>>,
    products: Vec<Weak<RefCell<Substrate<R, S, T>>>>,
    reaction: R,
}

impl<R: Reaction<T>, S: State<T>, T: Information> CatalyticCentre<R, S, T> {
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
    pub fn new(educts: Vec<Weak<RefCell<Substrate<R, S, T>>>>, products: Vec<Weak<RefCell<Substrate<R, S, T>>>>, reaction: R) -> Self {
        assert_eq!(educts.len(), reaction.get_educt_number(),
            "The number of required educts for reaction {:?} is {}, but {} educts were supplied.",
            reaction, reaction.get_educt_number(), educts.len());
        assert_eq!(products.len(), reaction.get_product_number(),
            "The number of required products for reaction {:?} is {}, but {} products were supplied.",
            reaction, reaction.get_product_number(), products.len());
        CatalyticCentre{educts, products, reaction}
    }

    /// Calculates the values of the products after performing the
    /// [`Reaction`] specific for this catalytic centre.
    ///
    /// # Parameters
    ///
    /// * `time_of_catalysis` - the timepoint at which the catalysis happens
    /// as [`Iteration`](oben::evolution::helper::Iteration)
    ///
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    fn calculate_product_values(&self, time_of_catalysis: Iteration) -> Vec<T> {
        // TODO: refactor this ugly code
        let strong: Vec<Rc<RefCell<Substrate<R, S, T>>>> = self.educts.iter()
            // This unwrap must succeed as the containing structure will always be dropped first
            // and no substrate references are leaked.
            .map(|weak| weak.upgrade().unwrap())
            .collect();
        let mut educts: Vec<RefMut<Substrate<R, S, T>>> = strong.iter()
            .map(|sub| sub.borrow_mut())
            .collect();
        let educts: Vec<&T> = educts.iter_mut()
            .map(|sub| sub.value(time_of_catalysis))
            .collect();
        self.reaction.react(&educts, time_of_catalysis)
    }

    /// Catalyses the [`Reaction`] specific for this catalytic centre.
    ///
    /// # Parameters
    ///
    /// * `time_of_catalysis` - the timepoint at which the catalysis happens
    /// as [`Iteration`](oben::evolution::helper::Iteration)
    ///
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    pub fn catalyse(&self, time_of_catalysis: Iteration) {
        let mut product_values = self.calculate_product_values(time_of_catalysis);
        for product in &self.products {
            // TODO: maybe switch to VecDeque and use pop_first()
            // This unwrap must succeed as the containing structure will always be dropped first
            // and no substrate references are leaked.
            product.upgrade()
                .unwrap()
                .borrow_mut()
                .set_value(product_values.remove(0));
        }
    }

    /// Returns all receptors that detect the product [`Substrate`]s
    /// of the catalytic centre.
    ///
    /// [`Substrate`]: ./struct.Substrate.html
    pub fn cascading_receptors(&self) -> Vec<Rc<Receptor<R, S, T>>> {
        // This unwrap must succeed as the containing structure will always be dropped first
        // and no substrate references are leaked.
        self.products.iter().flat_map(substrate_reference_to_receptors).collect()
    }
}

impl<R: Reaction<T>, S: State<T>, T: Information> PartialEq for CatalyticCentre<R, S, T>  {
    fn eq(&self, other: &Self) -> bool {
        if self.reaction != other.reaction { return false; }
        if self.educts.len() != other.educts.len() { return false; }
        if self.products.len() != other.products.len() { return false; }
        // Compare the raw pointers.
        let educt_pairs_not_matching = self.educts.iter()
            .zip(other.educts.iter())
            .any(|(a, b)| !a.ptr_eq(b));
        if educt_pairs_not_matching { return false; }
        self.products.iter()
            .zip(other.products.iter())
            .all(|(a, b)| a.ptr_eq(b))
    }
}

/// Returns the receptors of a weak substrate reference.
/// 
/// # Parameters
/// 
/// * `substrate` - the substrate reference
fn substrate_reference_to_receptors<R: Reaction<T>, S: State<T>, T: Information>(substrate: &Weak<RefCell<Substrate<R, S, T>>>) -> Vec<Rc<Receptor<R, S, T>>> {
    let upgrade: Rc<RefCell<Substrate<R, S, T>>> = substrate.upgrade().unwrap();
    let s: &RefCell<Substrate<R, S, T>> = upgrade.borrow();
    return s.borrow().receptors();
}