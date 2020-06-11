//! The `protein` module contains the executive part of the evolutionary network.
extern crate bitvec;

use std::cell::Ref;
use core::cell::RefCell;
use std::rc::Rc;
use bitvec::{boxed::BitBox, order::Local};
use super::chemistry::{Reaction, State};

/// A `Substrate` represents a chemical entity of a specific value. Additionally
/// a `Substrate` is aware of all [`Receptor`]s detecting its changes.
///
/// [`Receptor`]: ./struct.Receptor.html
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Substrate {
    value: BitBox<Local, u8>,
    receptors: Vec<Rc<Receptor>>,
}

impl Substrate {
    /// Creates a new `Substrate` with the specified binary value.
    pub fn new(value: BitBox<Local, u8>) -> Self {
            Substrate{value, receptors: vec!()}
    }

    /// Set the value of this substrate.
    /// This method will not notify the corresponding [`Receptor`]s.
    ///
    /// # Parameters
    ///
    /// * `value` - the new value of the substrate
    ///
    /// [`Receptor`]: ./struct.Receptor.html
    pub fn set_value(&mut self, value: BitBox<Local, u8>) {
        self.value = value;
    }

    /// Returns the binary value of this substrate.
    pub fn value(&self) -> &BitBox<Local, u8> {
        &self.value
    }

    /// Returns all receptors detecting this substrate.
    pub fn receptors(&self) -> Vec<Rc<Receptor>> {
        self.receptors.iter().map(Rc::clone).collect()
    }

    /// Adds a [`Receptor`] to the substrate that should be notified upon change.
    ///
    /// [`Receptor`]: ./struct.Receptor.html
    pub fn add_receptor(&mut self, receptor: Rc<Receptor>) {
        self.receptors.push(receptor);
    }
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Receptor {
    substrates: Vec<Rc<RefCell<Substrate>>>,
    state: State,
    enzyme: CatalyticCentre,
}

impl Receptor {
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
    pub fn new(substrates: Vec<Rc<RefCell<Substrate>>>, state: State, enzyme: CatalyticCentre) -> Self {
        assert_eq!(substrates.len(), state.get_substrate_number(),
            "The number of required substrates to check for state {:?} is {}, but {} substrates were supplied.",
            state, state.get_substrate_number(), substrates.len());
        Receptor{substrates, state, enzyme}
    }

    /// Detects the [`State`] of its substrates and determines if triggering the
    /// [`CatalyticCentre`]'s reaction is appropriate.
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    fn should_trigger(&self) -> bool {
        // TODO: refactor this ugly code
        let substrates: Vec<Ref<Substrate>> = self.substrates.iter()
            .map(|sub| sub.borrow())
            .collect();
        let substrates: Vec<&BitBox<Local, u8>> = substrates.iter()
            .map(|sub| sub.value())
            .collect();
        self.state.detect(&substrates)
    }

    /// Detects the [`State`] of its substrates and triggers a
    /// [`CatalyticCentre`]'s reaction if appropriate. Subsequent ("cascading")
    /// receptors, which are supposed to be checked after the reaction was triggered,
    /// are returned.
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn detect(&self) -> Vec<Rc<Receptor>> {
        if self.should_trigger() {
            self.enzyme.catalyse();
            self.enzyme.cascading_receptors()
        } else {
            vec!()
        }
    }
}

/// A `CatalyticCentre` produces products from educt [`Substrate`]s
/// by performing a [`Reaction`].
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CatalyticCentre {
    educts: Vec<Rc<RefCell<Substrate>>>,
    products: Vec<Rc<RefCell<Substrate>>>,
    reaction: Reaction,
}

impl CatalyticCentre {
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
    pub fn new(educts: Vec<Rc<RefCell<Substrate>>>, products: Vec<Rc<RefCell<Substrate>>>, reaction: Reaction) -> Self {
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
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    fn calculate_product_values(&self) -> Vec<BitBox<Local, u8>> {
        // TODO: refactor this ugly code
        let educts: Vec<Ref<Substrate>> = self.educts.iter()
            .map(|sub| sub.borrow())
            .collect();
        let educts: Vec<&BitBox<Local, u8>> = educts.iter()
            .map(|sub| sub.value())
            .collect();
        self.reaction.react(&educts)
    }

    /// Catalyses the [`Reaction`] specific for this catalytic centre.
    ///
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    pub fn catalyse(&self) {
        let mut product_values = self.calculate_product_values();
        for product in &self.products {
            // TODO: maybe switch to VecDeque and use pop_first()
            product.borrow_mut().set_value(product_values.remove(0));
        }
    }

    /// Returns all receptors that detect the product [`Substrate`]s
    /// of the catalytic centre.
    ///
    /// [`Substrate`]: ./struct.Substrate.html
    pub fn cascading_receptors(&self) -> Vec<Rc<Receptor>> {
        self.products.iter().flat_map(|product| product.borrow().receptors()).collect()
    }
}
