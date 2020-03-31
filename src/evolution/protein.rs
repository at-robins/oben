//! The `protein` module contains the executive part of the evolutionary network.
extern crate bitvec;

use std::cell::Ref;
use core::cell::RefCell;
use std::rc::Rc;
use bitvec::boxed::BitBox;
use super::chemistry::{Reaction, State};

/// A `Substrate` represents a chemical entity of a specific value. Additionally
/// a `Substrate` is aware of all [`Receptor`]s detecting its changes.
///
/// [`Receptor`]: ./struct.Receptor.html
#[derive(Clone)]
pub struct Substrate {
    value: BitBox,
    receptors: Vec<Rc<Receptor>>,
}

impl Substrate {
    /// Creates a new `Substrate` with the specified binary value.
    pub fn new(value: BitBox) -> Self {
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
    pub fn set_value(&mut self, value: BitBox) {
        self.value = value;
    }

    /// Returns the binary value of this substrate.
    pub fn value(&self) -> &BitBox {
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
#[derive(Clone)]
pub struct Receptor {
    substrates: Vec<Rc<RefCell<Substrate>>>,
    state: State,
    enzyme: CatalyticCentre,
}

impl Receptor {
    /// Creates a `Receptor` detecting the specified [`State`] of its substrates and triggering
    /// the [`CatalyticCentre`]'s reaction if appropriate.
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn new(state: State, enzyme: CatalyticCentre) -> Self {
        Receptor{substrates: vec!(), state, enzyme}
    }

    /// Detects the [`State`] of its substrates and determines if triggering the
    /// [`CatalyticCentre`]'s reaction is appropriate.
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    fn should_trigger(&self) -> bool {
        // TODO: insert check of substrate number
        // TODO: refactor this ugly code
        let substrates: Vec<Ref<Substrate>> = self.substrates.iter()
            .map(|sub| sub.borrow())
            .collect();
        let substrates: Vec<&BitBox> = substrates.iter()
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

    /// Adds a [`Substrate`] to the substrate that should be notified upon change.
    ///
    /// [`Substrate`]: ./struct.Substrate.html
    pub fn add_substrate(&mut self, substrate: Rc<RefCell<Substrate>>) {
        self.substrates.push(substrate);
    }
}

/// A `CatalyticCentre` produces products from educt [`Substrate`]s
/// by performing a [`Reaction`].
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Clone)]
pub struct CatalyticCentre {
    // TODO: implement `new` and remove public modifier
    pub educts: Vec<Rc<RefCell<Substrate>>>,
    pub products: Vec<Rc<RefCell<Substrate>>>,
    pub reaction: Reaction,
}

impl CatalyticCentre {
    ///
    pub fn new(reaction: Reaction) -> Self {
        CatalyticCentre{educts: vec!(), products: vec!(), reaction}
    }

    /// Calculates the values of the products after performing the
    /// [`Reaction`] specific for this catalytic centre.
    ///
    /// [`Reaction`]: ../chemistry/struct.Reaction.html
    fn calculate_product_values(&self) -> Vec<BitBox> {
        // TODO: insert check of product and educt number
        // TODO: refactor this ugly code
        let educts: Vec<Ref<Substrate>> = self.educts.iter()
            .map(|sub| sub.borrow())
            .collect();
        let educts: Vec<&BitBox> = educts.iter()
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
