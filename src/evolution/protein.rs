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
    pub value: BitBox,
    pub receptors: Vec<Rc<Receptor>>,
}

impl Substrate {
    /*pub fn new(value: BitBox) -> Self {

    }*/

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
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Reaction.html
#[derive(Clone)]
pub struct Receptor {
    substrates: Vec<Rc<Substrate>>,
    state: State,
    enzym: CatalyticCentre,
}

impl Receptor {
    /// Detects the [`State`] of its substrates and triggers a
    /// [`CatalyticCentre`]'s reaction if appropriate. Subsequent ("cascading")
    /// receptors, which are supposed to be checked after the reaction was triggered,
    /// are returned.
    ///
    /// [`State`]: ../chemistry/struct.State.html
    /// [`CatalyticCentre`]: ./struct.CatalyticCentre.html
    pub fn detect(&self) -> Vec<Rc<Receptor>> {
        // TODO: insert check of substrate number
        let substrates: Vec<&BitBox> = self.substrates.iter()
            .map(|sub| sub.value())
            .collect();
        if self.state.detect(&substrates) {
            self.enzym.catalyse();
            self.enzym.cascading_receptors()
        } else {
            vec!()
        }
    }
}

/// A `CatalyticCentre` produces products from educt [`Substrate`]s
/// by performing a [`Reaction`].
///
/// [`Substrate`]: ./struct.Substrate.html
/// [`Reaction`]: ../chemistry/struct.Substrate.html
#[derive(Clone)]
pub struct CatalyticCentre {
    pub educts: Vec<Rc<RefCell<Substrate>>>,
    pub products: Vec<Rc<RefCell<Substrate>>>,
    pub reaction: Reaction,
}

impl CatalyticCentre {
    /// Catalyses the [`Reaction`] specific for this catalytic centre.
    ///
    /// [`Reaction`]: ../chemistry/struct.Substrate.html
    pub fn catalyse(&self) {
        // TODO: insert check of product and educt number
        let mut product_values;
        {
            // TODO: refactor this ugly code
            let educts: Vec<Ref<Substrate>> = self.educts.iter()
                .map(|sub| sub.borrow())
                .collect();
            let educts: Vec<&BitBox> = educts.iter()
                .map(|sub| sub.value())
                .collect();
            product_values = self.reaction.react(&educts);
        }
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
