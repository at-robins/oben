//! The `protein` module contains the executive part of the evolutionary network.
extern crate bitvec;

use core::cell::RefCell;
use std::rc::Rc;
use bitvec::boxed::BitBox;
use super::chemistry::{Reaction, State};

/// A `Substrate` represents a chemical entity of a specific value. Upon change
/// [`Receptor`]s of the substrate will be notified.
/// 
/// [`Receptor`]: ./struct.Receptor.html
pub struct Substrate {
    value: BitBox,
    receptors: Vec<Rc<Receptor>>,
}

impl Substrate {
    /// Set the value of this substrate and notify all corresponding receptors.
    /// 
    /// # Parameters
    /// 
    /// * `value` - the new value of the substrate
    pub fn set_value(&mut self, value: BitBox) {
        self.value = value;
        self.notify_receptors();
    }
    
    /// Returns the binary value of this substrate.
    pub fn value(&self) -> &BitBox {
        &self.value
    }
    
    /// Causes all receptors detecting this substrate to reevaluate.
    pub fn notify_receptors(&self) {
        for receptor in &self.receptors {
            receptor.detect();
        }
    }
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
/// 
/// [`Substrate`]: ./struct.Substrate.html 
/// [`Reaction`]: ../chemistry/struct.Substrate.html 
pub struct Receptor {
    substrates: Vec<Rc<Substrate>>,
    state: State,
    enzym: CatalyticCentre,
}

impl Receptor {
    pub fn detect(&self) {
        // TODO: insert check of substrate number
        let substrates: Vec<&BitBox> = self.substrates.iter()
            .map(|sub| sub.value())
            .collect();
        if self.state.detect(&substrates) {
            self.enzym.catalyse();
        }
    }
}

/// A `Protein` produces products from educt [`Substrate`]s 
/// by performing a [`Reaction`].
/// 
/// [`Substrate`]: ./struct.Substrate.html 
/// [`Reaction`]: ../chemistry/struct.Substrate.html 
pub struct CatalyticCentre {
    educts: Vec<Rc<Substrate>>,
    products: Vec<Rc<RefCell<Substrate>>>,
    reaction: Reaction,
}

impl CatalyticCentre {
    pub fn catalyse(&self) {
        // TODO: insert check of product and educt number
        let educts: Vec<&BitBox> = self.educts.iter()
            .map(|sub| sub.value())
            .collect();
        let mut product_values = self.reaction.react(&educts);
        for product in &self.products {
            // TODO: maybe switch to VecDeque and use pop_first()
            product.borrow_mut().set_value(product_values.remove(0));
        }
    }
}
