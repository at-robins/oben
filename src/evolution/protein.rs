extern crate bitvec;

use std::rc::Rc;
use bitvec::boxed::BitBox;
use super::chemistry::Reaction;

/// A `Substrate` represents a chemical entity of a specific value. Upon change
/// [`Receptor`]s of the substrate will be notified.
/// 
/// [`Receptor`]: ./struct.Receptor.html
pub struct Substrate {
    value: BitBox,
    receptors: Vec<Rc<Receptor>>,
}

impl Substrate {
    /// Set the value of this substrate.alloc
    /// 
    /// # Parameters
    /// 
    /// * `value` - the new value of the substrate
    pub fn set_value(&mut self, value: BitBox) {
        self.value = value;
    }
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
/// 
/// [`Substrate`]: ./struct.Substrate.html 
/// [`Reaction`]: ../chemistry/struct.Substrate.html 
pub struct Receptor {
    substrates: Vec<Rc<Substrate>>,
    enzym: Protein,
}

/// A `Protein` produces products from educt [`Substrate`]s 
/// by performing a [`Reaction`].
/// 
/// [`Substrate`]: ./struct.Substrate.html 
/// [`Reaction`]: ../chemistry/struct.Substrate.html 
pub struct Protein {
    educts: Vec<Rc<Substrate>>,
    products: Vec<Rc<Substrate>>,
    reaction: Reaction,
}

/*impl Protein {
    pub fn catalyse(&self) {
        self.reaction.react();
    }
}*/
