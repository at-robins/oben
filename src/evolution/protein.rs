extern crate bitvec;

use std::rc::Rc;
use bitvec::boxed::BitBox;

/// A `Substrate` represents a chemical entity of a specific value. Upon change
/// [`Receptor`]s of the substrate will be notified.
/// 
/// [`Receptor`]: ./struct.Receptor.html
pub struct Substrate {
    value: BitBox,
    receptors: Vec<Rc<Receptor>>,
}

/// A `Receptor` is a sensor for [`Substrate`] changes and a trigger
/// for [`Reaction`]s.
/// 
/// [`Substrate`]: ./struct.Substrate.html 
/// [`Reaction`]: ../chemistry/struct.Substrate.html 
pub struct Receptor {
    
}

pub struct Protein {
    
}