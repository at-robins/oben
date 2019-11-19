//! The `chemistry` module contains elementary bit actions 
//! for the evolutionary network.
extern crate bitvec;

use bitvec::{boxed::BitBox, vec::BitVec};

/// A `State` is an elementary operation for comparing two substrates.
pub enum State {
    /// A state operation comparing two substrates for equality.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 12, 34]);
    /// let b: BitBox = a.clone();
    /// assert!(State::Equals.detect(&[&a,&b]));
    /// ```
    Equals,
    /// A state operation comparing two substrates for inequality.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 12, 34]);
    /// let b: BitBox = a.clone();
    /// assert!(!State::Not.detect(&[&a,&b]));
    /// ```
    Not,
    /// A state operation checking if one substrate is greater than the other.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 12, 34]);
    /// let b: BitBox = BitBox::from_slice(&[254u8]);
    /// assert!(State::Greater.detect(&[&a,&b]));
    /// ```
    Greater,
    /// A state operation checking if one substrate is samller than the other.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 12, 34]);
    /// let b: BitBox = BitBox::from_slice(&[254u8]);
    /// assert!(State::Lesser.detect(&[&b,&a]));
    /// ```
    Lesser,
    /// A state operation checking if one substrate is greater or equal to the other.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 12, 34]);
    /// let b: BitBox = a.clone();
    /// assert!(State::GreaterOrEqual.detect(&[&a,&b]));
    /// ```
    GreaterOrEqual,
    /// A state operation checking if one substrate is smaller or equal to the other.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 12, 34]);
    /// let b: BitBox = a.clone();
    /// assert!(State::LesserOrEqual.detect(&[&a,&b]));
    /// ```
    LesserOrEqual,
    /// A state operation checking if a substrate is completely set.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[255u8, 255, 255]);
    /// assert!(State::All.detect(&[&a]));
    /// ```
    All,
    /// A state operation checking if a substrate is completely unset.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[0u8, 0, 0]);
    /// assert!(State::None.detect(&[&a]));
    /// ```
    Some,
    /// A state operation checking if a substrate is partially set.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[0u8, 1, 0]);
    /// assert!(State::NotAll.detect(&[&a]));
    /// ```
    NotAll,
    /// A state operation checking if a substrate is completely unset.
    /// 
    /// # Example
    /// ```
    /// extern crate bitvec;
    /// 
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    /// 
    /// let a: BitBox = BitBox::from_slice(&[0u8, 0, 0]);
    /// assert!(State::None.detect(&[&a]));
    /// ```
    None,
}

impl State {
    /// Compares a number of substrates for a logical property.
    /// 
    /// # Parameters
    /// 
    /// * `substrates` - the substrates to perform the logical operation on
    /// 
    /// # Panics
    /// 
    /// If the number of substrates is not exactly equal to the required one. 
    pub fn detect(&self, substrates: &[&BitBox]) -> bool {
        assert_eq!(substrates.len(), self.get_substrate_number(), 
            "The number of required substrates is {}, but {} substrates were supplied", 
            self.get_substrate_number(), substrates.len());

        match self {
            State::Equals => substrates[0] == substrates[1],
            State::Not => substrates[0] != substrates[1],
            State::Greater => substrates[0] > substrates[1],
            State::Lesser => substrates[0] < substrates[1],
            State::GreaterOrEqual => substrates[0] >= substrates[1],
            State::LesserOrEqual => substrates[0] <= substrates[1],
            State::All => substrates[0].all(),
            State::Some => substrates[0].some(),
            State::None => substrates[0].not_any(),
            State::NotAll => substrates[0].not_all(),
        }
    }
    
    /// Returns the number of substrates required to perform the logical comparison.
    pub fn get_substrate_number(&self) -> usize {
        match self {
            State::Equals => 2,
            State::Not => 2,
            State::Greater => 2,
            State::Lesser => 2,
            State::GreaterOrEqual => 2,
            State::LesserOrEqual => 2,
            State::All => 1,
            State::Some => 1,
            State::NotAll => 1,
            State::None => 1,
        }
    }
    
}

pub enum Reaction {
    
}
