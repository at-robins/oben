//! The `chemistry` module contains elementary bit actions
//! for the evolutionary network.
extern crate bitvec;
extern crate rand;
extern crate serde;

use bitvec::{prelude::BitStore, boxed::BitBox, order::BitOrder, vec::BitVec};
use rand::{distributions::{Distribution, Standard}, Rng};
use serde::{Deserialize, Serialize};

/// A `State` is an elementary operation for comparing binary substrates.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
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
    /// let a: BitBox = BitBox::from_slice(&[255usize, 12, 34]);
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
    /// let a: BitBox = BitBox::from_slice(&[255usize, 12, 34]);
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
    /// let a: BitBox = BitBox::from_slice(&[255usize, 12, 34]);
    /// let b: BitBox = BitBox::from_slice(&[254usize]);
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
    /// let a: BitBox = BitBox::from_slice(&[255usize, 12, 34]);
    /// let b: BitBox = BitBox::from_slice(&[254usize]);
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
    /// let a: BitBox = BitBox::from_slice(&[255usize, 12, 34]);
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
    /// let a: BitBox = BitBox::from_slice(&[255usize, 12, 34]);
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
    /// let a: BitBox = BitBox::from_slice(&[usize::max_value(), usize::max_value()]);
    /// assert!(State::All.detect(&[&a]));
    /// ```
    All,
    /// A state operation checking if a substrate is not completely unset.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BitBox = BitBox::from_slice(&[0usize, 1, 0]);
    /// assert!(State::Some.detect(&[&a]));
    /// ```
    Some,
    /// A state operation checking if a substrate is not completely set.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BitBox = BitBox::from_slice(&[0usize, 1, 0]);
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
    /// let a: BitBox = BitBox::from_slice(&[0usize, 0, 0]);
    /// assert!(State::None.detect(&[&a]));
    /// ```
    None,
    /// A state operation always returning `true`.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    ///
    /// assert!(State::Always.detect(&[]));
    /// ```
    Always,
    /// A state operation always returning `false`.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::State;
    ///
    /// assert!(!State::Never.detect(&[]));
    /// ```
    Never,
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
    pub fn detect<O,S>(&self, substrates: &[&BitBox<O,S>]) -> bool
        where O: BitOrder, S: BitStore {
        assert_eq!(substrates.len(), self.get_substrate_number(),
            "The number of required substrates is {}, but {} substrates were supplied.",
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
            State::Always => true,
            State::Never => false,
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
            State::Always => 0,
            State::Never => 0,
        }
    }

}

impl Distribution<State> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> State {
        match rng.gen_range(0u8, 12) {
            0 => State::Equals,
            1 => State::Not,
            2 => State::Greater,
            3 => State::Lesser,
            4 => State::GreaterOrEqual,
            5 => State::LesserOrEqual,
            6 => State::All,
            7 => State::Some,
            8 => State::NotAll,
            9 => State::None,
            10 => State::Always,
            11 => State::Never,
            _ => panic!("A random number with no matching state was created.")
        }
    }
}

/// A `Reaction` represents an elementary operation for modification of binary substrates.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub enum Reaction {
    /// A binary `AND` operation.
    And,
    /// A binary `OR` operation.
    Or,
    /// A binary `XOR` operation.
    XOr,
    /// A binary right shift operation.
    ShiftRight,
    /// A binary left shift operation.
    ShiftLeft,
    /// A binary addition.
    Add,
    /// An operation, appending a binary value to another one.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt1_value = BitBox::from_slice(&[255usize, 32]);
    /// let educt2_value = BitBox::from_slice(&[4usize, 35]);
    /// let educts = vec!(&educt1_value, &educt2_value);
    ///
    /// let product_value: BitBox = BitBox::from_slice(&[255usize, 32, 4, 35]);
    /// let product = vec!(product_value);
    ///
    /// assert_eq!(Reaction::Append.react(&educts[..]), product);
    /// ```
    Append,
    /// A binary inversion.
    Inverse,
    /// A binary reverse.
    Reverse,
    /// A duplication of a binary value.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt_values = BitBox::from_slice(&[255usize, 32]);
    /// let educts = vec!(&educt_values);
    ///
    /// assert_eq!(Reaction::Duplicate.react(&educts[..])[0], *educts[0]);
    /// ```
    Duplicate,
    /// No function.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educts: Vec<&BitBox> = vec!();
    ///
    /// assert!(Reaction::Misfunction.react(&educts[..]).is_empty());
    /// ```
    Misfunction,
    /// A random binary number.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt_values = BitBox::from_slice(&[255usize, 32]);
    /// let educts = vec!(&educt_values);
    ///
    /// assert_eq!(Reaction::Random.react(&educts[..])[0].len(), educts[0].len());
    /// ```
    Random,
}

impl Reaction {
    /// Performs the specified reaction and returns the products.
    ///
    /// # Parameters
    ///
    /// * `educts` - the educts to convert into products
    ///
    /// # Panics
    ///
    /// If the number of supplied educts and created products is not
    /// exactly equal to the required one.
    pub fn react<O, S>(&self, educts: &[&BitBox<O, S>]) -> Vec<BitBox<O, S>>
        where O: BitOrder, S: BitStore {
        assert_eq!(educts.len(), self.get_educt_number(),
            "The number of required educts is {}, but {} educts were supplied.",
            self.get_educt_number(), educts.len());

        let products = match self {
            Reaction::And => vec!(educts[0].clone() & educts[1].clone()),
            Reaction::Or => vec!(educts[0].clone() | educts[1].clone()),
            Reaction::XOr => vec!(educts[0].clone() ^ educts[1].clone()),
            Reaction::ShiftRight => vec!(educts[0].clone() >> educts[1].len()),
            Reaction::ShiftLeft => vec!(educts[0].clone() << educts[1].len()),
            Reaction::Add => vec!(educts[0].clone() + educts[1].clone()),
            Reaction::Append => {
                let mut a = BitVec::from(educts[0].clone());
                let mut b = BitVec::from(educts[1].clone());
                a.append(&mut b);
                vec!(a.into())
            },
            Reaction::Inverse => vec!(!educts[0].clone()),
            Reaction::Reverse => {
                let mut a = educts[0].clone();
                a.reverse();
                vec!(a)
            },
            Reaction::Duplicate => vec!(educts[0].clone()),
            Reaction::Misfunction => vec!(),
            Reaction::Random => {
                let mut random_bits = BitVec::with_capacity(educts[0].len());
                for _i in 0..educts[0].len() {
                    random_bits.push(rand::random())
                }
                vec!(random_bits.into())
            },
        };

        assert_eq!(products.len(), self.get_product_number(),
            "The number of required products is {}, but {} products were created.",
            self.get_product_number(), products.len());

        products
    }

    /// Returns the number of educts required to perform the reaction.
    pub fn get_educt_number(&self) -> usize {
        match self {
            Reaction::And => 2,
            Reaction::Or => 2,
            Reaction::XOr => 2,
            Reaction::ShiftRight => 2,
            Reaction::ShiftLeft => 2,
            Reaction::Add => 2,
            Reaction::Append => 2,
            Reaction::Inverse => 1,
            Reaction::Reverse => 1,
            Reaction::Duplicate => 1,
            Reaction::Misfunction => 0,
            Reaction::Random => 1,
        }
    }

    /// Returns the number of products required to perform the reaction.
    pub fn get_product_number(&self) -> usize {
        match self {
            Reaction::And => 1,
            Reaction::Or => 1,
            Reaction::XOr => 1,
            Reaction::ShiftRight => 1,
            Reaction::ShiftLeft => 1,
            Reaction::Add => 1,
            Reaction::Append => 1,
            Reaction::Inverse => 1,
            Reaction::Reverse => 1,
            Reaction::Duplicate => 1,
            Reaction::Misfunction => 0,
            Reaction::Random => 1,
        }
    }
}

impl Distribution<Reaction> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Reaction {
        match rng.gen_range(0u8, 12) {
            0 => Reaction::And,
            1 => Reaction::Or,
            2 => Reaction::XOr,
            3 => Reaction::ShiftRight,
            4 => Reaction::ShiftLeft,
            5 => Reaction::Add,
            6 => Reaction::Append,
            7 => Reaction::Inverse,
            8 => Reaction::Reverse,
            9 => Reaction::Duplicate,
            10 => Reaction::Misfunction,
            11 => Reaction::Random,
            _ => panic!("A random number with no matching reaction was created.")
        }
    }
}
