//! The `binary_chemistry` module contains elementary bit actions
//! for a binary evolutionary network.
extern crate bitvec;
extern crate rand;
extern crate serde;

use bitvec::vec::BitVec;
use rand::{distributions::{Distribution, Standard}, Rng, thread_rng};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use super::BinarySubstrate;
use super::super::chemistry::{Reaction, State};
use super::super::gene::CrossOver;
use super::super::helper::do_a_or_b;

/// A `BinaryState` is an elementary operation for comparing binary substrates.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub enum BinaryState {
    /// A state operation comparing two substrates for equality.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[255, 12, 34]);
    /// let b: BinarySubstrate = a.clone();
    /// assert!(BinaryState::Equals.detect(&[&a,&b]));
    /// ```
    Equals,
    /// A state operation comparing two substrates for inequality.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[255, 12, 34]);
    /// let b: BinarySubstrate = a.clone();
    /// assert!(!BinaryState::Not.detect(&[&a,&b]));
    /// ```
    Not,
    /// A state operation checking if one substrate is greater than the other.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[128, 12, 34]);
    /// let b: BinarySubstrate = BitBox::from_slice(&[127]);
    /// assert!(BinaryState::Greater.detect(&[&a, &b]));
    /// ```
    Greater,
    /// A state operation checking if one substrate is samller than the other.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[128, 12, 34]);
    /// let b: BinarySubstrate = BitBox::from_slice(&[127]);
    /// assert!(BinaryState::Lesser.detect(&[&b, &a]));
    /// ```
    Lesser,
    /// A state operation checking if one substrate is greater or equal to the other.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[255, 12, 34]);
    /// let b: BinarySubstrate = a.clone();
    /// assert!(BinaryState::GreaterOrEqual.detect(&[&a,&b]));
    /// ```
    GreaterOrEqual,
    /// A state operation checking if one substrate is smaller or equal to the other.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[255, 12, 34]);
    /// let b: BinarySubstrate = a.clone();
    /// assert!(BinaryState::LesserOrEqual.detect(&[&a,&b]));
    /// ```
    LesserOrEqual,
    /// A state operation checking if a substrate is completely set.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[u8::max_value(), u8::max_value()]);
    /// assert!(BinaryState::All.detect(&[&a]));
    /// ```
    All,
    /// A state operation checking if a substrate is not completely unset.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[0, 1, 0]);
    /// assert!(BinaryState::Some.detect(&[&a]));
    /// ```
    Some,
    /// A state operation checking if a substrate is not completely set.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[0, 1, 0]);
    /// assert!(BinaryState::NotAll.detect(&[&a]));
    /// ```
    NotAll,
    /// A state operation checking if a substrate is completely unset.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: BinarySubstrate = BitBox::from_slice(&[0, 0, 0]);
    /// assert!(BinaryState::None.detect(&[&a]));
    /// ```
    None,
    /// A state operation always returning `true`.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: Vec<&BinarySubstrate> = Vec::new();
    /// assert!(BinaryState::Always.detect(&a));
    /// ```
    Always,
    /// A state operation always returning `false`.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::State;
    ///
    /// let a: Vec<&BinarySubstrate> = Vec::new();
    /// assert!(!BinaryState::Never.detect(&a));
    /// ```
    Never,
}

impl State<BinarySubstrate> for BinaryState {
    /// Compares a number of substrates for a logical property.
    ///
    /// # Parameters
    ///
    /// * `substrates` - the substrates to perform the logical operation on
    ///
    /// # Panics
    ///
    /// If the number of substrates is not exactly equal to the required one.
    fn detect(&self, substrates: &[&BinarySubstrate]) -> bool {
        assert_eq!(substrates.len(), self.get_substrate_number(),
            "The number of required substrates is {}, but {} substrates were supplied.",
            self.get_substrate_number(), substrates.len());

        match self {
            BinaryState::Equals => substrates[0] == substrates[1],
            BinaryState::Not => substrates[0] != substrates[1],
            BinaryState::Greater => substrates[0] > substrates[1],
            BinaryState::Lesser => substrates[0] < substrates[1],
            BinaryState::GreaterOrEqual => substrates[0] >= substrates[1],
            BinaryState::LesserOrEqual => substrates[0] <= substrates[1],
            BinaryState::All => substrates[0].all(),
            BinaryState::Some => substrates[0].some(),
            BinaryState::None => substrates[0].not_any(),
            BinaryState::NotAll => substrates[0].not_all(),
            BinaryState::Always => true,
            BinaryState::Never => false,
        }
    }

    /// Returns the number of substrates required to perform the logical comparison.
    fn get_substrate_number(&self) -> usize {
        match self {
            BinaryState::Equals => 2,
            BinaryState::Not => 2,
            BinaryState::Greater => 2,
            BinaryState::Lesser => 2,
            BinaryState::GreaterOrEqual => 2,
            BinaryState::LesserOrEqual => 2,
            BinaryState::All => 1,
            BinaryState::Some => 1,
            BinaryState::NotAll => 1,
            BinaryState::None => 1,
            BinaryState::Always => 0,
            BinaryState::Never => 0,
        }
    }

    fn random() -> Self {
        thread_rng().gen()
    }
}

impl CrossOver for BinaryState {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        do_a_or_b(|| self.clone(), || other.clone())
    }
}

impl Distribution<BinaryState> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BinaryState {
        match rng.gen_range(0u8, 12) {
            0 => BinaryState::Equals,
            1 => BinaryState::Not,
            2 => BinaryState::Greater,
            3 => BinaryState::Lesser,
            4 => BinaryState::GreaterOrEqual,
            5 => BinaryState::LesserOrEqual,
            6 => BinaryState::All,
            7 => BinaryState::Some,
            8 => BinaryState::NotAll,
            9 => BinaryState::None,
            10 => BinaryState::Always,
            11 => BinaryState::Never,
            _ => panic!("A random number with no matching state was created.")
        }
    }
}

/// A `BinaryReaction` represents an elementary operation for modification of binary substrates.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub enum BinaryReaction {
    /// A binary `AND` operation.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt1_value: BinarySubstrate = BitBox::from_slice(&[255]);
    /// let educt2_value: BinarySubstrate = BitBox::from_slice(&[32]);
    /// let educts = vec!(&educt1_value, &educt2_value);
    ///
    /// let product_value: BinarySubstrate = BitBox::from_element(32);
    /// let product = vec!(product_value);
    ///
    /// assert_eq!(BinaryReaction::And.react(&educts[..]), product);
    /// ```
    And,
    /// A binary `OR` operation.
    Or,
    /// A binary `XOR` operation.
    XOr,
    /// A binary right shift operation.
    ShiftRight,
    /// A binary left shift operation.
    ShiftLeft,
    /// An operation, appending a binary value to another one.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt1_value: BinarySubstrate = BitBox::from_slice(&[255, 32]);
    /// let educt2_value: BinarySubstrate = BitBox::from_slice(&[4, 35]);
    /// let educts = vec!(&educt1_value, &educt2_value);
    ///
    /// let product_value: BinarySubstrate = BitBox::from_slice(&[255, 32, 4, 35]);
    /// let product = vec!(product_value);
    ///
    /// assert_eq!(BinaryReaction::Append.react(&educts[..]), product);
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
    /// use oben::evolution::binary::BinarySubstrate;;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt_values: BinarySubstrate = BitBox::from_slice(&[255, 32]);
    /// let educts = vec!(&educt_values);
    ///
    /// assert_eq!(BinaryReaction::Duplicate.react(&educts[..])[0], *educts[0]);
    /// ```
    Duplicate,
    /// No function.
    ///
    /// # Example
    /// ```
    /// use oben::evolution::binary::BinarySubstrate;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educts: Vec<&BinarySubstrate> = vec!();
    ///
    /// assert!(BinaryReaction::Misfunction.react(&educts[..]).is_empty());
    /// ```
    Misfunction,
    /// A random binary number.
    ///
    /// # Example
    /// ```
    /// extern crate bitvec;
    ///
    /// use bitvec::boxed::BitBox;
    /// use oben::evolution::binary::BinarySubstrate;;
    /// use oben::evolution::chemistry::Reaction;
    ///
    /// let educt_values: BinarySubstrate = BitBox::from_slice(&[255, 32]);
    /// let educts = vec!(&educt_values);
    ///
    /// assert_eq!(BinaryReaction::Random.react(&educts[..])[0].len(), educts[0].len());
    /// ```
    Random,
}

impl Reaction<BinarySubstrate> for BinaryReaction {
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
    fn react(&self, educts: &[&BinarySubstrate]) -> Vec<BinarySubstrate> {
        assert_eq!(educts.len(), self.get_educt_number(),
            "The number of required educts is {}, but {} educts were supplied.",
            self.get_educt_number(), educts.len());

        let products = match self {
            BinaryReaction::And => vec!(educts[0].clone() & educts[1].clone()),
            BinaryReaction::Or => vec!(educts[0].clone() | educts[1].clone()),
            BinaryReaction::XOr => vec!(educts[0].clone() ^ educts[1].clone()),
            // TODO: Re-implement the bit shift operations as soon as native BitBox shift gets fixed.
            // Currently casting to a BitVec is needed otherwise BitBoxes of length n * container
            // type (here u8 => n * 8) will panic upon shifting.
            BinaryReaction::ShiftRight => vec!((educts[0].to_vec() >> educts[1].len()).into()),
            BinaryReaction::ShiftLeft => vec!((educts[0].to_vec() << educts[1].len()).into()),
            BinaryReaction::Append => {
                let mut a = BitVec::from(educts[0].clone());
                let mut b = BitVec::from(educts[1].clone());
                a.append(&mut b);
                vec!(a.into())
            },
            BinaryReaction::Inverse => vec!(!educts[0].clone()),
            BinaryReaction::Reverse => {
                let mut a = educts[0].clone();
                a.reverse();
                vec!(a)
            },
            BinaryReaction::Duplicate => vec!(educts[0].clone()),
            BinaryReaction::Misfunction => vec!(),
            BinaryReaction::Random => {
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
    fn get_educt_number(&self) -> usize {
        match self {
            BinaryReaction::And => 2,
            BinaryReaction::Or => 2,
            BinaryReaction::XOr => 2,
            BinaryReaction::ShiftRight => 2,
            BinaryReaction::ShiftLeft => 2,
            BinaryReaction::Append => 2,
            BinaryReaction::Inverse => 1,
            BinaryReaction::Reverse => 1,
            BinaryReaction::Duplicate => 1,
            BinaryReaction::Misfunction => 0,
            BinaryReaction::Random => 1,
        }
    }

    /// Returns the number of products required to perform the reaction.
    fn get_product_number(&self) -> usize {
        match self {
            BinaryReaction::And => 1,
            BinaryReaction::Or => 1,
            BinaryReaction::XOr => 1,
            BinaryReaction::ShiftRight => 1,
            BinaryReaction::ShiftLeft => 1,
            BinaryReaction::Append => 1,
            BinaryReaction::Inverse => 1,
            BinaryReaction::Reverse => 1,
            BinaryReaction::Duplicate => 1,
            BinaryReaction::Misfunction => 0,
            BinaryReaction::Random => 1,
        }
    }

    fn random() -> Self {
        thread_rng().gen()
    }
}

impl CrossOver for BinaryReaction {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        do_a_or_b(|| self.clone(), || other.clone())
    }
}

impl Distribution<BinaryReaction> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BinaryReaction {
        match rng.gen_range(0u8, 11) {
            0 => BinaryReaction::And,
            1 => BinaryReaction::Or,
            2 => BinaryReaction::XOr,
            3 => BinaryReaction::ShiftRight,
            4 => BinaryReaction::ShiftLeft,
            5 => BinaryReaction::Append,
            6 => BinaryReaction::Inverse,
            7 => BinaryReaction::Reverse,
            8 => BinaryReaction::Duplicate,
            9 => BinaryReaction::Misfunction,
            10 => BinaryReaction::Random,
            _ => panic!("A random number with no matching reaction was created.")
        }
    }
}
