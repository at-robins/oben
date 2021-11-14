extern crate rand;

use rand::{Rng, thread_rng};
use super::*;


#[test]
/// Tests if the function `as_64` can correctly convert a 64 bit `BitBox` into a 64 bit array.
fn test_as_64() {
    let test_64: [u8; 8] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_64);
    assert_eq!(test_64, as_64(&test_substrate));
}

#[test]
#[should_panic]
/// Tests if the function `as_64` panics if the passed array contains less than 64 bit.
fn test_panic_less_as_64() {
    let test_64: [u8; 7] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_64);
    as_64(&test_substrate);
}

#[test]
#[should_panic]
/// Tests if the function `as_64` panics if the passed array contains more than 64 bit.
fn test_panic_greater_as_64() {
    let test_64: [u8; 9] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_64);
    as_64(&test_substrate);
}

#[test]
/// Tests if the function `as_f64` can correctly convert a 64 bit `BitBox` into a `f64`.
fn test_as_f64() {
    let test_f64: f64 = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_f64.to_be_bytes());
    assert_eq!(test_f64, as_f64(&test_substrate));
}

#[test]
#[should_panic]
/// Tests if the function `as_f64` panics if the passed substrate contains less than 64 bit.
fn test_panic_less_as_f64() {
    let test_f64: [u8; 7] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_f64);
    as_f64(&test_substrate);
}

#[test]
#[should_panic]
/// Tests if the function `as_f64` panics if the passed substrate contains more than 64 bit.
fn test_panic_greater_as_f64() {
    let test_f64: [u8; 9] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_f64);
    as_f64(&test_substrate);
}

#[test]
/// Tests if the function `f64_as_binary` can correctly convert a `f64` into a 64 bit `BitBox`.
fn test_f64_as_binary() {
    let test_f64: f64 = thread_rng().gen();
    let test_substrate: BinarySubstrate = f64_to_binary(test_f64);
    assert_eq!(test_f64, as_f64(&test_substrate));
}

#[test]
/// Tests if the function `as_u64` can correctly convert a 64 bit `BitBox` into a `u64`.
fn test_as_u64() {
    let test_u64: u64 = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_u64.to_be_bytes());
    assert_eq!(test_u64, as_u64(&test_substrate));
}

#[test]
#[should_panic]
/// Tests if the function `as_u64` panics if the passed substrate contains less than 64 bit.
fn test_panic_less_as_u64() {
    let test_u64: [u8; 7] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_u64);
    as_u64(&test_substrate);
}

#[test]
#[should_panic]
/// Tests if the function `as_u64` panics if the passed substrate contains more than 64 bit.
fn test_panic_greater_as_u64() {
    let test_u64: [u8; 9] = thread_rng().gen();
    let test_substrate: BinarySubstrate = BitBox::from_slice(&test_u64);
    as_u64(&test_substrate);
}

#[test]
/// Tests if the function `u64_as_binary` can correctly convert a `u64` into a 64 bit `BitBox`.
fn test_u64_as_binary() {
    let test_u64: u64 = thread_rng().gen();
    let test_substrate: BinarySubstrate = u64_to_binary(test_u64);
    assert_eq!(test_u64, as_u64(&test_substrate));
}

mod test_binary_substrate;
