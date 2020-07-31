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
/// Tests if the function `a_or_b` correctly returns the supplied values.
fn test_a_or_b() {
    // let mut first: u64 = 0;
    // let total: u32 = 100_000;
    // for _ in 0..total {
        let a: u64 = thread_rng().gen();
        let b: u64 = thread_rng().gen();
        let result = a_or_b(a, b);
        assert!(result == a || result == b);
    //     if result == a {
    //         first += 1;
    //     }
    // }
    // let percentage = first as f64 / total as f64;
    // assert!(percentage <= 0.75 && percentage >= 0.25);
}

#[test]
/// Tests if the function `do_a_or_b` correctly executes the supplied functions.
fn test_do_a_or_b() {
    // let mut first: u64 = 0;
    // let total: u32 = 100_000;
    // for _ in 0..total {
        let a: u64 = thread_rng().gen();
        let b: u64 = thread_rng().gen();
        let result = do_a_or_b(|| a, || b);
        assert!(result == a || result == b);
    //     if result == a {
    //         first += 1;
    //     }
    // }
    // let percentage = first as f64 / total as f64;
    // assert!(percentage <= 0.75 && percentage >= 0.25);
}

mod test_binary_substrate;
