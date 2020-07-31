use super::*;

#[test]
/// Tests if the function `is_similar` of the `CrossOver` trait can correctly detect similar
/// substrates.
fn test_cross_over_is_similar() {
    // Test similar substrates.
    let a: [u8; 8] = thread_rng().gen();
    let test_substrate_a: BinarySubstrate = BitBox::from_slice(&a);
    let test_substrate_b: BinarySubstrate = test_substrate_a.clone();
    assert!(test_substrate_a.is_similar(&test_substrate_b));
    assert!(test_substrate_b.is_similar(&test_substrate_a));
    // Test non-similar substrates.
    let b: [u8; 16] = thread_rng().gen();
    let test_substrate_b: BinarySubstrate = BitBox::from_slice(&b);
    assert!(!test_substrate_a.is_similar(&test_substrate_b));
    assert!(!test_substrate_b.is_similar(&test_substrate_a));
}

#[test]
/// Tests if the function `cross_over` of the `CrossOver` trait can correctly recombine
/// substrates.
fn test_cross_over_cross_over() {
    // Test similar substrates.
    let a: [u8; 16] = thread_rng().gen();
    let b: [u8; 16] = thread_rng().gen();
    let test_substrate_a: BinarySubstrate = BitBox::from_slice(&a);
    let test_substrate_b: BinarySubstrate = BitBox::from_slice(&b);
    let recombined = test_substrate_a.cross_over(&test_substrate_b);
    for i in 0.. test_substrate_a.len() {
        assert!(recombined.get(i) == test_substrate_a.get(i)
            || recombined.get(i) == test_substrate_b.get(i));
    }
    // Test non-similar substrates.
    let a: [u8; 8] = thread_rng().gen();
    let b: [u8; 16] = thread_rng().gen();
    let test_substrate_a: BinarySubstrate = BitBox::from_slice(&a);
    let test_substrate_b: BinarySubstrate = BitBox::from_slice(&b);
    let recombined = test_substrate_a.cross_over(&test_substrate_b);
    assert!(recombined == test_substrate_a
        || recombined == test_substrate_b);
}
