use super::*;

#[test]
/// Tests if the function `cross_over` of the `CrossOver` trait can correctly recombine
/// substrates.
fn test_vec_cross_over() {
    // Test vectors of the same length.
    {
        let a: [u8; 16] = thread_rng().gen();
        let b: [u8; 16] = thread_rng().gen();
        let vec_a: Vec<BinarySubstrate> = vec!(BitBox::from_slice(&a));
        let vec_b: Vec<BinarySubstrate> = vec!(BitBox::from_slice(&b));
        let recombined = vec_a.cross_over(&vec_b);
        for i in 0..vec_a[0].len() {
            assert!(recombined[0].get(i) == vec_a[0].get(i)
                || recombined[0].get(i) == vec_b[0].get(i));
        }
    }
    // Test vectors of different length.
    {
        let a: [u8; 16] = thread_rng().gen();
        let b: [u8; 16] = thread_rng().gen();
        let c: [u8; 16] = thread_rng().gen();
        let vec_a: Vec<BinarySubstrate> = vec!(BitBox::from_slice(&a));
        let vec_b: Vec<BinarySubstrate> = vec!(BitBox::from_slice(&b), BitBox::from_slice(&c));
        let recombined = vec_a.cross_over(&vec_b);
        for i in 0.. vec_a[0].len() {
            assert!(recombined[0].get(i) == vec_a[0].get(i)
                || recombined[0].get(i) == vec_b[0].get(i));
        }
        assert!(recombined.len() > 0 && recombined.len() < 3);
        if recombined.len() == 2 {
            assert!(recombined[1] == vec_b[1]);
        }
    }
}
