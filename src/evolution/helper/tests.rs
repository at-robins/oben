use super::*;

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

#[test]
/// Tests if equality of `Iteration` works as expected.
fn test_iteration_equality() {
    let mut a = Iteration::new();
    let mut b = Iteration::new();
    for _ in 0..100 {
        assert_eq!(a, b);
        a = a.increment();
        assert_ne!(a, b);
        b = b.increment();
    }
}

#[test]
/// Tests if ordering of `Iteration` works as expected.
fn test_iteration_ordering() {
    let mut a = Iteration::new();
    let mut b = Iteration::new();
    for _ in 0..100 {
        assert!(a == b);
        assert!(!(a > b));
        assert!(!(a < b));
        assert!(!(b > a));
        assert!(!(b < a));
        a = a.increment();
        assert!(a != b);
        assert!(a > b);
        assert!(!(a < b));
        assert!(!(b > a));
        assert!(b < a);
        b = b.increment();
    }
}
