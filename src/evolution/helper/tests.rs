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
/// Tests if the function 'increment' of `Iteration` works as expected.
fn test_iteration_increment() {
    let start = Iteration::new();
    let mut last;
    let mut current = Iteration::new();
    for i in 0..1000 {
        assert_eq!(current - start, i);
        last = current;
        current = current.increment();
        assert_eq!(current - last, 1);
    }
}

#[test]
/// Tests if the function 'difference' of `Iteration` works as expected.
fn test_iteration_difference() {
    // Test difference in range.
    {
        assert_eq!(Iteration::new().difference(Iteration::new()), 0);
        let test_difference = 240924246i32;
        let a = Iteration::new();
        let b = Iteration{current_iteration: test_difference.into()};
        assert_eq!(b.difference(a), test_difference);
        assert_eq!(a.difference(b), -test_difference);
    }
    // Test difference out of range.
    {
        let a = Iteration::new();
        let b = Iteration{current_iteration: i32::MAX as i128 + 429828};
        assert_eq!(b.difference(a), i32::MAX);
        assert_eq!(a.difference(b), i32::MIN);
    }
}

#[test]
/// Tests if the operation 'sub' of `Iteration` works as expected.
fn test_iteration_sub() {
    // Test difference in range.
    {
        assert_eq!(Iteration::new() - Iteration::new(), 0);
        let test_difference = 240924246i32;
        let a = Iteration::new();
        let b = Iteration{current_iteration: test_difference.into()};
        assert_eq!(b - a, test_difference);
        assert_eq!(a - b, -test_difference);
    }
    // Test difference out of range.
    {
        let a = Iteration::new();
        let b = Iteration{current_iteration: i32::MAX as i128 + 429828};
        assert_eq!(b - a, i32::MAX);
        assert_eq!(a - b, i32::MIN);
    }
}

#[test]
/// Tests if equality of `Iteration` works as expected.
fn test_iteration_equality() {
    let mut a = Iteration::new();
    let mut b = Iteration::new();
    for _ in 0..1000 {
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
    for _ in 0..1000 {
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
