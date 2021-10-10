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
/// Tests if the function `nonlinear_normal_positve` correctly returns the supplied value.
fn test_nonlinear_normal_positve() {
    assert_eq!(nonlinear_normal_positve(0.0), 0.0);
    assert_eq!(nonlinear_normal_positve(23.45), 23.45);
    assert_eq!(nonlinear_normal_positve(-0.0), 0.0);
    assert_eq!(nonlinear_normal_positve(-23.45), 0.0);
    assert_eq!(nonlinear_normal_positve(f64::INFINITY), 0.0);
    assert_eq!(nonlinear_normal_positve(f64::NEG_INFINITY), 0.0);
    assert_eq!(nonlinear_normal_positve(f64::NAN), 0.0);
    assert_eq!(nonlinear_normal_positve(f64::MIN_POSITIVE / 2.0), 0.0);
    assert_eq!(nonlinear_normal_positve(f64::MAX * 2.0), 0.0);
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

#[test]
/// Tests the basic functionality of `ActionChain`.
fn test_action_chain() {
    let mut action_chain = ActionChain::new();
    let start_iter = action_chain.current_iteration();
    assert!(action_chain.push_action(1));
    assert!(action_chain.push_action(2));
    assert!(action_chain.push_action(3));
    assert!(!action_chain.push_action(1));
    assert_eq!(vec!(1, 2, 3), action_chain.pop_actions());
    assert_eq!(action_chain.current_iteration() - start_iter, 1);
    let empty: Vec<i32> = Vec::new();
    assert_eq!(empty, action_chain.pop_actions());
    assert_eq!(action_chain.current_iteration() - start_iter, 2);
}

#[test]
/// Tests the creation of an `ActionChain` from an `Iteration`.
fn test_action_chain_from_iteration() {
    let mut iter = Iteration::new();
    for _ in 0..100 {
        iter = iter.increment();
    }
    let action_chain: ActionChain<i32> = iter.into();
    assert_eq!(action_chain.current_iteration(), iter);
    assert_eq!(action_chain.current_iteration() - Iteration::new(), 100);
}

#[test]
/// Tests the `increment` function of the `ScalingFactor` struct.
fn test_scaling_factor_increment() {
    let mut factor = ScalingFactor::new(2.0);
    assert_ulps_eq!(factor.value(), 1.0);
    assert_eq!(factor.exponent(), 0);
    factor.increment();
    assert_ulps_eq!(factor.value(), 2.0);
    assert_eq!(factor.exponent(), 1);
    factor.increment();
    factor.increment();
    factor.decrement();
    factor.increment();
    assert_ulps_eq!(factor.value(), 8.0);
    assert_eq!(factor.exponent(), 3);
}

#[test]
/// Tests the `decrement` function of the `ScalingFactor` struct.
fn test_scaling_factor_decrement() {
    let mut factor = ScalingFactor::new(2.0);
    assert_ulps_eq!(factor.value(), 1.0);
    assert_eq!(factor.exponent(), 0);
    factor.decrement();
    assert_ulps_eq!(factor.value(), 0.5);
    assert_eq!(factor.exponent(), -1);
    factor.decrement();
    factor.decrement();
    factor.increment();
    factor.decrement();
    assert_ulps_eq!(factor.value(), 0.125);
    assert_eq!(factor.exponent(), -3);
}

#[test]
/// Tests the `ScalingFactor` struct with a different base.
fn test_scaling_factor_base() {
    let mut factor = ScalingFactor::new(1.1);
    assert_ulps_eq!(factor.value(), 1.0);
    assert_eq!(factor.exponent(), 0);
    factor.increment();
    assert_ulps_eq!(factor.value(), 1.1);
    assert_eq!(factor.exponent(), 1);
    factor.increment();
    assert_ulps_eq!(factor.value(), 1.21);
    assert_eq!(factor.exponent(), 2);
    factor.decrement();
    factor.decrement();
    factor.decrement();
    assert_ulps_eq!(factor.value(), 1.0/1.1);
    assert_eq!(factor.exponent(), -1);
}