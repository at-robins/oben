extern crate approx;
extern crate rand;

use super::*;
use rand::{Rng, thread_rng};

#[test]
#[should_panic]
/// Tests if the function `new` correctly panics on invalid values.
fn test_new_negative_total() {
    Resource::new(-1.0, 1.0);
}

#[test]
#[should_panic]
/// Tests if the function `new` correctly panics on invalid values.
fn test_new_negative_half_life() {
    Resource::new(1.0, -1.0);
}

#[test]
#[should_panic]
/// Tests if the function `new` correctly panics on invalid values.
fn test_new_nan_total() {
    Resource::new(f64::NAN, 1.0);
}

#[test]
#[should_panic]
/// Tests if the function `new` correctly panics on invalid values.
fn test_new_nan_half_life() {
    Resource::new(1.0, f64::NAN);
}

#[test]
/// Tests if the function `new` correctly creates resources.
fn test_new() {
    for _ in 1..1000 {
        // Test random valid values.
        let total: f64 = thread_rng().gen();
        let half_time: f64 = thread_rng().gen();
        let resource = Resource::new(total, half_time);
        assert!(ulps_eq!(resource.total(), total));
        assert!(ulps_eq!(resource.available(), total));
        assert!(ulps_eq!(resource.recycling(), 0.0));
        assert!(ulps_eq!(resource.half_life(), half_time));
    } {
        // Test zero.
        let resource = Resource::new(0.0, 0.0);
        assert!(ulps_eq!(resource.total(), 0.0));
        assert!(ulps_eq!(resource.available(), 0.0));
        assert!(ulps_eq!(resource.recycling(), 0.0));
        assert!(ulps_eq!(resource.half_life(), 0.0));
    }
}

#[test]
#[should_panic]
/// Tests if the function `claim_resources` correctly panics on invalid values.
fn test_claim_resources_negative() {
    Resource::default().claim_resources(-1.0);
}

#[test]
#[should_panic]
/// Tests if the function `claim_resources` correctly panics on invalid values.
fn test_claim_resources_nan() {
    Resource::default().claim_resources(f64::NAN);
}

#[test]
/// Tests if the function `claim_resources` works correctly.
fn test_claim_resources() {
    {
        // Test some resources being claimed.
        let total = 113982.4219;
        let try_to_claim = 45967.9248;
        let half_life = 29.39;
        let mut resource = Resource::new(total, half_life);
        assert!(try_to_claim < total);
        let claimed = resource.claim_resources(try_to_claim);
        assert!(ulps_eq!(try_to_claim, claimed));
        assert!(ulps_eq!(resource.total(), total - claimed));
    } {
        // Test all resources being claimed.
        let total = 113982.4219;
        let try_to_claim = 213982.936963;
        let half_life = 229.379;
        let mut resource = Resource::new(total, half_life);
        assert!(try_to_claim > total);
        let claimed = resource.claim_resources(try_to_claim);
        assert!(ulps_eq!(total, claimed));
        assert!(ulps_eq!(resource.total(), 0.0));
    } {
        // Test no resources being claimed.
        let total = 113982.4219;
        let try_to_claim = 0.0;
        let half_life = 229.379;
        let mut resource = Resource::new(total, half_life);
        let claimed = resource.claim_resources(try_to_claim);
        assert!(ulps_eq!(try_to_claim, claimed));
        assert!(ulps_eq!(resource.total(), total));
    }
}

#[test]
#[should_panic]
/// Tests if the function `claim_resources` correctly panics on invalid values.
fn test_repatriate_resources_negative() {
    Resource::default().repatriate_resources(-1.0);
}

#[test]
#[should_panic]
/// Tests if the function `claim_resources` correctly panics on invalid values.
fn test_repatriate_claim_resources_nan() {
    Resource::default().repatriate_resources(f64::NAN);
}

#[test]
/// Tests if the function `repatriate_resources` works correctly.
fn test_repatriate_resources() {
    for _ in 1..1000 {
        // Test random resources being repatriated.
        let total: f64 = thread_rng().gen();
        let repatriated: f64 = thread_rng().gen();
        let half_time: f64 = thread_rng().gen();
        let mut resource = Resource::new(total, half_time);
        assert!(ulps_eq!(resource.total(), total));
        assert!(ulps_eq!(resource.available(), total));
        assert!(ulps_eq!(resource.recycling(), 0.0));
        resource.repatriate_resources(repatriated);
        assert!(ulps_eq!(resource.total(), total + repatriated));
        assert!(ulps_eq!(resource.available(), total));
        assert!(ulps_eq!(resource.recycling(), repatriated));
    } {
        // Test no resources being repatriated.
        let total: f64 = thread_rng().gen();
        let half_time: f64 = thread_rng().gen();
        let mut resource = Resource::new(total, half_time);
        assert!(ulps_eq!(resource.total(), total));
        assert!(ulps_eq!(resource.available(), total));
        assert!(ulps_eq!(resource.recycling(), 0.0));
        resource.repatriate_resources(0.0);
        assert!(ulps_eq!(resource.total(), total));
        assert!(ulps_eq!(resource.available(), total));
        assert!(ulps_eq!(resource.recycling(), 0.0));
    }
}

#[test]
/// Tests if the function `default` does not panic.
fn test_default() {
    Resource::default();
}

#[test]
/// Tests if the function `recycle` works correctly.
fn test_recycle_resources() {
    {
        // Test half life of 1 generation.
        let total: f64 = 0.0;
        let repatriated: f64 = 4.0;
        let half_time: f64 = 1.0;
        let mut resource = Resource::new(total, half_time);
        resource.repatriate_resources(repatriated);
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 0.0));
        assert!(ulps_eq!(resource.recycling(), 4.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 2.0));
        assert!(ulps_eq!(resource.recycling(), 2.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.0));
        assert!(ulps_eq!(resource.recycling(), 1.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.5));
        assert!(ulps_eq!(resource.recycling(), 0.5));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.75));
        assert!(ulps_eq!(resource.recycling(), 0.25));
        
    } {
        // Test half life of 3 generation.
        let total: f64 = 0.0;
        let repatriated: f64 = 4.0;
        let half_time: f64 = 3.0;
        let mut resource = Resource::new(total, half_time);
        resource.repatriate_resources(repatriated);
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 0.0));
        assert!(ulps_eq!(resource.recycling(), 4.0));
        resource.recycle();
        resource.recycle();
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 2.0));
        assert!(ulps_eq!(resource.recycling(), 2.0));
        resource.recycle();
        resource.recycle();
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.0));
        assert!(ulps_eq!(resource.recycling(), 1.0));
        resource.recycle();
        resource.recycle();
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.5));
        assert!(ulps_eq!(resource.recycling(), 0.5));
        resource.recycle();
        resource.recycle();
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.75));
        assert!(ulps_eq!(resource.recycling(), 0.25));
        
    } {
        // Test half life of 0.5 generation.
        let total: f64 = 0.0;
        let repatriated: f64 = 4.0;
        let half_time: f64 = 0.5;
        let mut resource = Resource::new(total, half_time);
        resource.repatriate_resources(repatriated);
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 0.0));
        assert!(ulps_eq!(resource.recycling(), 4.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.0));
        assert!(ulps_eq!(resource.recycling(), 1.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 3.75));
        assert!(ulps_eq!(resource.recycling(), 0.25));
        
    } {
        // Test half life of 0.0 generation.
        let total: f64 = 0.0;
        let repatriated: f64 = 4.0;
        let half_time: f64 = 0.0;
        let mut resource = Resource::new(total, half_time);
        resource.repatriate_resources(repatriated);
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 0.0));
        assert!(ulps_eq!(resource.recycling(), 4.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 4.0));
        assert!(ulps_eq!(resource.recycling(), 0.0));
        resource.recycle();
        assert!(ulps_eq!(resource.total(), 4.0));
        assert!(ulps_eq!(resource.available(), 4.0));
        assert!(ulps_eq!(resource.recycling(), 0.0));
        
    }
}