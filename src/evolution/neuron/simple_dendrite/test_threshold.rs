use super::*;

#[test]
/// Tests if the function `new` correctly handles all inputs.
fn test_new() {
    {
        let threshold = 0.85258248;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), threshold);
    }
    {
        let threshold = 0.0;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), threshold);
    }
    {
        let threshold = -0.0;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 0.0);
    }
    {
        let threshold = -0.2940984;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 0.0);
    }
    {
        let threshold = f64::INFINITY;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 1.0);
    }
    {
        let threshold = f64::NEG_INFINITY;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 0.0);
    }
    {
        let threshold = f64::NAN;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 0.0);
    }
    {
        let threshold = f64::MAX * 2.0;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 1.0);
    }
    {
        let threshold = f64::MIN / 2.0;
        let dendrite = SimpleDendriteThreshold::new(threshold);
        assert_ulps_eq!(dendrite.threshold().value(), 0.0);
    }
}

#[test]
/// Tests if the function `is_similar` correctly detects similarity.
fn test_is_similar() {
    let dendrite_a = SimpleDendriteThreshold::random();
    let dendrite_b = SimpleDendriteThreshold::random();
    assert!(dendrite_a.is_similar(&dendrite_b));
    assert!(dendrite_b.is_similar(&dendrite_a));
}

#[test]
/// Tests if the function `cross_over` correctly recombines dendrites.
fn test_cross_over() {
    let dendrite_a: SimpleDendriteThreshold = thread_rng().gen();
    let dendrite_b: SimpleDendriteThreshold = thread_rng().gen();
    let dendrite_recombined = dendrite_a.cross_over(&dendrite_b);
    assert!(dendrite_recombined.threshold() <= dendrite_a.threshold() + dendrite_b.threshold());
}

#[test]
/// Tests if the function `detect` correctly detects neuron states.
fn test_detect() {
    let neuron = SimpleNeuron::new(0.5);
    {
        // Test lower threshold.
        let dendrite = SimpleDendriteThreshold::new(0.4);
        assert_eq!(dendrite.get_substrate_number(), 1);
        assert!(dendrite.detect(&[&neuron], Iteration::new()))
    }
    {
        // Test equal threshold.
        let dendrite = SimpleDendriteThreshold::new(0.5);
        assert_eq!(dendrite.get_substrate_number(), 1);
        assert!(dendrite.detect(&[&neuron], Iteration::new()))
    }
    {
        // Test higher threshold.
        let dendrite = SimpleDendriteThreshold::new(0.6);
        assert_eq!(dendrite.get_substrate_number(), 1);
        assert!(!dendrite.detect(&[&neuron], Iteration::new()))
    }
}
