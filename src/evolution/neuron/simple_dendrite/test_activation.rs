use super::*;

#[test]
/// Tests if the function `new` correctly handles all inputs.
fn test_new() {
    {
        let weight = 0.85258248;
        let dendrite = SimpleDendriteActivationPotential::new(weight, false);
        assert_ulps_eq!(dendrite.weight().value(), weight);
        assert!(!dendrite.is_inhibitory());
    }
    {
        let weight = 0.0;
        let dendrite = SimpleDendriteActivationPotential::new(weight, true);
        assert_ulps_eq!(dendrite.weight().value(), weight);
        assert!(dendrite.is_inhibitory());
    }
    {
        let weight = -0.0;
        let dendrite = SimpleDendriteActivationPotential::new(weight, false);
        assert_ulps_eq!(dendrite.weight().value(), 0.0);
        assert!(!dendrite.is_inhibitory());
    }
    {
        let weight = -0.2940984;
        let dendrite = SimpleDendriteActivationPotential::new(weight, false);
        assert_ulps_eq!(dendrite.weight().value(), 0.0);
        assert!(!dendrite.is_inhibitory());
    }
    {
        let weight = f64::INFINITY;
        let dendrite = SimpleDendriteActivationPotential::new(weight, true);
        assert_ulps_eq!(dendrite.weight().value(), 1.0);
        assert!(dendrite.is_inhibitory());
    }
    {
        let weight = f64::NEG_INFINITY;
        let dendrite = SimpleDendriteActivationPotential::new(weight, true);
        assert_ulps_eq!(dendrite.weight().value(), 0.0);
        assert!(dendrite.is_inhibitory());
    }
    {
        let weight = f64::NAN;
        let dendrite = SimpleDendriteActivationPotential::new(weight, false);
        assert_ulps_eq!(dendrite.weight().value(), 0.0);
        assert!(!dendrite.is_inhibitory());
    }
    {
        let weight = f64::MAX * 2.0;
        let dendrite = SimpleDendriteActivationPotential::new(weight, false);
        assert_ulps_eq!(dendrite.weight().value(), 1.0);
        assert!(!dendrite.is_inhibitory());
    }
    {
        let weight = f64::MIN / 2.0;
        let dendrite = SimpleDendriteActivationPotential::new(weight, true);
        assert_ulps_eq!(dendrite.weight().value(), 0.0);
        assert!(dendrite.is_inhibitory());
    }
}

#[test]
/// Tests if the function `is_similar` correctly detects similarity.
fn test_is_similar() {
    let dendrite_a = SimpleDendriteActivationPotential::random();
    let dendrite_b = SimpleDendriteActivationPotential::random();
    assert!(dendrite_a.is_similar(&dendrite_b));
    assert!(dendrite_b.is_similar(&dendrite_a));
}

#[test]
/// Tests if the function `cross_over` correctly recombines dendrites.
fn test_cross_over() {
    let dendrite_a: SimpleDendriteActivationPotential = thread_rng().gen();
    let dendrite_b: SimpleDendriteActivationPotential = thread_rng().gen();
    let dendrite_recombined = dendrite_a.cross_over(&dendrite_b);
    assert!(dendrite_recombined.weight() <= dendrite_a.weight() + dendrite_b.weight());
    assert!(
        (dendrite_recombined.is_inhibitory() == dendrite_a.is_inhibitory())
            || (dendrite_recombined.is_inhibitory() == dendrite_b.is_inhibitory())
    )
}

#[test]
/// Tests if the function `react` correctly influences neurons when non inhibitory.
fn test_react_non_inhibitory() {
    let neuron = SimpleNeuron::new(0.5);
    {
        // No addition.
        let dendrite = SimpleDendriteActivationPotential::new(0.0, false);
        assert_eq!(dendrite.get_educt_number(), 1);
        assert_eq!(dendrite.get_product_number(), 1);
        let products = dendrite.react(&[&neuron], Iteration::new());
        assert_eq!(products.len(), 1);
        assert_ulps_eq!(products[0].base_potential().value(), 0.5);
        assert_ulps_eq!(products[0].current_potential().value(), 0.5);
    }
    {
        // Normal addition.
        let dendrite = SimpleDendriteActivationPotential::new(0.2, false);
        assert_eq!(dendrite.get_educt_number(), 1);
        assert_eq!(dendrite.get_product_number(), 1);
        let products = dendrite.react(&[&neuron], Iteration::new());
        assert_eq!(products.len(), 1);
        assert_ulps_eq!(products[0].base_potential().value(), 0.5);
        assert_ulps_eq!(products[0].current_potential().value(), 0.7);
    }
    {
        // Overflow.
        let dendrite = SimpleDendriteActivationPotential::new(0.7, false);
        assert_eq!(dendrite.get_educt_number(), 1);
        assert_eq!(dendrite.get_product_number(), 1);
        let products = dendrite.react(&[&neuron], Iteration::new());
        assert_eq!(products.len(), 1);
        assert_ulps_eq!(products[0].base_potential().value(), 0.5);
        assert_ulps_eq!(products[0].current_potential().value(), 1.0);
    }
}

#[test]
/// Tests if the function `react` correctly influences neurons when inhibitory.
fn test_react_inhibitory() {
    let neuron = SimpleNeuron::new(0.5);
    {
        // No subtraction.
        let dendrite = SimpleDendriteActivationPotential::new(0.0, true);
        assert_eq!(dendrite.get_educt_number(), 1);
        assert_eq!(dendrite.get_product_number(), 1);
        let products = dendrite.react(&[&neuron], Iteration::new());
        assert_eq!(products.len(), 1);
        assert_ulps_eq!(products[0].base_potential().value(), 0.5);
        assert_ulps_eq!(products[0].current_potential().value(), 0.5);
    }
    {
        // Normal subtraction.
        let dendrite = SimpleDendriteActivationPotential::new(0.2, true);
        assert_eq!(dendrite.get_educt_number(), 1);
        assert_eq!(dendrite.get_product_number(), 1);
        let products = dendrite.react(&[&neuron], Iteration::new());
        assert_eq!(products.len(), 1);
        assert_ulps_eq!(products[0].base_potential().value(), 0.5);
        assert_ulps_eq!(products[0].current_potential().value(), 0.3);
    }
    {
        // Underflow.
        let dendrite = SimpleDendriteActivationPotential::new(0.7, true);
        assert_eq!(dendrite.get_educt_number(), 1);
        assert_eq!(dendrite.get_product_number(), 1);
        let products = dendrite.react(&[&neuron], Iteration::new());
        assert_eq!(products.len(), 1);
        assert_ulps_eq!(products[0].base_potential().value(), 0.5);
        assert_ulps_eq!(products[0].current_potential().value(), 0.0);
    }
}
