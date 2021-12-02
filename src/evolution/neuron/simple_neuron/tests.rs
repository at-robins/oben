use rand::thread_rng;

use super::*;

#[test]
/// Tests if the function `new` correctly handles all inputs.
fn test_new() {
    {
        let potential = 0.85258248;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), potential);
        assert_ulps_eq!(neuron.current_potential().value(), potential);
    }
    {
        let potential = 0.0;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), potential);
        assert_ulps_eq!(neuron.current_potential().value(), potential);
    }
    {
        let potential = -0.0;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = -0.2940984;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = f64::INFINITY;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 1.0);
        assert_ulps_eq!(neuron.current_potential().value(), 1.0);
    }
    {
        let potential = f64::NEG_INFINITY;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = f64::NAN;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = f64::MAX * 2.0;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 1.0);
        assert_ulps_eq!(neuron.current_potential().value(), 1.0);
    }
    {
        let potential = f64::MIN / 2.0;
        let neuron = SimpleNeuron::new(potential);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
}

#[test]
/// Tests if the function `with_new_current_potential` correctly handles all inputs.
fn test_with_new_current_potential() {
    {
        let potential = 0.85258248;
        let new_current_potential = 0.3509390;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(
            offspring_neuron.current_potential().value(),
            new_current_potential
        );
    }
    {
        let potential = 0.85258248;
        let new_current_potential = 0.0;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(
            offspring_neuron.current_potential().value(),
            new_current_potential
        );
    }
    {
        let potential = 0.85258248;
        let new_current_potential = -0.0;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = -0.9046406490;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::INFINITY;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 1.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::NEG_INFINITY;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::NAN;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::MAX * 2.0;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 1.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::MIN / 2.0;
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
}

#[test]
/// Tests if the function `add` correctly handles neuron addition.
fn test_add() {
    let potential_a = 0.12456789;
    let potential_b = 0.76504321;
    let neuron_a = SimpleNeuron::new(potential_a);
    let neuron_b = SimpleNeuron::new(potential_b);
    let neuron_add_a_to_b = neuron_a + neuron_b;
    let neuron_add_b_to_a = neuron_b + neuron_a;
    assert_ulps_eq!(neuron_add_a_to_b.base_potential().value(), potential_a);
    assert_ulps_eq!(
        neuron_add_a_to_b.current_potential().value(),
        potential_a + potential_b
    );
    assert_ulps_eq!(neuron_add_b_to_a.base_potential().value(), potential_b);
    assert_ulps_eq!(
        neuron_add_b_to_a.current_potential().value(),
        potential_a + potential_b
    );
}

#[test]
/// Tests if the function `is_similar` correctly detects similarity.
fn test_is_similar() {
    let neuron_a = SimpleNeuron::random();
    let neuron_b = SimpleNeuron::random();
    assert!(neuron_a.is_similar(&neuron_b));
    assert!(neuron_b.is_similar(&neuron_a));
}

#[test]
/// Tests if the function `cross_over` correctly recombines neurons.
fn test_cross_over() {
    let neuron_a: SimpleNeuron = thread_rng().gen();
    let neuron_b: SimpleNeuron = thread_rng().gen();
    let neuron_recombined = neuron_a.cross_over(&neuron_b);
    assert_eq!(
        neuron_recombined.base_potential(),
        neuron_recombined.current_potential()
    );
    assert!(neuron_recombined <= neuron_a + neuron_b);
}
