use super::*;

#[test]
/// Tests if the function `new` correctly handles all inputs.
fn test_new() {
    {
        let potential = 85258.8248;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), potential);
        assert_eq!(neuron.current_potential(), potential);
    } {
        let potential = 0.0;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), potential);
        assert_eq!(neuron.current_potential(), potential);
    } {
        let potential = -0.0;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    } {
        let potential = -29.40984;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    } {
        let potential = f64::INFINITY;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    } {
        let potential = f64::NEG_INFINITY;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    } {
        let potential = f64::NAN;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    } {
        let potential = f64::MAX * 2.0;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    } {
        let potential = f64::MIN / 2.0;
        let neuron = SimpleNeuron::new(potential);
        assert_eq!(neuron.base_potential(), 0.0);
        assert_eq!(neuron.current_potential(), 0.0);
    }
}

#[test]
/// Tests if the function `with_new_current_potential` correctly handles all inputs.
fn test_with_new_current_potential() {
    {
        let potential = 85258.8248;
        let new_current_potential = 350.9390; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), new_current_potential);
    } {
        let potential = 85258.8248;
        let new_current_potential = 0.0; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), new_current_potential);
    } {
        let potential = 85258.8248;
        let new_current_potential = -0.0; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    } {
        let potential = 85258.8248;
        let new_current_potential = -90464064.90; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    } {
        let potential = 85258.8248;
        let new_current_potential = f64::INFINITY; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    } {
        let potential = 85258.8248;
        let new_current_potential = f64::NEG_INFINITY; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    } {
        let potential = 85258.8248;
        let new_current_potential = f64::NAN; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    } {
        let potential = 85258.8248;
        let new_current_potential = f64::MAX * 2.0; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    } {
        let potential = 85258.8248;
        let new_current_potential = f64::MIN / 2.0; 
        let base_neuron = SimpleNeuron::new(potential);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_eq!(offspring_neuron.base_potential(), potential);
        assert_eq!(offspring_neuron.current_potential(), 0.0);
    }
}

#[test]
/// Tests if the function `add` correctly handles neuron addition.
fn test_add() {
    let potential_a = 12456789.0;
    let potential_b = 987650.4321;
    let neuron_a = SimpleNeuron::new(potential_a);
    let neuron_b = SimpleNeuron::new(potential_b);
    let neuron_add_a_to_b = neuron_a + neuron_b;
    let neuron_add_b_to_a = neuron_b + neuron_a;
    assert_eq!(neuron_add_a_to_b.base_potential(), potential_a);
    assert_ulps_eq!(neuron_add_a_to_b.current_potential(), potential_a + potential_b);
    assert_eq!(neuron_add_b_to_a.base_potential(), potential_b);
    assert_ulps_eq!(neuron_add_b_to_a.current_potential(), potential_a + potential_b);
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
    let neuron_a = SimpleNeuron::random();
    let neuron_b = SimpleNeuron::random();
    let neuron_recombined= neuron_a.cross_over(&neuron_b);
    assert_eq!(neuron_recombined.base_potential(), neuron_recombined.current_potential());
    let binary_potential_a = f64_to_binary(neuron_a.base_potential());
    let binary_potential_b = f64_to_binary(neuron_b.base_potential());
    let binary_potential_recombined = f64_to_binary(neuron_recombined.base_potential()); 
    for i in 0..binary_potential_recombined.len() {
        assert!(binary_potential_recombined[i] == binary_potential_a[i] || binary_potential_recombined[i] == binary_potential_b[i]);
    }
}