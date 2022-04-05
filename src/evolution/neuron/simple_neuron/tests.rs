use rand::thread_rng;

use super::*;

const POTENTIAL_HALFLIFE_TIME: f64 = 16.0;

#[test]
/// Tests if the function `new` correctly handles all inputs.
fn test_new() {
    {
        let potential = 0.85258248;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), potential);
        assert_ulps_eq!(neuron.current_potential().value(), potential);
    }
    {
        let potential = 0.0;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), potential);
        assert_ulps_eq!(neuron.current_potential().value(), potential);
    }
    {
        let potential = -0.0;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = -0.2940984;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = f64::INFINITY;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), 1.0);
        assert_ulps_eq!(neuron.current_potential().value(), 1.0);
    }
    {
        let potential = f64::NEG_INFINITY;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = f64::NAN;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), 0.0);
        assert_ulps_eq!(neuron.current_potential().value(), 0.0);
    }
    {
        let potential = f64::MAX * 2.0;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        assert_ulps_eq!(neuron.base_potential().value(), 1.0);
        assert_ulps_eq!(neuron.current_potential().value(), 1.0);
    }
    {
        let potential = f64::MIN / 2.0;
        let neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
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
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), new_current_potential);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = 0.0;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), new_current_potential);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = -0.0;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = -0.9046406490;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::INFINITY;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 1.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::NEG_INFINITY;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::NAN;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 0.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::MAX * 2.0;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
        let offspring_neuron = base_neuron.with_new_current_potential(new_current_potential);
        assert_ulps_eq!(offspring_neuron.base_potential().value(), potential);
        assert_ulps_eq!(offspring_neuron.current_potential().value(), 1.0);
    }
    {
        let potential = 0.85258248;
        let new_current_potential = f64::MIN / 2.0;
        let base_neuron = SimpleNeuron::new(potential, POTENTIAL_HALFLIFE_TIME);
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
    let neuron_a = SimpleNeuron::new(potential_a, POTENTIAL_HALFLIFE_TIME);
    let neuron_b = SimpleNeuron::new(potential_b, POTENTIAL_HALFLIFE_TIME);
    let neuron_add_a_to_b = neuron_a + neuron_b;
    let neuron_add_b_to_a = neuron_b + neuron_a;
    assert_ulps_eq!(neuron_add_a_to_b.base_potential().value(), potential_a);
    assert_ulps_eq!(neuron_add_a_to_b.current_potential().value(), potential_a + potential_b);
    assert_ulps_eq!(neuron_add_b_to_a.base_potential().value(), potential_b);
    assert_ulps_eq!(neuron_add_b_to_a.current_potential().value(), potential_a + potential_b);
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
    assert_eq!(neuron_recombined.base_potential(), neuron_recombined.current_potential());
    assert!(neuron_recombined <= neuron_a + neuron_b);
    assert!(neuron_recombined.potential_halflife_time() >= POTENTIAL_HALFLIFE_TIME_MINIMUM);
    assert!(neuron_recombined.potential_halflife_time() <= POTENTIAL_HALFLIFE_TIME_MAXIMUM);
}

#[test]
/// Tests if the function `update_value` correctly updates neuron values.
fn test_update_value() {
    {
        // Test same base and current potential.
        let potential: Nlbf64 = 0.4336654.into();
        let neuron_original: SimpleNeuron = SimpleNeuron {
            current_potential: potential,
            base_potential: potential,
            potential_halflife_time: POTENTIAL_HALFLIFE_TIME,
        };
        let mut neuron_updated: SimpleNeuron = neuron_original.clone();
        assert_eq!(neuron_original, neuron_updated);
        neuron_updated.update_value(5869);
        assert_eq!(neuron_original, neuron_updated);
    }
    {
        // Test positive difference between base and current potential.
        let base_potential: Nlbf64 = 0.4.into();
        let current_potential: Nlbf64 = 0.9.into();
        let neuron_original: SimpleNeuron = SimpleNeuron {
            current_potential,
            base_potential,
            potential_halflife_time: POTENTIAL_HALFLIFE_TIME,
        };
        let mut neuron_updated: SimpleNeuron = neuron_original.clone();
        assert_eq!(neuron_original, neuron_updated);
        neuron_updated.update_value(POTENTIAL_HALFLIFE_TIME as i32);
        assert_ne!(neuron_original, neuron_updated);
        assert_eq!(neuron_updated.base_potential(), neuron_original.base_potential());
        assert_ulps_eq!(neuron_updated.current_potential().value(), 0.65);
        neuron_updated.update_value((POTENTIAL_HALFLIFE_TIME * 1.5) as i32);
        assert_ne!(neuron_original, neuron_updated);
        assert_eq!(neuron_updated.base_potential(), neuron_original.base_potential());
        assert_ulps_eq!(neuron_updated.current_potential().value(), 0.48838834764831845);
    }
    {
        // Test negative difference between base and current potential.
        let base_potential: Nlbf64 = 0.9.into();
        let current_potential: Nlbf64 = 0.1.into();
        let neuron_original: SimpleNeuron = SimpleNeuron {
            current_potential,
            base_potential,
            potential_halflife_time: POTENTIAL_HALFLIFE_TIME,
        };
        let mut neuron_updated: SimpleNeuron = neuron_original.clone();
        assert_eq!(neuron_original, neuron_updated);
        neuron_updated.update_value((POTENTIAL_HALFLIFE_TIME as i32) * 2);
        assert_ne!(neuron_original, neuron_updated);
        assert_eq!(neuron_updated.base_potential(), neuron_original.base_potential());
        assert_ulps_eq!(neuron_updated.current_potential().value(), 0.7);
        neuron_updated.update_value((POTENTIAL_HALFLIFE_TIME * 1.5) as i32);
        assert_ne!(neuron_original, neuron_updated);
        assert_eq!(neuron_updated.base_potential(), neuron_original.base_potential());
        assert_ulps_eq!(neuron_updated.current_potential().value(), 0.8292893218813453);
    }
    {
        // Test no time passed.
        let base_potential: Nlbf64 = 0.5.into();
        let current_potential: Nlbf64 = 0.4.into();
        let neuron_original: SimpleNeuron = SimpleNeuron {
            current_potential,
            base_potential,
            potential_halflife_time: POTENTIAL_HALFLIFE_TIME,
        };
        let mut neuron_updated: SimpleNeuron = neuron_original.clone();
        assert_eq!(neuron_original, neuron_updated);
        neuron_updated.update_value(0);
        assert_eq!(neuron_original, neuron_updated);
    }
}

#[test]
/// Tests if the function `potential_half_time` correctly caps the underlying value.
fn test_potential_half_time() {
    {
        let neuron: SimpleNeuron = SimpleNeuron::new(0.5, POTENTIAL_HALFLIFE_TIME_MINIMUM - 100.0);
        assert!(neuron.potential_halflife_time() >= POTENTIAL_HALFLIFE_TIME_MINIMUM);
        assert!(neuron.potential_halflife_time() <= POTENTIAL_HALFLIFE_TIME_MAXIMUM);
    }
    {
        let neuron: SimpleNeuron = SimpleNeuron::new(0.5, 0.0);
        assert!(neuron.potential_halflife_time() >= POTENTIAL_HALFLIFE_TIME_MINIMUM);
        assert!(neuron.potential_halflife_time() <= POTENTIAL_HALFLIFE_TIME_MAXIMUM);
    }
    {
        let neuron: SimpleNeuron = SimpleNeuron::new(0.5, POTENTIAL_HALFLIFE_TIME_MAXIMUM + 100.0);
        assert!(neuron.potential_halflife_time() >= POTENTIAL_HALFLIFE_TIME_MINIMUM);
        assert!(neuron.potential_halflife_time() <= POTENTIAL_HALFLIFE_TIME_MAXIMUM);
    }
}
