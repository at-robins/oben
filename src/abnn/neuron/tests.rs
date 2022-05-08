use super::*;

#[test]
/// Tests if the function `new` correctly handles all inputs.
fn test_new() {
    let neuron = Neuron::new();
    assert_ulps_eq!(neuron.value(), NEURON_BASE_VALUE);
    assert_eq!(neuron.last_value_update(), Iteration::new());
    assert!(neuron.dendrites.lock().is_empty());
}

#[test]
/// Tests if the function `set_value` correctly handles all inputs.
fn test_set_value() {
    {
        let new_value = 0.3509390;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), new_value);
    }
    {
        let new_value = 0.0;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), new_value);
    }
    {
        let new_value = -0.0;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 0.0);
    }
    {
        let new_value = -0.9046406490;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 0.0);
    }
    {
        let new_value = f64::INFINITY;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 1.0);
    }
    {
        let new_value = f64::NEG_INFINITY;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 0.0);
    }
    {
        let new_value = f64::NAN;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 0.0);
    }
    {
        let new_value = f64::MAX * 2.0;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 1.0);
    }
    {
        let new_value = f64::MIN / 2.0;
        let neuron = Neuron::new();
        neuron.set_value(new_value);
        assert_ulps_eq!(neuron.value(), 0.0);
    }
}

#[test]
/// Tests if the function `update_value` correctly updates neuron values.
fn test_update_value() {
    {
        // Test no difference from base value.
        let neuron = Neuron::new();
        assert_eq!(neuron.value(), NEURON_BASE_VALUE);
        neuron.update_value(time_at_point(16));
        assert_eq!(neuron.value(), NEURON_BASE_VALUE);
    }
    {
        // Test positive difference between base and current potential.
        let current_value = 0.9;
        let neuron = Neuron::new();
        neuron.set_value(current_value);
        assert_eq!(neuron.value(), current_value);
        neuron.update_value(time_at_point(NEURON_VALUE_HALFLIFE as usize));
        assert_ulps_eq!(
            neuron.value(),
            ((current_value - NEURON_BASE_VALUE) / 2.0) + NEURON_BASE_VALUE
        );
        neuron.update_value(time_at_point((2.0 * NEURON_VALUE_HALFLIFE) as usize));
        assert_ulps_eq!(
            neuron.value(),
            ((current_value - NEURON_BASE_VALUE) / 4.0) + NEURON_BASE_VALUE
        );
    }
    {
        // Test negative difference between base and current potential.
        let current_value = 0.0;
        let neuron = Neuron::new();
        neuron.set_value(current_value);
        assert_eq!(neuron.value(), current_value);
        neuron.update_value(time_at_point(NEURON_VALUE_HALFLIFE as usize));
        assert_ulps_eq!(
            neuron.value(),
            ((current_value - NEURON_BASE_VALUE) / 2.0) + NEURON_BASE_VALUE
        );
        neuron.update_value(time_at_point((2.0 * NEURON_VALUE_HALFLIFE) as usize));
        assert_ulps_eq!(
            neuron.value(),
            ((current_value - NEURON_BASE_VALUE) / 4.0) + NEURON_BASE_VALUE
        );
    }
    {
        // Test no time passed.
        let current_value = 0.9;
        let neuron = Neuron::new();
        neuron.set_value(current_value);
        assert_eq!(neuron.value(), current_value);
        neuron.update_value(time_at_point(0));
        assert_ulps_eq!(neuron.value(), current_value);
    }
}

#[test]
#[should_panic]
/// Tests if the function `update_value` correctly panics for invalid values.
fn test_update_value_panic() {
    let neuron = Neuron::new();
    neuron.update_value(time_at_point(5));
    neuron.update_value(time_at_point(3));
}

fn time_at_point(timepoint: usize) -> Iteration {
    let mut time = Iteration::new();
    if timepoint > 0 {
        for _ in 0..timepoint {
            time = time.increment();
        }
    }
    time
}
