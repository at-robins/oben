use rand::Rng;
use rand_distr::{Distribution, Standard};
use serde::{Deserialize, Serialize};

use crate::evolution::{
    binary::{as_u64, flip_random_bit, u64_to_binary},
    gene::CrossOver,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
/// A non linear, equally spaced floating point representation bound between
/// lower bound `L` and upper bound `U`.
///
/// # Warning 1
///
/// Due to the current limitation of not being able to use floating point types
/// as const generics the const generics constitue floating point numbers
/// by specifying a number and an exponent (power of 10) for upper and lower bound.
/// lower bound: `L * 10 ^ LE`
/// upper bound: `U * 10 ^ UE`
///
/// # Warning 2
///
/// If the lower bound is greater than or equal to the upper bound or any of the bounds is non-finite
/// the lower bound is returned for all operations.
pub struct Nlgbf64<const L: i64, const LE: i32, const U: i64, const UE: i32> {
    value: u64,
}

impl<const L: i64, const LE: i32, const U: i64, const UE: i32> Nlgbf64<L, LE, U, UE> {
    /// The minimum value of L.
    pub const MIN: Nlgbf64<L, LE, U, UE> = Nlgbf64 { value: u64::MIN };
    /// The maximum value of U.
    pub const MAX: Nlgbf64<L, LE, U, UE> = Nlgbf64 { value: u64::MAX };

    /// Returns the computed lower bound based on the supplied const generics.
    fn lower_bound() -> f64 {
        (L as f64) * 10.0f64.powi(LE)
    }

    /// Returns the computed upper bound based on the supplied const generics.
    fn upper_bound() -> f64 {
        (U as f64) * 10.0f64.powi(UE)
    }

    /// Retruns the value as floating point number.
    pub fn value(&self) -> f64 {
        if Self::lower_bound() >= Self::upper_bound() {
            Self::lower_bound()
        } else if self.value == 0 {
            Self::lower_bound()
        } else if self.value == u64::MAX {
            Self::upper_bound()
        } else {
            Self::lower_bound()
                + (Self::upper_bound() - Self::lower_bound()) * (self.value as f64)
                    / (u64::MAX as f64)
        }
    }

    /// Retruns the value as internal u64 representation.
    pub fn value_as_u64(&self) -> u64 {
        self.value
    }

    /// Flips a random bit of the [`Nlgbf64`] and returns the result.
    ///
    /// # Parameters
    ///
    /// * `base` - the number to mutate
    pub fn flip_random_bit(base: Self) -> Self {
        let binary_base = u64_to_binary(base.value);
        Self {
            value: as_u64(&flip_random_bit(&binary_base)),
        }
    }
}

impl<const L: i64, const LE: i32, const U: i64, const UE: i32> From<f64> for Nlgbf64<L, LE, U, UE> {
    fn from(value: f64) -> Self {
        if Self::lower_bound() >= Self::upper_bound() {
            Self::MIN
        } else if value <= Self::lower_bound() {
            Self::MIN
        } else if value >= Self::upper_bound() {
            Self::MAX
        } else {
            Self {
                value: (((value - Self::lower_bound())
                    / (Self::upper_bound() - Self::lower_bound()))
                    * u64::MAX as f64) as u64,
            }
        }
    }
}

impl<const L: i64, const LE: i32, const U: i64, const UE: i32> From<&f64> for Nlgbf64<L, LE, U, UE> {
    fn from(value: &f64) -> Self {
        Self::from(*value)
    }
}

impl<const L: i64, const LE: i32, const U: i64, const UE: i32> Default for Nlgbf64<L, LE, U, UE> {
    fn default() -> Self {
        Self { value: 0 }
    }
}

impl<const L: i64, const LE: i32, const U: i64, const UE: i32> Distribution<Nlgbf64<L, LE, U, UE>>
    for Standard
{
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Nlgbf64<L, LE, U, UE> {
        if Nlgbf64::<L, LE, U, UE>::lower_bound() >= Nlgbf64::<L, LE, U, UE>::upper_bound() {
            Nlgbf64::<L, LE, U, UE>::MIN
        } else {
            Nlgbf64 { value: rng.gen() }
        }
    }
}

impl<const L: i64, const LE: i32, const U: i64, const UE: i32> CrossOver for Nlgbf64<L, LE, U, UE> {
    fn is_similar(&self, _other: &Self) -> bool {
        true
    }

    fn cross_over(&self, other: &Self) -> Self {
        Nlgbf64 {
            value: as_u64(&u64_to_binary(self.value).cross_over(&u64_to_binary(other.value))),
        }
    }
}

#[cfg(test)]
mod tests {
    use rand::thread_rng;

    use super::*;

    #[test]
    /// Tests the `value` function of the `Nlgbf64` struct.
    fn test_nlbf64_value() {
        {
            assert_ulps_eq!(Nlgbf64::<0, 0, 1, 0> { value: 0 }.value(), 0.0);
            assert_ulps_eq!(Nlgbf64::<-3, 2, 1, 0> { value: 0 }.value(), -300.0);
            assert_ulps_eq!(Nlgbf64::<12345, -5, 1, 0> { value: 0 }.value(), 0.12345);
        }
        {
            assert_ulps_eq!(
                Nlgbf64::<0, 0, 1, 0> {
                    value: u64::MAX / 4,
                }
                .value(),
                0.25
            );
            assert_ulps_eq!(
                Nlgbf64::<-50, 0, 50, 0> {
                    value: u64::MAX / 4,
                }
                .value(),
                -25.0
            );
            assert_ulps_eq!(
                Nlgbf64::<12, -2, 24, -2> {
                    value: u64::MAX / 4,
                }
                .value(),
                0.15
            );
        }
        {
            assert_ulps_eq!(
                Nlgbf64::<0, 0, 1, 0> {
                    value: u64::MAX / 2,
                }
                .value(),
                0.5
            );
            assert_ulps_eq!(
                Nlgbf64::<-50, 0, 50, 0> {
                    value: u64::MAX / 2,
                }
                .value(),
                0.0
            );
            assert_ulps_eq!(
                Nlgbf64::<12, -2, 24, -2> {
                    value: u64::MAX / 2,
                }
                .value(),
                0.18
            );
        }
        {
            assert_ulps_eq!(
                Nlgbf64::<0, 0, 1, 0> {
                    value: (u64::MAX / 4) * 3,
                }
                .value(),
                0.75
            );
            assert_ulps_eq!(
                Nlgbf64::<-50, 0, 50, 0> {
                    value: (u64::MAX / 4) * 3,
                }
                .value(),
                25.0
            );
            assert_ulps_eq!(
                Nlgbf64::<12, -2, 24, -2> {
                    value: (u64::MAX / 4) * 3,
                }
                .value(),
                0.21
            );
        }
        {
            assert_ulps_eq!(Nlgbf64::<0, 0, 1, 0> { value: u64::MAX }.value(), 1.0);
            assert_ulps_eq!(Nlgbf64::<-50, 0, 50, 0> { value: u64::MAX }.value(), 50.0);
            assert_ulps_eq!(Nlgbf64::<12, -2, 24, -2> { value: u64::MAX }.value(), 0.24);
        }
        {
            assert_ulps_eq!(Nlgbf64::<1, 0, 0, 0> { value: u64::MAX }.value(), 1.0);
            assert_ulps_eq!(Nlgbf64::<50, 0, -50, 0> { value: u64::MAX }.value(), 50.0);
            assert_ulps_eq!(Nlgbf64::<24, -2, 12, -2> { value: u64::MAX }.value(), 0.24);
        }
    }

    #[test]
    /// Tests the `from` function for conversion of `f64` to `Nlgbf64`.
    fn test_nlbf64_from_f64() {
        {
            // Default conversion.
            let initial_value_x = 0.3584503953569;
            let initial_value_y = 12.345;
            let initial_value_z = 0.18;
            let x: Nlgbf64<0, 0, 1, 0> = initial_value_x.into();
            let y: Nlgbf64<-5, 1, 5, 1> = initial_value_y.into();
            let z: Nlgbf64<12, -2, 24, -2> = initial_value_z.into();
            assert_ulps_eq!(x.value(), initial_value_x);
            // Relative testing as the equal spacing interferes with
            // the ulps variant here.
            assert_relative_eq!(y.value(), initial_value_y, epsilon = 0.000000001);
            assert_ulps_eq!(z.value(), initial_value_z);
        }
        {
            // Smaller than lower bound.
            let x: Nlgbf64<0, 0, 1, 0> = (-0.3584503953569).into();
            let y: Nlgbf64<-5, 1, 5, 1> = (-51.234).into();
            let z: Nlgbf64<12, -2, 24, -2> = 0.0.into();
            assert_ulps_eq!(x.value(), 0.0);
            assert_ulps_eq!(y.value(), -50.0);
            assert_ulps_eq!(z.value(), 0.12);
        }
        {
            // Larger than upper bound.
            let x: Nlgbf64<0, 0, 1, 0> = 1.3584503953569.into();
            let y: Nlgbf64<-5, 1, 5, 1> = 492.1.into();
            let z: Nlgbf64<12, -2, 24, -2> = 0.3.into();
            assert_ulps_eq!(x.value(), 1.0);
            assert_ulps_eq!(y.value(), 50.0);
            assert_ulps_eq!(z.value(), 0.24);
        }
        {
            // NaN.
            let x: Nlgbf64<0, 0, 1, 0> = f64::NAN.into();
            let y: Nlgbf64<-5, 1, 5, 1> = f64::NAN.into();
            let z: Nlgbf64<12, -2, 24, -2> = f64::NAN.into();
            assert_ulps_eq!(x.value(), 0.0);
            assert_ulps_eq!(y.value(), -50.0);
            assert_ulps_eq!(z.value(), 0.12);
        }
        {
            // Incorrect bounds.
            let initial_value_x = 0.3584503953569;
            let initial_value_y = 12.345;
            let initial_value_z = 0.18;
            let x: Nlgbf64<1, 0, 0, 0> = initial_value_x.into();
            let y: Nlgbf64<5, 1, -5, 1> = initial_value_y.into();
            let z: Nlgbf64<24, -2, 12, -2> = initial_value_z.into();
            assert_ulps_eq!(x.value(), 1.0);
            assert_ulps_eq!(y.value(), 50.0);
            assert_ulps_eq!(z.value(), 0.24);
        }
    }

    #[test]
    /// Tests the `from` function for conversion of `f64` to `Nlgbf64`.
    fn test_nlbf64_from_f64_ref() {
        {
            // Default conversion.
            let initial_value_x = &0.3584503953569;
            let initial_value_y = &12.345;
            let initial_value_z = &0.18;
            let x: Nlgbf64<0, 0, 1, 0> = initial_value_x.into();
            let y: Nlgbf64<-5, 1, 5, 1> = initial_value_y.into();
            let z: Nlgbf64<12, -2, 24, -2> = initial_value_z.into();
            assert_ulps_eq!(x.value(), initial_value_x);
            // Relative testing as the equal spacing interferes with
            // the ulps variant here.
            assert_relative_eq!(y.value(), initial_value_y, epsilon = 0.000000001);
            assert_ulps_eq!(z.value(), initial_value_z);
        }
        {
            // Smaller than lower bound.
            let x: Nlgbf64<0, 0, 1, 0> = (&-0.3584503953569).into();
            let y: Nlgbf64<-5, 1, 5, 1> = (&-51.234).into();
            let z: Nlgbf64<12, -2, 24, -2> = (&0.0).into();
            assert_ulps_eq!(x.value(), 0.0);
            assert_ulps_eq!(y.value(), -50.0);
            assert_ulps_eq!(z.value(), 0.12);
        }
        {
            // Larger than upper bound.
            let x: Nlgbf64<0, 0, 1, 0> = (&1.3584503953569).into();
            let y: Nlgbf64<-5, 1, 5, 1> = (&492.1).into();
            let z: Nlgbf64<12, -2, 24, -2> = (&0.3).into();
            assert_ulps_eq!(x.value(), 1.0);
            assert_ulps_eq!(y.value(), 50.0);
            assert_ulps_eq!(z.value(), 0.24);
        }
        {
            // NaN.
            let x: Nlgbf64<0, 0, 1, 0> = (&f64::NAN).into();
            let y: Nlgbf64<-5, 1, 5, 1> = (&f64::NAN).into();
            let z: Nlgbf64<12, -2, 24, -2> = (&f64::NAN).into();
            assert_ulps_eq!(x.value(), 0.0);
            assert_ulps_eq!(y.value(), -50.0);
            assert_ulps_eq!(z.value(), 0.12);
        }
        {
            // Incorrect bounds.
            let initial_value_x = &0.3584503953569;
            let initial_value_y = &12.345;
            let initial_value_z = &0.18;
            let x: Nlgbf64<1, 0, 0, 0> = initial_value_x.into();
            let y: Nlgbf64<5, 1, -5, 1> = initial_value_y.into();
            let z: Nlgbf64<24, -2, 12, -2> = initial_value_z.into();
            assert_ulps_eq!(x.value(), 1.0);
            assert_ulps_eq!(y.value(), 50.0);
            assert_ulps_eq!(z.value(), 0.24);
        }
    }

    #[test]
    /// Tests if the function `is_similar` of the `CrossOver` trait can correctly detect similar
    /// substrates.
    fn test_nlbf64_cross_over_is_similar() {
        // Test similar substrates.
        let a: Nlgbf64<0, 0, 1, 0> = 0.1.into();
        let b: Nlgbf64<0, 0, 1, 0> = 0.6.into();
        assert!(a.is_similar(&b));
    }

    #[test]
    /// Tests if the function `cross_over` of the `CrossOver` trait can correctly recombine
    /// substrates.
    fn test_nlbf64_cross_over_cross_over() {
        let a: Nlgbf64<0, 0, 1, 0> = thread_rng().gen();
        let b: Nlgbf64<0, 0, 1, 0> = thread_rng().gen();
        let binary_a = u64_to_binary(a.value);
        let binary_b = u64_to_binary(b.value);
        let recombined = a.cross_over(&b);
        let binary_recombined = u64_to_binary(recombined.value);
        for i in 0..binary_a.len() {
            assert!(
                binary_recombined.get(i) == binary_a.get(i)
                    || binary_recombined.get(i) == binary_b.get(i)
            );
        }
    }

    #[test]
    /// Tests if the constant `MIN` of the `Nlgbf64` struct returns the correct value.
    fn test_nlbf64_min() {
        assert_ulps_eq!(Nlgbf64::<0, 0, 1, 0>::MIN.value(), 0.0);
        assert_ulps_eq!(Nlgbf64::<-5, 1, 5, 1>::MIN.value(), -50.0);
        assert_ulps_eq!(Nlgbf64::<12, -2, 24, -2>::MIN.value(), 0.12);
    }

    #[test]
    /// Tests if the constant `MAX` of the `Nlgbf64` struct returns the correct value.
    fn test_nlbf64_max() {
        assert_ulps_eq!(Nlgbf64::<0, 0, 1, 0>::MAX.value(), 1.0);
        assert_ulps_eq!(Nlgbf64::<-5, 1, 5, 1>::MAX.value(), 50.0);
        assert_ulps_eq!(Nlgbf64::<12, -2, 24, -2>::MAX.value(), 0.24);
    }

    #[test]
    /// Tests if the function `flip_random_bit` correctly mutates a single bit.
    fn test_nlbf64_flip_random_bit() {
        let original_value = 49495049;
        let original = Nlgbf64::<12, -2, 24, -2> {
            value: original_value,
        };
        let mutated = Nlgbf64::flip_random_bit(original);
        let difference = u64_to_binary(mutated.value ^ original.value);
        assert_eq!(difference.count_ones(), 1);
    }
}
