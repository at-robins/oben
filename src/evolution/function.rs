use std::{borrow::Borrow, fmt::Display};

use rand::Rng;
use serde::{Deserialize, Serialize};

use super::{
    binary::{as_f64, f64_to_binary, flip_random_bit},
    chemistry::Information,
    gene::CrossOver,
    helper::do_a_or_b,
};

pub const EXTERNAL_PARAMETER_DEFAULT: f64 = f64::NAN;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
/// A representation of a one-dimensional mathematical function.
pub enum MathematicalFunction {
    Absolute(Box<MathematicalFunction>),
    Add(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Constant(f64),
    Cosinus(Box<MathematicalFunction>),
    CosinusH(Box<MathematicalFunction>),
    Divide(Box<MathematicalFunction>, Box<MathematicalFunction>),
    ExternalParameter(usize),
    LimitLowerBound(Box<MathematicalFunction>, Box<MathematicalFunction>),
    LimitUpperBound(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Logarithm(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Maximum(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Minimum(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Multiply(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Power(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Sinus(Box<MathematicalFunction>),
    SinusH(Box<MathematicalFunction>),
    Subtract(Box<MathematicalFunction>, Box<MathematicalFunction>),
    Tangens(Box<MathematicalFunction>),
    TangensH(Box<MathematicalFunction>),
}

impl MathematicalFunction {
    /// Evaluates the function based on the specified external parameters.
    ///
    /// # Parameters
    ///
    /// * `external_parameters` - the specified external parameters to use during evaluation
    pub fn evaluate<V: Borrow<Vec<f64>>>(&self, external_parameters: V) -> f64 {
        match &self {
            MathematicalFunction::Absolute(value) => value.evaluate(external_parameters).abs(),
            MathematicalFunction::Add(param_a, param_b) => {
                param_a.evaluate(external_parameters.borrow())
                    + param_b.evaluate(external_parameters.borrow())
            },
            MathematicalFunction::Constant(value) => *value,
            MathematicalFunction::Cosinus(value) => value.evaluate(external_parameters).cos(),
            MathematicalFunction::CosinusH(value) => value.evaluate(external_parameters).cosh(),
            MathematicalFunction::Divide(param_a, param_b) => {
                param_a.evaluate(external_parameters.borrow())
                    / param_b.evaluate(external_parameters.borrow())
            },
            MathematicalFunction::ExternalParameter(identifier) => {
                Self::identifier_to_external_parameter(*identifier, external_parameters)
            },
            MathematicalFunction::LimitLowerBound(value, lower_bound) => {
                let val = value.evaluate(external_parameters.borrow());
                let bound = lower_bound.evaluate(external_parameters.borrow());
                if val < bound {
                    bound
                } else {
                    val
                }
            },
            MathematicalFunction::LimitUpperBound(value, upper_bound) => {
                let val = value.evaluate(external_parameters.borrow());
                let bound = upper_bound.evaluate(external_parameters.borrow());
                if val > bound {
                    bound
                } else {
                    val
                }
            },
            MathematicalFunction::Logarithm(value, base) => value
                .evaluate(external_parameters.borrow())
                .log(base.evaluate(external_parameters.borrow())),
            MathematicalFunction::Maximum(param_a, param_b) => {
                let a = param_a.evaluate(external_parameters.borrow());
                let b = param_b.evaluate(external_parameters.borrow());
                if a >= b {
                    a
                } else {
                    b
                }
            },
            MathematicalFunction::Minimum(param_a, param_b) => {
                let a = param_a.evaluate(external_parameters.borrow());
                let b = param_b.evaluate(external_parameters.borrow());
                if a <= b {
                    a
                } else {
                    b
                }
            },
            MathematicalFunction::Multiply(param_a, param_b) => {
                param_a.evaluate(external_parameters.borrow())
                    * param_b.evaluate(external_parameters.borrow())
            },
            MathematicalFunction::Power(value, exponent) => value
                .evaluate(external_parameters.borrow())
                .powf(exponent.evaluate(external_parameters.borrow())),
            MathematicalFunction::Sinus(value) => value.evaluate(external_parameters).sin(),
            MathematicalFunction::SinusH(value) => value.evaluate(external_parameters).sinh(),
            MathematicalFunction::Subtract(param_a, param_b) => {
                param_a.evaluate(external_parameters.borrow())
                    - param_b.evaluate(external_parameters.borrow())
            },
            MathematicalFunction::Tangens(value) => value.evaluate(external_parameters).tan(),
            MathematicalFunction::TangensH(value) => value.evaluate(external_parameters).tanh(),
        }
    }

    /// Performes the lookup of the internally saved external parameter identifier in the
    /// specified list of parameters.
    ///
    /// # Parameters
    ///
    /// * `identifier` - the identifier to perform the lookup for
    /// * `external_parameters` - the list to use for the lookup
    fn identifier_to_external_parameter<V: Borrow<Vec<f64>>>(
        identifier: usize,
        external_parameters: V,
    ) -> f64 {
        *(external_parameters
            .borrow()
            .get(identifier)
            .unwrap_or(&EXTERNAL_PARAMETER_DEFAULT))
    }

    /// Samples a mathematical function with a maximum depth as specified.
    ///
    /// # Parameters
    ///
    /// * `rng` - the random number generator to use for sampling
    /// * `max_depth` - the maximum number of function evaluation layers
    /// * `parameters` - the external parameter identifiers that can be used
    pub fn sample_with_max_depth<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        rng: &mut R,
        max_depth: usize,
        parameters: V,
    ) -> Self {
        Self::sample_with_max_depth_internal(rng, max_depth, parameters, 0)
    }

    fn sample_with_max_depth_internal<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        rng: &mut R,
        max_depth: usize,
        parameters: V,
        current_depth: usize,
    ) -> Self {
        if current_depth >= max_depth {
            match rng.gen_range(0usize..2) {
                0 => Self::sample_constant(rng),
                1 => Self::sample_external_parameter(rng, parameters),
                _ => panic!("Sampling out of bounds."),
            }
        } else {
            let next_depth = current_depth + 1;
            match rng.gen_range(0usize..19) {
                0 => Self::sample_constant(rng),
                1 => Self::sample_external_parameter(rng, parameters),
                2 => Self::Absolute(Self::sample_boxed(
                    rng,
                    max_depth,
                    parameters.borrow(),
                    next_depth,
                )),
                3 => Self::Add(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                4 => Self::Cosinus(Self::sample_boxed(
                    rng,
                    max_depth,
                    parameters.borrow(),
                    next_depth,
                )),
                5 => Self::CosinusH(Self::sample_boxed(
                    rng,
                    max_depth,
                    parameters.borrow(),
                    next_depth,
                )),
                6 => Self::Divide(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                7 => Self::LimitLowerBound(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                8 => Self::LimitUpperBound(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                9 => Self::Logarithm(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                10 => Self::Maximum(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                11 => Self::Minimum(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                12 => Self::Multiply(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                13 => Self::Power(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                14 => {
                    Self::Sinus(Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth))
                },
                15 => Self::SinusH(Self::sample_boxed(
                    rng,
                    max_depth,
                    parameters.borrow(),
                    next_depth,
                )),
                16 => Self::Subtract(
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                    Self::sample_boxed(rng, max_depth, parameters.borrow(), next_depth),
                ),
                17 => Self::Tangens(Self::sample_boxed(
                    rng,
                    max_depth,
                    parameters.borrow(),
                    next_depth,
                )),
                18 => Self::TangensH(Self::sample_boxed(
                    rng,
                    max_depth,
                    parameters.borrow(),
                    next_depth,
                )),
                _ => panic!("Sampling out of bounds."),
            }
        }
    }

    fn sample_boxed<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        rng: &mut R,
        max_depth: usize,
        parameters: V,
        current_depth: usize,
    ) -> Box<Self> {
        Box::new(Self::sample_with_max_depth_internal(rng, max_depth, parameters, current_depth))
    }

    fn sample_constant<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::Constant(f64::from_be_bytes(rng.gen::<u64>().to_be_bytes()))
    }

    fn sample_external_parameter<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        rng: &mut R,
        parameters: V,
    ) -> Self {
        let parameters: &Vec<usize> = parameters.borrow();
        if parameters.is_empty() {
            Self::sample_constant(rng)
        } else {
            Self::ExternalParameter(*parameters.get(rng.gen_range(0..parameters.len())).unwrap())
        }
    }

    /// Returns a mutated version of the `MathematicalFunction`.
    /// 
    /// # Parameters
    /// 
    /// * `rng` - the pseudorandom number generator to use
    /// * `external_parameters` - the parameter identifiers supplied by external sources 
    /// * `mutation_chance` - the chance of a mutation at an given branch of the function
    pub fn mutate<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        &self,
        rng: &mut R,
        external_parameters: V,
        mutation_chance: f64,
    ) -> Self {
        if mutation_chance < 0.0 || mutation_chance > 1.0 {
            panic!("The chance {} is out of bounds.", mutation_chance);
        }
        if rng.gen_range(0.0..=1.0f64) <= mutation_chance {
            match rng.gen_range(0..2) {
                0 => self.mutate_unwrap_function_or_change_value(rng, external_parameters),
                1 => self.mutate_wrap_function(rng, external_parameters),
                _ => panic!("Sampling out of bounds."),
            }
        } else {
            match &self {
                MathematicalFunction::Absolute(value) => Self::Absolute(Box::new(value.mutate(
                    rng,
                    external_parameters,
                    mutation_chance,
                ))),
                MathematicalFunction::Add(augend, addend) => Self::Add(
                    Box::new(augend.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(addend.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Constant(_) => self.clone(),
                MathematicalFunction::Cosinus(value) => {
                    Self::Cosinus(Box::new(value.mutate(rng, external_parameters, mutation_chance)))
                },
                MathematicalFunction::CosinusH(value) => Self::CosinusH(Box::new(value.mutate(
                    rng,
                    external_parameters,
                    mutation_chance,
                ))),
                MathematicalFunction::Divide(dividend, divisor) => Self::Divide(
                    Box::new(dividend.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(divisor.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::ExternalParameter(_) => self.clone(),
                MathematicalFunction::LimitLowerBound(value, lower_bound) => Self::LimitLowerBound(
                    Box::new(value.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(lower_bound.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::LimitUpperBound(value, upper_bound) => Self::LimitUpperBound(
                    Box::new(value.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(upper_bound.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Logarithm(value, base) => Self::Logarithm(
                    Box::new(value.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(base.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Maximum(a, b) => Self::Maximum(
                    Box::new(a.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(b.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Minimum(a, b) => Self::Minimum(
                    Box::new(a.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(b.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Multiply(multiplier, multiplicand) => Self::Multiply(
                    Box::new(multiplier.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(multiplicand.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Power(value, exponent) => Self::Power(
                    Box::new(value.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(exponent.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Sinus(value) => {
                    Self::Sinus(Box::new(value.mutate(rng, external_parameters, mutation_chance)))
                },
                MathematicalFunction::SinusH(value) => {
                    Self::SinusH(Box::new(value.mutate(rng, external_parameters, mutation_chance)))
                },
                MathematicalFunction::Subtract(minuend, subtrahend) => Self::Subtract(
                    Box::new(minuend.mutate(rng, external_parameters.borrow(), mutation_chance)),
                    Box::new(subtrahend.mutate(rng, external_parameters.borrow(), mutation_chance)),
                ),
                MathematicalFunction::Tangens(value) => {
                    Self::Tangens(Box::new(value.mutate(rng, external_parameters, mutation_chance)))
                },
                MathematicalFunction::TangensH(value) => Self::TangensH(Box::new(value.mutate(
                    rng,
                    external_parameters,
                    mutation_chance,
                ))),
            }
        }
    }

    fn mutate_unwrap_function_or_change_value<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        &self,
        rng: &mut R,
        external_parameters: V,
    ) -> Self {
        match &self {
            MathematicalFunction::Absolute(value) => *value.clone(),
            MathematicalFunction::Add(augend, addend) => {
                do_a_or_b(|| *augend.clone(), || *addend.clone())
            },
            MathematicalFunction::Constant(value) => {
                MathematicalFunction::Constant(as_f64(&flip_random_bit(&f64_to_binary(*value))))
            },
            MathematicalFunction::Cosinus(value) => *value.clone(),
            MathematicalFunction::CosinusH(value) => *value.clone(),
            MathematicalFunction::Divide(dividend, divisor) => {
                do_a_or_b(|| *dividend.clone(), || *divisor.clone())
            },
            MathematicalFunction::ExternalParameter(_) => {
                MathematicalFunction::sample_external_parameter(rng, external_parameters)
            },
            MathematicalFunction::LimitLowerBound(value, lower_bound) => {
                do_a_or_b(|| *value.clone(), || *lower_bound.clone())
            },
            MathematicalFunction::LimitUpperBound(value, upper_bound) => {
                do_a_or_b(|| *value.clone(), || *upper_bound.clone())
            },
            MathematicalFunction::Logarithm(value, base) => {
                do_a_or_b(|| *value.clone(), || *base.clone())
            },
            MathematicalFunction::Maximum(a, b) => do_a_or_b(|| *a.clone(), || *b.clone()),
            MathematicalFunction::Minimum(a, b) => do_a_or_b(|| *a.clone(), || *b.clone()),
            MathematicalFunction::Multiply(multiplier, multiplicand) => {
                do_a_or_b(|| *multiplier.clone(), || *multiplicand.clone())
            },
            MathematicalFunction::Power(value, exponent) => {
                do_a_or_b(|| *value.clone(), || *exponent.clone())
            },
            MathematicalFunction::Sinus(value) => *value.clone(),
            MathematicalFunction::SinusH(value) => *value.clone(),
            MathematicalFunction::Subtract(minuend, subtrahend) => {
                do_a_or_b(|| *minuend.clone(), || *subtrahend.clone())
            },
            MathematicalFunction::Tangens(value) => *value.clone(),
            MathematicalFunction::TangensH(value) => *value.clone(),
        }
    }

    fn mutate_wrap_function<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        &self,
        rng: &mut R,
        parameters: V,
    ) -> Self {
        match rng.gen_range(0usize..17) {
            0 => Self::Absolute(Box::new(self.clone())),
            1 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Add(left, right)
            },
            2 => Self::Cosinus(Box::new(self.clone())),
            3 => Self::CosinusH(Box::new(self.clone())),
            4 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Divide(left, right)
            },
            5 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::LimitLowerBound(left, right)
            },
            6 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::LimitUpperBound(left, right)
            },
            7 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Logarithm(left, right)
            },
            8 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Maximum(left, right)
            },
            9 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Minimum(left, right)
            },
            10 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Multiply(left, right)
            },
            11 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Power(left, right)
            },
            12 => Self::Sinus(Box::new(self.clone())),
            13 => Self::SinusH(Box::new(self.clone())),
            14 => {
                let (left, right) = self.wrap_two_arguments(rng, parameters);
                Self::Subtract(left, right)
            },
            15 => Self::Tangens(Box::new(self.clone())),
            16 => Self::TangensH(Box::new(self.clone())),
            _ => panic!("Sampling out of bounds."),
        }
    }

    fn wrap_two_arguments<R: Rng + ?Sized, V: Borrow<Vec<usize>>>(
        &self,
        rng: &mut R,
        parameters: V,
    ) -> (Box<MathematicalFunction>, Box<MathematicalFunction>) {
        if rng.gen() {
            (Box::new(self.clone()), Self::sample_boxed(rng, 0, parameters.borrow(), 0))
        } else {
            (Self::sample_boxed(rng, 0, parameters.borrow(), 0), Box::new(self.clone()))
        }
    }
}

impl CrossOver for MathematicalFunction {
    fn is_similar(&self, other: &Self) -> bool {
        match (self, other) {
            (MathematicalFunction::Absolute(_), MathematicalFunction::Absolute(_)) => true,
            (MathematicalFunction::Add(_, _), MathematicalFunction::Add(_, _)) => true,
            (MathematicalFunction::Constant(_), MathematicalFunction::Constant(_)) => true,
            (MathematicalFunction::Cosinus(_), MathematicalFunction::Cosinus(_)) => true,
            (MathematicalFunction::CosinusH(_), MathematicalFunction::CosinusH(_)) => true,
            (MathematicalFunction::Divide(_, _), MathematicalFunction::Divide(_, _)) => true,
            (
                MathematicalFunction::ExternalParameter(_),
                MathematicalFunction::ExternalParameter(_),
            ) => true,
            (
                MathematicalFunction::LimitLowerBound(_, _),
                MathematicalFunction::LimitLowerBound(_, _),
            ) => true,
            (
                MathematicalFunction::LimitUpperBound(_, _),
                MathematicalFunction::LimitUpperBound(_, _),
            ) => true,
            (MathematicalFunction::Logarithm(_, _), MathematicalFunction::Logarithm(_, _)) => true,
            (MathematicalFunction::Maximum(_, _), MathematicalFunction::Maximum(_, _)) => true,
            (MathematicalFunction::Minimum(_, _), MathematicalFunction::Minimum(_, _)) => true,
            (MathematicalFunction::Multiply(_, _), MathematicalFunction::Multiply(_, _)) => true,
            (MathematicalFunction::Power(_, _), MathematicalFunction::Power(_, _)) => true,
            (MathematicalFunction::Sinus(_), MathematicalFunction::Sinus(_)) => true,
            (MathematicalFunction::SinusH(_), MathematicalFunction::SinusH(_)) => true,
            (MathematicalFunction::Subtract(_, _), MathematicalFunction::Subtract(_, _)) => true,
            (MathematicalFunction::Tangens(_), MathematicalFunction::Tangens(_)) => true,
            (MathematicalFunction::TangensH(_), MathematicalFunction::TangensH(_)) => true,
            _ => false,
        }
    }

    fn cross_over(&self, other: &Self) -> Self {
        match (self, other) {
            (MathematicalFunction::Absolute(a), MathematicalFunction::Absolute(b)) => {
                MathematicalFunction::Absolute(Box::new(a.cross_over(b)))
            },
            (
                MathematicalFunction::Add(augend_a, addend_a),
                MathematicalFunction::Add(augend_b, addend_b),
            ) => MathematicalFunction::Add(
                Box::new(augend_a.cross_over(augend_b)),
                Box::new(addend_a.cross_over(addend_b)),
            ),
            (MathematicalFunction::Constant(a), MathematicalFunction::Constant(b)) => {
                MathematicalFunction::Constant(as_f64(
                    &f64_to_binary(*a).cross_over(&f64_to_binary(*b)),
                ))
            },
            (MathematicalFunction::Cosinus(a), MathematicalFunction::Cosinus(b)) => {
                MathematicalFunction::Cosinus(Box::new(a.cross_over(b)))
            },
            (MathematicalFunction::CosinusH(a), MathematicalFunction::CosinusH(b)) => {
                MathematicalFunction::CosinusH(Box::new(a.cross_over(b)))
            },
            (
                MathematicalFunction::Divide(dividend_a, divisor_a),
                MathematicalFunction::Divide(dividend_b, divisor_b),
            ) => MathematicalFunction::Divide(
                Box::new(dividend_a.cross_over(dividend_b)),
                Box::new(divisor_a.cross_over(divisor_b)),
            ),
            (
                MathematicalFunction::LimitLowerBound(value_a, bound_a),
                MathematicalFunction::LimitLowerBound(value_b, bound_b),
            ) => MathematicalFunction::LimitLowerBound(
                Box::new(value_a.cross_over(value_b)),
                Box::new(bound_a.cross_over(bound_b)),
            ),
            (
                MathematicalFunction::LimitUpperBound(value_a, bound_a),
                MathematicalFunction::LimitUpperBound(value_b, bound_b),
            ) => MathematicalFunction::LimitUpperBound(
                Box::new(value_a.cross_over(value_b)),
                Box::new(bound_a.cross_over(bound_b)),
            ),
            (
                MathematicalFunction::Logarithm(value_a, base_a),
                MathematicalFunction::Logarithm(value_b, base_b),
            ) => MathematicalFunction::Logarithm(
                Box::new(value_a.cross_over(value_b)),
                Box::new(base_a.cross_over(base_b)),
            ),
            (
                MathematicalFunction::Maximum(left_a, right_a),
                MathematicalFunction::Maximum(left_b, right_b),
            ) => MathematicalFunction::Maximum(
                Box::new(left_a.cross_over(left_b)),
                Box::new(right_a.cross_over(right_b)),
            ),
            (
                MathematicalFunction::Minimum(left_a, right_a),
                MathematicalFunction::Minimum(left_b, right_b),
            ) => MathematicalFunction::Minimum(
                Box::new(left_a.cross_over(left_b)),
                Box::new(right_a.cross_over(right_b)),
            ),
            (
                MathematicalFunction::Multiply(multiplier_a, multiplicand_a),
                MathematicalFunction::Multiply(multiplier_b, multiplicand_b),
            ) => MathematicalFunction::Multiply(
                Box::new(multiplier_a.cross_over(multiplier_b)),
                Box::new(multiplicand_a.cross_over(multiplicand_b)),
            ),
            (
                MathematicalFunction::Power(value_a, exponent_a),
                MathematicalFunction::Power(value_b, exponent_b),
            ) => MathematicalFunction::Power(
                Box::new(value_a.cross_over(value_b)),
                Box::new(exponent_a.cross_over(exponent_b)),
            ),
            (MathematicalFunction::Sinus(a), MathematicalFunction::Sinus(b)) => {
                MathematicalFunction::Sinus(Box::new(a.cross_over(b)))
            },
            (MathematicalFunction::SinusH(a), MathematicalFunction::SinusH(b)) => {
                MathematicalFunction::SinusH(Box::new(a.cross_over(b)))
            },
            (
                MathematicalFunction::Subtract(minuend_a, subtrahend_a),
                MathematicalFunction::Subtract(minuend_b, subtrahend_b),
            ) => MathematicalFunction::Subtract(
                Box::new(minuend_a.cross_over(minuend_b)),
                Box::new(subtrahend_a.cross_over(subtrahend_b)),
            ),
            (MathematicalFunction::Tangens(a), MathematicalFunction::Tangens(b)) => {
                MathematicalFunction::Tangens(Box::new(a.cross_over(b)))
            },
            (MathematicalFunction::TangensH(a), MathematicalFunction::TangensH(b)) => {
                MathematicalFunction::TangensH(Box::new(a.cross_over(b)))
            },
            (a, b) => do_a_or_b(|| a.clone(), || b.clone()),
        }
    }
}

impl Display for MathematicalFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self {
            MathematicalFunction::Absolute(value) => write!(f, "abs({})", value),
            MathematicalFunction::Add(augend, addend) => {
                write!(f, "({}) + ({})", augend, addend)
            },
            MathematicalFunction::Constant(value) => write!(f, "{}", value),
            MathematicalFunction::Cosinus(value) => write!(f, "cos({})", value),
            MathematicalFunction::CosinusH(value) => write!(f, "cosh({})", value),
            MathematicalFunction::Divide(dividend, divisor) => {
                write!(f, "({}) / ({})", dividend, divisor)
            },
            MathematicalFunction::ExternalParameter(identifier) => {
                write!(f, "param({})", identifier)
            },
            MathematicalFunction::LimitLowerBound(value, lower_bound) => {
                write!(f, "limit_lower({}, {})", value, lower_bound)
            },
            MathematicalFunction::LimitUpperBound(value, upper_bound) => {
                write!(f, "limit_upper({}, {})", value, upper_bound)
            },
            MathematicalFunction::Logarithm(value, base) => {
                write!(f, "log({}, {})", value, base)
            },
            MathematicalFunction::Maximum(a, b) => {
                write!(f, "max({}, {})", a, b)
            },
            MathematicalFunction::Minimum(a, b) => {
                write!(f, "min({}, {})", a, b)
            },
            MathematicalFunction::Multiply(multiplier, multiplicand) => {
                write!(f, "({}) * ({})", multiplier, multiplicand)
            },
            MathematicalFunction::Power(value, exponent) => {
                write!(f, "({})^({})", value, exponent)
            },
            MathematicalFunction::Sinus(value) => write!(f, "sin({})", value),
            MathematicalFunction::SinusH(value) => write!(f, "sinh({})", value),
            MathematicalFunction::Subtract(minuend, subtrahend) => {
                write!(f, "({}) - ({})", minuend, subtrahend)
            },
            MathematicalFunction::Tangens(value) => write!(f, "tan({})", value),
            MathematicalFunction::TangensH(value) => write!(f, "tanh({})", value),
        }
    }
}

impl Default for MathematicalFunction {
    fn default() -> Self {
        Self::Constant(0.0)
    }
}

impl Information for MathematicalFunction {
    fn update_value(&mut self, _time_passed: i32) {}
}

unsafe impl Send for MathematicalFunction {}

unsafe impl Sync for MathematicalFunction {}

pub mod function_mutation;
