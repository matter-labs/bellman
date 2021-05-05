use crate::pairing::ff::{Field, PrimeField, PrimeFieldRepr};
use crate::pairing::{Engine, CurveAffine, CurveProjective};
use crate::bit_vec::BitVec;

use crate::{SynthesisError};
use std::marker::PhantomData;

use crate::worker::Worker;
use crate::plonk::domains::*;
use crate::plonk::polynomials::*;

pub use crate::plonk::cs::variable::*;
use crate::plonk::better_cs::utils::*;
pub use super::lookup_tables::{self, *};

use crate::plonk::fft::cooley_tukey_ntt::*;

use super::utils::*;
pub use super::data_structures::*;
pub use super::setup::*;

pub trait SynthesisMode: Clone + Send + Sync + std::fmt::Debug {
    const PRODUCE_WITNESS: bool;
    const PRODUCE_SETUP: bool;
}
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct SynthesisModeGenerateSetup;

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct SynthesisModeProve;

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct SynthesisModeTesting;

impl SynthesisMode for SynthesisModeGenerateSetup {
    const PRODUCE_WITNESS: bool = false;
    const PRODUCE_SETUP: bool = true;
}

impl SynthesisMode for SynthesisModeProve {
    const PRODUCE_WITNESS: bool = true;
    const PRODUCE_SETUP: bool = false;
}

impl SynthesisMode for SynthesisModeTesting {
    const PRODUCE_WITNESS: bool = true;
    const PRODUCE_SETUP: bool = true;
}

pub trait Circuit<E: Engine> {
    type MainGate: MainGate<E>;
    fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError>;
    fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
        Ok(
            vec![Self::MainGate::default().into_internal()]
        )
    }
}

#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub enum Coefficient {
    PlusOne,
    MinusOne,
    Other
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct PolynomialMultiplicativeTerm(pub Coefficient, pub Vec<PolynomialInConstraint>);

impl PolynomialMultiplicativeTerm {
    fn degree(&self) -> usize {
        self.1.len()
    }
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub(crate) struct PolynomialOpeningRequest {
    pub(crate) id: PolyIdentifier,
    pub(crate) dilation: TimeDilation,
}

pub trait GateInternal<E: Engine>: Send 
    + Sync 
    + 'static 
    + std::any::Any 
    + std::fmt::Debug
{
    fn name(&self) -> &'static str;
    fn degree(&self) -> usize;
    fn can_include_public_inputs(&self) -> bool;
    fn all_queried_polynomials(&self) -> Vec<PolynomialInConstraint>;
    fn setup_polynomials(&self) -> Vec<PolyIdentifier>;
    fn variable_polynomials(&self) -> Vec<PolyIdentifier>;
    fn witness_polynomials(&self) -> Vec<PolyIdentifier> {
        vec![]
    }
    fn benefits_from_linearization(&self) -> bool;
    fn linearizes_over(&self) -> Vec<PolynomialInConstraint>;
    fn needs_opened_for_linearization(&self) -> Vec<PolynomialInConstraint>;
    fn num_quotient_terms(&self) -> usize;
    fn verify_on_row<'a>(&self, row: usize, poly_storage: &AssembledPolynomialStorage<'a, E>, last_row: bool) -> E::Fr;
    fn contribute_into_quotient<'a, 'b>(
        &self, 
        domain_size: usize,
        poly_storage: &mut AssembledPolynomialStorage<'a, E>,
        monomials_storage: & AssembledPolynomialStorageForMonomialForms<'b, E>,
        challenges: &[E::Fr],
        omegas_bitreversed: &BitReversedOmegas<E::Fr>,
        omegas_inv_bitreversed: &OmegasInvBitreversed<E::Fr>,
        worker: &Worker
    ) -> Result<Polynomial<E::Fr, Values>, SynthesisError>;
    fn contribute_into_linearization<'a>(
        &self, 
        domain_size: usize,
        at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        monomials_storage: & AssembledPolynomialStorageForMonomialForms<'a, E>,
        challenges: &[E::Fr],
        worker: &Worker
    ) -> Result<Polynomial<E::Fr, Coefficients>, SynthesisError>;
    fn contribute_into_verification_equation(
        &self, 
        domain_size: usize,
        at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        challenges: &[E::Fr],
    ) -> Result<E::Fr, SynthesisError>;
    fn put_public_inputs_into_selector_id(&self) -> Option<usize>;
    fn box_clone(&self) -> Box<dyn GateInternal<E>>;
    fn contribute_into_linearization_commitment(
        &self, 
        domain_size: usize,
        at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        commitments_storage: &std::collections::HashMap<PolyIdentifier, E::G1Affine>,
        challenges: &[E::Fr],
    ) -> Result<E::G1, SynthesisError>;
}

pub trait Gate<E: Engine>: GateInternal<E>
    + Sized
    + Clone
    + std::hash::Hash
    + std::default::Default 
{
    fn as_internal(&self) -> &dyn GateInternal<E> {
        self as &dyn GateInternal<E>
    }

    fn into_internal(self) -> Box<dyn GateInternal<E>> {
        Box::from(self) as Box<dyn GateInternal<E>>
    }
}

pub trait MainGate<E: Engine>: Gate<E> {
    const NUM_LINEAR_TERMS: usize;
    const NUM_VARIABLES: usize;
    const NUM_VARIABLES_ON_NEXT_STEP: usize;

    // fn num_linear_terms(&self) -> usize;
    // fn num_variables(&self) -> usize;
    // fn num_variables_of_next_step(&self) -> usize;

    fn range_of_multiplicative_term() -> std::ops::Range<usize>;
    fn range_of_linear_terms() -> std::ops::Range<usize>;
    fn index_for_constant_term() -> usize;
    fn range_of_next_step_linear_terms() -> std::ops::Range<usize>;
    fn format_term(instance: MainGateTerm<E>, padding: Variable) -> Result<(Vec<Variable>, Vec<E::Fr>), SynthesisError>;
    fn format_linear_term_with_duplicates(instance: MainGateTerm<E>, padding: Variable) -> Result<(Vec<Variable>, Vec<E::Fr>), SynthesisError>;
    fn dummy_vars_to_inscribe(dummy: Variable) -> Vec<Variable>;
    fn empty_coefficients() -> Vec<E::Fr>;
    fn contribute_into_quotient_for_public_inputs<'a, 'b>(
        &self, 
        domain_size: usize,
        public_inputs: &[E::Fr],
        poly_storage: &mut AssembledPolynomialStorage<'b, E>,
        monomial_storage: & AssembledPolynomialStorageForMonomialForms<'a, E>,
        challenges: &[E::Fr],
        omegas_bitreversed: &BitReversedOmegas<E::Fr>,
        omegas_inv_bitreversed: &OmegasInvBitreversed<E::Fr>,
        worker: &Worker
    ) -> Result<Polynomial<E::Fr, Values>, SynthesisError>;
    fn contribute_into_linearization_for_public_inputs<'a>(
        &self, 
        domain_size: usize,
        public_inputs: &[E::Fr],
        at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        monomials_storage: & AssembledPolynomialStorageForMonomialForms<'a, E>,
        challenges: &[E::Fr],
        worker: &Worker
    ) -> Result<Polynomial<E::Fr, Coefficients>, SynthesisError>;
    fn add_inputs_into_quotient(
        &self, 
        domain_size: usize,
        public_inputs: &[E::Fr],
        at: E::Fr,
        challenges: &[E::Fr],
    ) -> Result<E::Fr, SynthesisError>;
    // fn contribute_into_verification_equation_for_public_inputs(
    //     &self, 
    //     domain_size: usize,
    //     public_inputs: &[E::Fr],
    //     at: E::Fr,
    //     queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
    //     challenges: &[E::Fr],
    // ) -> Result<E::Fr, SynthesisError>;
    fn contribute_into_linearization_commitment_for_public_inputs(
        &self, 
        domain_size: usize,
        public_inputs: &[E::Fr],
        at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        commitments_storage: &std::collections::HashMap<PolyIdentifier, E::G1Affine>,
        challenges: &[E::Fr],
    ) -> Result<E::G1, SynthesisError>;
}

impl<E: Engine> std::hash::Hash for dyn GateInternal<E> {
    fn hash<H>(&self, state: &mut H) where H: std::hash::Hasher {
        self.type_id().hash(state);
        self.name().hash(state);
        self.degree().hash(state);
    }
}

impl<E: Engine> PartialEq for dyn GateInternal<E> {
    fn eq(&self, other: &Self) -> bool {
        self.type_id() == other.type_id() &&
        self.name() == other.name() &&
        self.degree() == other.degree()
    }
}

impl<E: Engine> Eq for dyn GateInternal<E> {}

impl<E: Engine> Clone for Box<dyn GateInternal<E>> {
    fn clone(&self) -> Self {
        self.box_clone()
    }
}
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct LinearCombinationOfTerms(pub Vec<PolynomialMultiplicativeTerm>);

impl LinearCombinationOfTerms {
    fn terms(&self) -> &[PolynomialMultiplicativeTerm] {
       &self.0[..]
    }
}

#[derive(Clone, Debug)]
pub enum ArithmeticTerm<E: Engine>{
    Product(Vec<Variable>, E::Fr),
    SingleVariable(Variable, E::Fr),
    Constant(E::Fr),
}

impl<E: Engine> ArithmeticTerm<E> {
    pub fn from_variable(var: Variable) -> Self {
        ArithmeticTerm::SingleVariable(var, E::Fr::one())
    }

    pub fn from_variable_and_coeff(var: Variable, coeff: E::Fr) -> Self {
        ArithmeticTerm::SingleVariable(var, coeff)
    }

    pub fn constant(coeff: E::Fr) -> Self {
        ArithmeticTerm::Constant(coeff)
    }

    pub fn mul_by_variable(self, other: Variable) -> Self {
        match self {
            ArithmeticTerm::Product(mut terms, coeff) => {
                terms.push(other);

                ArithmeticTerm::Product(terms, coeff)
            },
            ArithmeticTerm::SingleVariable(this, coeff) => {
                let terms = vec![this, other];

                ArithmeticTerm::Product(terms, coeff)
            },
            ArithmeticTerm::Constant(coeff) => {
                let terms = vec![other];

                ArithmeticTerm::Product(terms, coeff)
            },
        }
    }

    pub fn scale(&mut self, by: &E::Fr) {
        match self {
            ArithmeticTerm::Product(_, ref mut coeff) => {
                coeff.mul_assign(by);
            },
            ArithmeticTerm::SingleVariable(_, ref mut coeff) => {
                coeff.mul_assign(by);
            },
            ArithmeticTerm::Constant(ref mut coeff) => {
                coeff.mul_assign(by);
            },
        }
    }
}

#[derive(Clone, Debug)]
pub struct MainGateTerm<E: Engine>{
    terms: Vec<ArithmeticTerm<E>>,
    vars_scratch: std::collections::HashMap<Variable, usize>,
    num_multiplicative_terms: usize,
    num_constant_terms: usize
}

impl<E: Engine> MainGateTerm<E> {
    pub fn new() -> Self {
        Self {
            terms: Vec::with_capacity(8),
            vars_scratch: std::collections::HashMap::with_capacity(8),
            num_multiplicative_terms: 0,
            num_constant_terms: 0
        }
    }

    pub fn len_without_constant(&self) -> usize {
        self.terms.len()
    }

    pub fn add_assign(&mut self, other: ArithmeticTerm<E>) {
        match other {
            ArithmeticTerm::Product(_, _) => {
                self.num_multiplicative_terms += 1;
                self.terms.push(other);
            },
            ArithmeticTerm::SingleVariable(var, coeff) => {
                // deduplicate 
                if self.vars_scratch.get(&var).is_some() {
                    let index = *self.vars_scratch.get(&var).unwrap();
                    match &mut self.terms[index] {
                        ArithmeticTerm::SingleVariable(_, ref mut c) => {
                            c.add_assign(&coeff);
                        },
                        _ => {
                            unreachable!()
                        }
                    }
                } else {
                    // just push
                    self.vars_scratch.insert(var, self.terms.len());
                    self.terms.push(other);
                }
            },
            ArithmeticTerm::Constant(_) => {
                self.num_constant_terms += 1;
                self.terms.push(other);
            },
        }
        
        debug_assert!(self.num_constant_terms <= 1, "must duplicate constants");        
    }

    pub fn add_assign_allowing_duplicates(&mut self, other: ArithmeticTerm<E>) {
        match other {
            ArithmeticTerm::Product(_, _) => {
                self.num_multiplicative_terms += 1;
                self.terms.push(other);
            },
            ArithmeticTerm::SingleVariable(_, _) => {
                // we just push and don't even count this variable as duplicatable
                self.terms.push(other);
            },
            ArithmeticTerm::Constant(_) => {
                self.num_constant_terms += 1;
                self.terms.push(other);
            },
        }
        
        debug_assert!(self.num_constant_terms <= 1, "must duplicate constants");        
    }

    pub fn sub_assign(&mut self, mut other: ArithmeticTerm<E>) {
        match &mut other {
            ArithmeticTerm::Product(_, ref mut coeff) => {
                coeff.negate();
            },
            ArithmeticTerm::SingleVariable(_, ref mut coeff) => {
                coeff.negate();
            },
            ArithmeticTerm::Constant(ref mut coeff) => {
                coeff.negate(); 
            },
        }

        self.add_assign(other);

        debug_assert!(self.num_constant_terms <= 1, "must not duplicate constants");        
    }

    pub fn sub_assign_allowing_duplicates(&mut self, mut other: ArithmeticTerm<E>) {
        match &mut other {
            ArithmeticTerm::Product(_, ref mut coeff) => {
                coeff.negate();
            },
            ArithmeticTerm::SingleVariable(_, ref mut coeff) => {
                coeff.negate();
            },
            ArithmeticTerm::Constant(ref mut coeff) => {
                coeff.negate(); 
            },
        }

        self.add_assign_allowing_duplicates(other);

        debug_assert!(self.num_constant_terms <= 1, "must not duplicate constants");        
    }
}

#[derive(Clone, Debug, Hash, Default)]
pub struct Width4MainGateWithDNext;

impl<E: Engine> GateInternal<E> for Width4MainGateWithDNext {
    fn name(&self) -> &'static str {
        "main gate of width 4 with D_next"
    }

    fn degree(&self) -> usize {
        3
    }

    fn can_include_public_inputs(&self) -> bool {
        true
    }

    fn all_queried_polynomials(&self) -> Vec<PolynomialInConstraint> {
        let name = <Self as GateInternal<E>>::name(&self);

        vec![
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 0)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 1)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 2)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 3)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 4)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 5)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 6)),

            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)),
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)),
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(3)),

            PolynomialInConstraint::from_id_and_dilation(PolyIdentifier::VariablesPolynomial(3), 1),
        ]
    }

    fn setup_polynomials(&self) -> Vec<PolyIdentifier> {
        let name = <Self as GateInternal<E>>::name(&self);

        vec![
            PolyIdentifier::GateSetupPolynomial(name, 0),
            PolyIdentifier::GateSetupPolynomial(name, 1),
            PolyIdentifier::GateSetupPolynomial(name, 2),
            PolyIdentifier::GateSetupPolynomial(name, 3),
            PolyIdentifier::GateSetupPolynomial(name, 4),
            PolyIdentifier::GateSetupPolynomial(name, 5),
            PolyIdentifier::GateSetupPolynomial(name, 6),
        ]
    }

    fn variable_polynomials(&self) -> Vec<PolyIdentifier> {
        vec![
            PolyIdentifier::VariablesPolynomial(0),
            PolyIdentifier::VariablesPolynomial(1),
            PolyIdentifier::VariablesPolynomial(2),
            PolyIdentifier::VariablesPolynomial(3),
        ]
    }

    fn benefits_from_linearization(&self) -> bool {
        true
    }

    fn linearizes_over(&self) -> Vec<PolynomialInConstraint> {
        let name = <Self as GateInternal<E>>::name(&self);

        vec![
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 0)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 1)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 2)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 3)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 4)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 5)),
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 6)),
        ]
    }

    fn needs_opened_for_linearization(&self) -> Vec<PolynomialInConstraint> {
        vec![
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)),
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)),
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(3)),

            PolynomialInConstraint::from_id_and_dilation(PolyIdentifier::VariablesPolynomial(3), 1),
        ]
    }

    fn num_quotient_terms(&self) -> usize {
        1
    }

    fn verify_on_row(&self, row: usize, poly_storage: &AssembledPolynomialStorage<E>, last_row: bool) -> E::Fr {
        let name = <Self as GateInternal<E>>::name(&self);

        let q_a = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 0), row);
        let q_b = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 1), row);
        let q_c = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 2), row);
        let q_d = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 3), row);
        let q_m = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 4), row);
        let q_const = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 5), row);
        let q_d_next = poly_storage.get_poly_at_step(PolyIdentifier::GateSetupPolynomial(name, 6), row);

        // println!("{}*A + {}*B + {}*C + {}*D + {} + {}*A*A + {}*D_next", q_a, q_b, q_c, q_d, q_const, q_m, q_d_next);
        let a_value = poly_storage.get_poly_at_step(PolyIdentifier::VariablesPolynomial(0), row);
        let b_value = poly_storage.get_poly_at_step(PolyIdentifier::VariablesPolynomial(1), row);
        let c_value = poly_storage.get_poly_at_step(PolyIdentifier::VariablesPolynomial(2), row);
        let d_value = poly_storage.get_poly_at_step(PolyIdentifier::VariablesPolynomial(3), row);
        let d_next_value = if last_row == false {
            Some(poly_storage.get_poly_at_step(PolyIdentifier::VariablesPolynomial(3), row+1))
        } else {
            None
        }; 

        // println!("A = {}, B = {}, C = {}, D = {}, D_Next = {:?}", a_value, b_value, c_value, d_value, d_next_value);

        let mut total = E::Fr::zero();

        for (q, v) in [a_value, b_value, c_value, d_value].iter()
                    .zip([q_a, q_b, q_c, q_d].iter())
        {
            let mut tmp = *q;
            tmp.mul_assign(v);
            total.add_assign(&tmp);
        }

        total.add_assign(&q_const);

        let mut tmp = q_m;
        tmp.mul_assign(&a_value);
        tmp.mul_assign(&b_value);
        total.add_assign(&tmp);

        if last_row == false {
            let mut tmp = d_next_value.expect("must be able to get d_next");
            tmp.mul_assign(&q_d_next);

            total.add_assign(&tmp);
        } else {
            assert!(q_d_next.is_zero());
        }

        total
    }

    fn contribute_into_quotient(
        &self, 
        _domain_size: usize,
        _poly_storage: &mut AssembledPolynomialStorage<E>,
        _monomial_storage: & AssembledPolynomialStorageForMonomialForms<E>,
        _challenges: &[E::Fr],
        _omegas_bitreversed: &BitReversedOmegas<E::Fr>,
        _omegas_inv_bitreversed: &OmegasInvBitreversed<E::Fr>,
        _worker: &Worker
    ) -> Result<Polynomial<E::Fr, Values>, SynthesisError> {
        unreachable!("this type of gate can only be used as a main gate");
    }

    fn contribute_into_linearization(
        &self, 
        _domain_size: usize,
        _at: E::Fr,
        _queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        _monomials_storage: & AssembledPolynomialStorageForMonomialForms<E>,
        _challenges: &[E::Fr],
        _worker: &Worker
    ) -> Result<Polynomial<E::Fr, Coefficients>, SynthesisError> {
        unreachable!("this gate is indended to be the main gate and should use main gate functions")
    }
    fn contribute_into_verification_equation(
        &self, 
        _domain_size: usize,
        _at: E::Fr,
        _queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        _challenges: &[E::Fr],
    ) -> Result<E::Fr, SynthesisError> {
        unreachable!("this gate is indended to be the main gate and should use main gate functions")
    }

    fn put_public_inputs_into_selector_id(&self) -> Option<usize> {
        Some(5)
    }

    fn box_clone(&self) -> Box<dyn GateInternal<E>> {
        Box::from(self.clone())
    }
    fn contribute_into_linearization_commitment(
        &self, 
        domain_size: usize,
        at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        commitments_storage: &std::collections::HashMap<PolyIdentifier, E::G1Affine>,
        challenges: &[E::Fr],
    ) -> Result<E::G1, SynthesisError> {
        unreachable!("this gate is indended to be the main gate and should use main gate functions")
    }
}

impl<E: Engine> Gate<E> for Width4MainGateWithDNext {

}

impl<E: Engine> MainGate<E> for Width4MainGateWithDNext {
    const NUM_LINEAR_TERMS: usize = 4;
    const NUM_VARIABLES: usize = 4;
    const NUM_VARIABLES_ON_NEXT_STEP: usize = 1;

    fn range_of_multiplicative_term() -> std::ops::Range<usize> {
        0..2
    }
    fn range_of_linear_terms() -> std::ops::Range<usize> {
        0..4
    }

    fn index_for_constant_term() -> usize {
        5
    }

    fn range_of_next_step_linear_terms() -> std::ops::Range<usize> {
        6..7
    }

    fn format_term(mut instance: MainGateTerm<E>, padding: Variable) -> Result<(Vec<Variable>, Vec<E::Fr>), SynthesisError> {
        let mut flattened_variables = vec![padding; 4];
        let mut flattened_coefficients = vec![E::Fr::zero(); 7];
        let mut bitmap = SimpleBitmap::new();

        let allowed_linear = 4;
        let allowed_multiplications = 1;
        let allowed_constants = 1;

        let mut used_in_multiplication = [padding; 2];

        debug_assert!(instance.num_constant_terms <= allowed_constants, "must not containt more constants than allowed");  
        debug_assert!(instance.num_multiplicative_terms <= allowed_multiplications, "must not containt more multiplications than allowed"); 
        debug_assert!(instance.terms.len() <= allowed_constants + allowed_multiplications + allowed_linear, "gate can not fit that many terms"); 

        if instance.num_multiplicative_terms != 0 {
            let index = instance.terms.iter().position(
                |t| {
                    match t {
                        ArithmeticTerm::Product(_, _) => true,
                        _ => false,
                    }
                }
            ).unwrap();

            let term = instance.terms.swap_remove(index);
            match term {
                ArithmeticTerm::Product(vars, coeff) => {
                    debug_assert_eq!(vars.len(), 2, "multiplicative terms must contain two variables");

                    flattened_variables[0] = vars[0];
                    flattened_variables[1] = vars[1];
                    used_in_multiplication[0] = vars[0];
                    used_in_multiplication[1] = vars[1];
                    flattened_coefficients[4] = coeff;
                    bitmap.set(0);
                    bitmap.set(1);
                },
                _ => {
                    unreachable!("must be multiplicative term");
                }
            }
        }

        if instance.num_constant_terms != 0 {
            let index = instance.terms.iter().position(
                |t| {
                    match t {
                        ArithmeticTerm::Constant(_) => true,
                        _ => false,
                    }
                }
            ).unwrap();

            let term = instance.terms.swap_remove(index);
            match term {
                ArithmeticTerm::Constant(coeff) => {
                    flattened_coefficients[5] = coeff;
                },
                _ => {
                    unreachable!("must be constant term");
                }
            }
        }

        // only additions left
        for term in instance.terms.into_iter() {
            match term {
                ArithmeticTerm::SingleVariable(var, coeff) => {
                    let index = flattened_variables.iter().position(
                        |&t| t == var
                    );
                    if let Some(index) = index {
                        // there is some variable there already, so it must have come from multiplication
                        assert!(used_in_multiplication[0] == var || used_in_multiplication[1] == var,
                            "variable in linear term must only happen already if it was in multiplication");
                        flattened_coefficients[index] = coeff;
                    } else {
                        let idx = bitmap.get_next_unused();
                        flattened_variables[idx] = var;
                        flattened_coefficients[idx] = coeff;
                        bitmap.set(idx);
                    }
                },
                _ => {
                    unreachable!("must be additive term");
                }
            }
        }

        Ok((flattened_variables, flattened_coefficients))
    }

    fn format_linear_term_with_duplicates(mut instance: MainGateTerm<E>, padding: Variable) -> Result<(Vec<Variable>, Vec<E::Fr>), SynthesisError> {
        let mut flattened_variables = vec![padding; 4];
        let mut flattened_coefficients = vec![E::Fr::zero(); 7];
        let mut bitmap = SimpleBitmap::new();

        let allowed_linear = 4;
        let allowed_multiplications = 0;
        let allowed_constants = 1;

        debug_assert!(instance.num_constant_terms <= allowed_constants, "must not containt more constants than allowed");  
        assert!(instance.num_multiplicative_terms <= allowed_multiplications, "must not containt multiplications"); 
        debug_assert!(instance.terms.len() <= allowed_constants + allowed_multiplications + allowed_linear, "gate can not fit that many terms"); 

        if instance.num_constant_terms != 0 {
            let index = instance.terms.iter().position(
                |t| {
                    match t {
                        ArithmeticTerm::Constant(_) => true,
                        _ => false,
                    }
                }
            ).unwrap();

            let term = instance.terms.swap_remove(index);
            match term {
                ArithmeticTerm::Constant(coeff) => {
                    flattened_coefficients[5] = coeff;
                },
                _ => {
                    unreachable!("must be multiplicative term");
                }
            }
        }

        for term in instance.terms.into_iter() {
            match term {
                ArithmeticTerm::SingleVariable(var, coeff) => {
                    let idx = bitmap.get_next_unused();
                    flattened_variables[idx] = var;
                    flattened_coefficients[idx] = coeff;
                    bitmap.set(idx);
                },
                _ => {
                    unreachable!("must be multiplicative term");
                }
            }
        }

        Ok((flattened_variables, flattened_coefficients))
    }

    fn dummy_vars_to_inscribe(dummy: Variable) -> Vec<Variable> {
        vec![dummy; 4]
    }

    fn empty_coefficients() -> Vec<E::Fr> {
        vec![E::Fr::zero(); 7]
    }

    fn contribute_into_quotient_for_public_inputs<'a, 'b>(
        &self, 
        domain_size: usize,
        public_inputs: &[E::Fr],
        poly_storage: &mut AssembledPolynomialStorage<'a, E>,
        monomials_storage: & AssembledPolynomialStorageForMonomialForms<'b, E>,
        challenges: &[E::Fr],
        omegas_bitreversed: &BitReversedOmegas<E::Fr>,
        omegas_inv_bitreversed: &OmegasInvBitreversed<E::Fr>,
        worker: &Worker
    ) -> Result<Polynomial<E::Fr, Values>, SynthesisError> {
        assert!(domain_size.is_power_of_two());
        assert_eq!(challenges.len(), <Self as GateInternal<E>>::num_quotient_terms(&self));

        let lde_factor = poly_storage.lde_factor;
        assert!(lde_factor.is_power_of_two());

        assert!(poly_storage.is_bitreversed);

        let coset_factor = E::Fr::multiplicative_generator();
        // Include the public inputs
        let mut inputs_poly = Polynomial::<E::Fr, Values>::new_for_size(domain_size)?;
        for (idx, &input) in public_inputs.iter().enumerate() {
            inputs_poly.as_mut()[idx] = input;
        }
        // go into monomial form

        let mut inputs_poly = inputs_poly.ifft_using_bitreversed_ntt(&worker, omegas_inv_bitreversed, &E::Fr::one())?;

        // add constants selectors vector
        let name = <Self as GateInternal<E>>::name(&self);

        let key = PolyIdentifier::GateSetupPolynomial(name, 5);
        let constants_poly_ref = monomials_storage.get_poly(key);
        inputs_poly.add_assign(&worker, constants_poly_ref);
        drop(constants_poly_ref);

        // LDE
        let mut t_1 = inputs_poly.bitreversed_lde_using_bitreversed_ntt(
            &worker, 
            lde_factor, 
            omegas_bitreversed, 
            &coset_factor
        )?;

        for p in <Self as GateInternal<E>>::all_queried_polynomials(&self).into_iter() {
            // skip public constants poly (was used in public inputs)
            if p == PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 5)) {
                continue;
            }
            ensure_in_map_or_create(&worker, 
                p, 
                domain_size, 
                omegas_bitreversed, 
                lde_factor, 
                coset_factor, 
                monomials_storage, 
                poly_storage
            )?;
        }

        let ldes_storage = &*poly_storage;

        // Q_A * A
        let q_a_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 0)),
            ldes_storage
        );
        let a_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            ldes_storage
        );
        let mut tmp = q_a_ref.clone();
        tmp.mul_assign(&worker, a_ref);
        t_1.add_assign(&worker, &tmp);
        drop(q_a_ref);
        drop(a_ref);

        // Q_B * B
        let q_b_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 1)),
            ldes_storage
        );
        let b_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)),
            ldes_storage
        );
        tmp.reuse_allocation(q_b_ref);
        tmp.mul_assign(&worker, b_ref);
        t_1.add_assign(&worker, &tmp);
        drop(q_b_ref);
        drop(b_ref);

        // // Q_C * C
        let q_c_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 2)),
            ldes_storage
        );
        let c_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)),
            ldes_storage
        );
        tmp.reuse_allocation(q_c_ref);
        tmp.mul_assign(&worker, c_ref);
        t_1.add_assign(&worker, &tmp);
        drop(q_c_ref);
        drop(c_ref);

        // // Q_D * D
        let q_d_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 3)),
            ldes_storage
        );
        let d_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(3)),
            ldes_storage
        );
        tmp.reuse_allocation(q_d_ref);
        tmp.mul_assign(&worker, d_ref);
        t_1.add_assign(&worker, &tmp);
        drop(q_d_ref);
        drop(d_ref);

        // Q_M * A * B
        let q_m_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 4)),
            ldes_storage
        );
        let a_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            ldes_storage
        );
        let b_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)),
            ldes_storage
        );
        tmp.reuse_allocation(q_m_ref);
        tmp.mul_assign(&worker, a_ref);
        tmp.mul_assign(&worker, b_ref);
        t_1.add_assign(&worker, &tmp);
        drop(q_m_ref);
        drop(a_ref);
        drop(b_ref);

        // Q_D_next * D_next
        let q_d_next_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id(PolyIdentifier::GateSetupPolynomial(name, 6)),
            ldes_storage
        );
        let d_next_ref = get_from_map_unchecked(
            PolynomialInConstraint::from_id_and_dilation(PolyIdentifier::VariablesPolynomial(3), 1),
            ldes_storage
        );
        tmp.reuse_allocation(q_d_next_ref);
        tmp.mul_assign(&worker, d_next_ref);
        t_1.add_assign(&worker, &tmp);
        drop(q_d_next_ref);
        drop(d_next_ref);

        t_1.scale(&worker, challenges[0]);

        Ok(t_1)
    }

    fn contribute_into_linearization_for_public_inputs(
        &self, 
        _domain_size: usize,
        _public_inputs: &[E::Fr],
        _at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        monomials_storage: & AssembledPolynomialStorageForMonomialForms<E>,
        challenges: &[E::Fr],
        worker: &Worker
    ) -> Result<Polynomial<E::Fr, Coefficients>, SynthesisError> {
        // we actually do not depend on public inputs, but we use this form for consistency
        assert_eq!(challenges.len(), 1);

        let a_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let b_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let c_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let d_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(3)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let d_next_value = *queried_values.get(&PolynomialInConstraint::from_id_and_dilation(PolyIdentifier::VariablesPolynomial(3), 1))
        .ok_or(SynthesisError::AssignmentMissing)?;

        let name = <Self as GateInternal<E>>::name(&self);

        // Q_a * A
        let mut result = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 0)).clone();
        result.scale(&worker, a_value);

        // Q_b * B
        let poly_ref = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 1));
        result.add_assign_scaled(&worker, poly_ref, &b_value);

        // Q_c * C
        let poly_ref = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 2));
        result.add_assign_scaled(&worker, poly_ref, &c_value);

        // Q_d * D
        let poly_ref = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 3));
        result.add_assign_scaled(&worker, poly_ref, &d_value);

        // Q_m * A*B
        let mut tmp = a_value;
        tmp.mul_assign(&b_value);
        let poly_ref = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 4));
        result.add_assign_scaled(&worker, poly_ref, &tmp);

        // Q_const
        let poly_ref = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 5));
        result.add_assign(&worker, poly_ref);

        // Q_dNext * D_next
        let poly_ref = monomials_storage.get_poly(PolyIdentifier::GateSetupPolynomial(name, 6));
        result.add_assign_scaled(&worker, poly_ref, &d_next_value);

        result.scale(&worker, challenges[0]);

        Ok(result)
    }
    fn add_inputs_into_quotient(
        &self, 
        domain_size: usize,
        public_inputs: &[E::Fr],
        at: E::Fr,
        challenges: &[E::Fr],
    ) -> Result<E::Fr, SynthesisError> {
        if public_inputs.len() == 0 {
            return Ok(E::Fr::zero());
        }
        assert_eq!(challenges.len(), 1);
        // just evaluate L_{i}(z) * value
        let mut contribution = E::Fr::zero();
        let domain = Domain::<E::Fr>::new_for_size(domain_size as u64)?;
        for (idx, inp) in public_inputs.iter().enumerate() {
            let mut tmp = evaluate_lagrange_poly_at_point(idx, &domain, at)?;
            tmp.mul_assign(&inp);

            contribution.add_assign(&tmp);
        }

        contribution.mul_assign(&challenges[0]);

        Ok(contribution)
    }
    fn contribute_into_linearization_commitment_for_public_inputs<'a>(
        &self, 
        _domain_size: usize,
        _public_inputs: &[E::Fr],
        _at: E::Fr,
        queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
        commitments_storage: &std::collections::HashMap<PolyIdentifier, E::G1Affine>,
        challenges: &[E::Fr],
    ) -> Result<E::G1, SynthesisError> {
        // we actually do not depend on public inputs, but we use this form for consistency
        assert_eq!(challenges.len(), 1);

        let mut aggregate = E::G1::zero();

        let a_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let b_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let c_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let d_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(3)))
        .ok_or(SynthesisError::AssignmentMissing)?;
        let d_next_value = *queried_values.get(&PolynomialInConstraint::from_id_and_dilation(PolyIdentifier::VariablesPolynomial(3), 1))
        .ok_or(SynthesisError::AssignmentMissing)?;

        let name = <Self as GateInternal<E>>::name(&self);

        // Q_a * A
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 0)).ok_or(SynthesisError::AssignmentMissing)?;
        let scaled = commitment.mul(a_value.into_repr());
        aggregate.add_assign(&scaled);

        // Q_b * B
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 1)).ok_or(SynthesisError::AssignmentMissing)?;
        let scaled = commitment.mul(b_value.into_repr());
        aggregate.add_assign(&scaled);

        // Q_c * C
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 2)).ok_or(SynthesisError::AssignmentMissing)?;
        let scaled = commitment.mul(c_value.into_repr());
        aggregate.add_assign(&scaled);

        // Q_d * D
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 3)).ok_or(SynthesisError::AssignmentMissing)?;
        let scaled = commitment.mul(d_value.into_repr());
        aggregate.add_assign(&scaled);

        // Q_m * A*B
        let mut tmp = a_value;
        tmp.mul_assign(&b_value);
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 4)).ok_or(SynthesisError::AssignmentMissing)?;
        let scaled = commitment.mul(tmp.into_repr());
        aggregate.add_assign(&scaled);

        // Q_const
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 5)).ok_or(SynthesisError::AssignmentMissing)?;
        aggregate.add_assign_mixed(&commitment);

        // Q_dNext * D_next
        let commitment = commitments_storage.get(&PolyIdentifier::GateSetupPolynomial(name, 6)).ok_or(SynthesisError::AssignmentMissing)?;
        let scaled = commitment.mul(d_next_value.into_repr());
        aggregate.add_assign(&scaled);

        aggregate.mul_assign(challenges[0]);

        Ok(aggregate)
    }
}

pub fn get_from_map_unchecked<'a, 'b: 'a, E: Engine>(
    key_with_dilation: PolynomialInConstraint,
    ldes_map: &'a AssembledPolynomialStorage<'b, E>
) -> &'a Polynomial<E::Fr, Values> {

    let (key, dilation_value) = key_with_dilation.into_id_and_raw_dilation();

    let r = if dilation_value == 0 {
        match key {
            k @ PolyIdentifier::VariablesPolynomial(..) => {
                ldes_map.state_map.get(&k).expect(&format!("Must get poly {:?} from ldes storage", &k)).as_ref()
            },
            k @ PolyIdentifier::WitnessPolynomial(..) => {
                ldes_map.witness_map.get(&k).expect(&format!("Must get poly {:?} from ldes storage", &k)).as_ref()
            },           
            k @ PolyIdentifier::GateSetupPolynomial(..) => {
                ldes_map.setup_map.get(&k).expect(&format!("Must get poly {:?} from ldes storage", &k)).as_ref()
            },
            _ => {
                unreachable!();
            }
        }
    } else {
        ldes_map.scratch_space.get(&key_with_dilation).expect(&format!("Must get poly {:?} from lde storage", &key_with_dilation)).as_ref()
    };

    r
}

pub fn ensure_in_map_or_create<'a, 'b, E: Engine>(
    worker: &Worker,
    key_with_dilation: PolynomialInConstraint,
    domain_size: usize,
    omegas_bitreversed: &BitReversedOmegas<E::Fr>,
    lde_factor: usize,
    coset_factor: E::Fr,
    monomials_map: & AssembledPolynomialStorageForMonomialForms<'a, E>,
    ldes_map: &mut AssembledPolynomialStorage<'b, E>
) -> Result<(), SynthesisError> {
    assert!(ldes_map.is_bitreversed);
    assert_eq!(ldes_map.lde_factor, lde_factor);

    let (key, dilation_value) = key_with_dilation.into_id_and_raw_dilation();

    let mut contains_in_scratch_or_maps = false;

    if dilation_value == 0 {
        match key {
            k @ PolyIdentifier::VariablesPolynomial(..) => {
                if ldes_map.state_map.get(&k).is_some() {
                    contains_in_scratch_or_maps = true;
                }
            },
            k @ PolyIdentifier::WitnessPolynomial(..) => {
                if ldes_map.witness_map.get(&k).is_some() {
                    contains_in_scratch_or_maps = true;
                }
            },           
            k @ PolyIdentifier::GateSetupPolynomial(..) => {
                if ldes_map.setup_map.get(&k).is_some() {
                    contains_in_scratch_or_maps = true;
                }
            },
            _ => {
                unreachable!();
            }
        }
    } else {
        if ldes_map.scratch_space.get(&key_with_dilation).is_some() {
            contains_in_scratch_or_maps = true;
        }
    };

    if !contains_in_scratch_or_maps {
        // optimistic case: we have already calculated value without dilation
        // but now need to just rotate
        let lde_without_dilation = match key {
            k @ PolyIdentifier::VariablesPolynomial(..) => {
                ldes_map.state_map.get(&k)
            },
            k @ PolyIdentifier::WitnessPolynomial(..) => {
                ldes_map.witness_map.get(&k)
            },           
            k @ PolyIdentifier::GateSetupPolynomial(..) => {
                ldes_map.setup_map.get(&k)
            },
            _ => {
                unreachable!();
            }
        };

        let mut done = false;

        let rotated = if let Some(lde) = lde_without_dilation.as_ref() {
            let rotation_factor = dilation_value * lde_factor;
            let f = lde.as_ref().clone_shifted_assuming_bitreversed(rotation_factor, worker)?;
            drop(lde);

            Some(f)
        } else {
            None
        };

        drop(lde_without_dilation);

        if let Some(f) = rotated {
            let proxy = PolynomialProxy::from_owned(f);
            ldes_map.scratch_space.insert(key_with_dilation, proxy);

            done = true;
        };

        if !done {
            // perform LDE and push

            let monomial = match key {
                k @ PolyIdentifier::VariablesPolynomial(..) => {
                    monomials_map.state_map.get(&k).unwrap().as_ref()
                },
                k @ PolyIdentifier::WitnessPolynomial(..) => {
                    monomials_map.witness_map.get(&k).unwrap().as_ref()
                },           
                k @ PolyIdentifier::GateSetupPolynomial(..) => {
                    monomials_map.setup_map.get(&k).unwrap().as_ref()
                },
                _ => {
                    unreachable!();
                }
            };
        
            let lde = monomial.clone().bitreversed_lde_using_bitreversed_ntt(
                &worker, 
                lde_factor, 
                omegas_bitreversed, 
                &coset_factor
            )?;
        
            let final_lde = if dilation_value != 0 {
                let rotation_factor = dilation_value * lde_factor;
                let f = lde.clone_shifted_assuming_bitreversed(rotation_factor, worker)?;
                drop(lde);
        
                f
            } else {
                lde
            };

            // insert back

            let proxy = PolynomialProxy::from_owned(final_lde);

            if dilation_value == 0 {
                match key {
                    k @ PolyIdentifier::VariablesPolynomial(..) => {
                        ldes_map.state_map.insert(k, proxy);
                    },
                    k @ PolyIdentifier::WitnessPolynomial(..) => {
                        ldes_map.witness_map.insert(k, proxy);
                    },           
                    k @ PolyIdentifier::GateSetupPolynomial(..) => {
                        ldes_map.setup_map.insert(k, proxy);
                    },
                    _ => {
                        unreachable!();
                    }
                }
            } else {
                ldes_map.scratch_space.insert(key_with_dilation, proxy);
            };

            done = true;
        }

        assert!(done);
    }

    Ok(())
}

struct SimpleBitmap(u64, usize);
impl SimpleBitmap {
    fn new() -> Self {
        Self(0u64, 0)
    }

    fn get_next_unused(&mut self) -> usize {
        for i in 0..64 {
            if self.get(i) == false {
                return i;
            }
        }

        unreachable!()
    }

    fn get(&self, idx: usize) -> bool {
        1u64 << idx & self.0 > 0
    }

    fn set(&mut self, idx: usize) {
        self.0 |= 1u64 << idx;
    }
}

pub trait PlonkConstraintSystemParams<E: Engine>: Sized + Copy + Clone + Send + Sync {
    const STATE_WIDTH: usize;
    const WITNESS_WIDTH: usize;
    const HAS_WITNESS_POLYNOMIALS: bool;
    const HAS_CUSTOM_GATES: bool;
    const CAN_ACCESS_NEXT_TRACE_STEP: bool;
}

use std::sync::Arc;

pub trait ConstraintSystem<E: Engine> {
    type Params: PlonkConstraintSystemParams<E>;
    type MainGate: MainGate<E>;

    // allocate a variable
    fn alloc<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>;

    // allocate an input variable
    fn alloc_input<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>;

    fn new_single_gate_for_trace_step<G: Gate<E>>(&mut self, 
        equation: &G,
        coefficients_assignments: &[E::Fr],
        variables_assignments: &[Variable],
        witness_assignments: &[E::Fr]
    ) -> Result<(), SynthesisError> {
        self.begin_gates_batch_for_step()?;
        self.new_gate_in_batch(
            equation,
            coefficients_assignments,
            variables_assignments,
            witness_assignments
        )?;
        self.end_gates_batch_for_step()
    }

    fn get_main_gate(&self) -> &Self::MainGate;

    fn allocate_main_gate(&mut self, term: MainGateTerm<E>) -> Result<(), SynthesisError> {
        let (vars, coeffs) = Self::MainGate::format_term(term, Self::get_dummy_variable())?;

        let mg = Self::MainGate::default();

        self.new_single_gate_for_trace_step(
            &mg,
            &coeffs,
            &vars,
            &[]
        )
    }

    fn begin_gates_batch_for_step(&mut self) -> Result<(), SynthesisError>;
    fn new_gate_in_batch<G: Gate<E>>(&mut self, 
        equation: &G,
        coefficients_assignments: &[E::Fr],
        variables_assignments: &[Variable],
        witness_assignments: &[E::Fr]
    ) -> Result<(), SynthesisError>;
    fn end_gates_batch_for_step(&mut self) -> Result<(), SynthesisError>;

    fn allocate_variables_without_gate(&mut self, 
        variables_assignments: &[Variable],
        witness_assignments: &[E::Fr]
    ) -> Result<(), SynthesisError>;

    fn get_value(&self, _variable: Variable) -> Result<E::Fr, SynthesisError> { 
        Err(SynthesisError::AssignmentMissing)
    }

    fn get_dummy_variable() -> Variable;

    fn get_explicit_zero(&mut self) -> Result<Variable, SynthesisError>;
    fn get_explicit_one(&mut self) -> Result<Variable, SynthesisError>;

    fn add_table(&mut self, table: LookupTableApplication<E>) -> Result<Arc<LookupTableApplication<E>>, SynthesisError>;
    fn get_table(&self, functional_name: &str) -> Result<Arc<LookupTableApplication<E>>, SynthesisError>; 

    fn add_multitable(&mut self, table: MultiTableApplication<E>) -> Result<(), SynthesisError>;
    fn get_multitable(&self, functional_name: &str) -> Result<Arc<MultiTableApplication<E>>, SynthesisError>; 

    fn apply_single_lookup_gate(&mut self, variables: &[Variable], gate: Arc<LookupTableApplication<E>>) -> Result<(), SynthesisError>;
    fn apply_multi_lookup_gate(&mut self, variables: &[Variable], gate: Arc<MultiTableApplication<E>>) -> Result<(), SynthesisError>;

    fn get_current_step_number(&self) -> usize;
    fn get_current_aux_gate_number(&self) -> usize;
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PlonkCsWidth4WithNextStepParams;
impl<E: Engine> PlonkConstraintSystemParams<E> for PlonkCsWidth4WithNextStepParams {
    const STATE_WIDTH: usize =  4;
    const WITNESS_WIDTH: usize = 0;
    const HAS_WITNESS_POLYNOMIALS: bool = false;
    const HAS_CUSTOM_GATES: bool =  false;
    const CAN_ACCESS_NEXT_TRACE_STEP: bool = true;
}

#[derive(Clone, Copy, Debug)]
pub struct PlonkCsWidth4WithNextStepAndCustomGatesParams;

impl<E: Engine> PlonkConstraintSystemParams<E> for PlonkCsWidth4WithNextStepAndCustomGatesParams {
    const STATE_WIDTH: usize =  4;
    const WITNESS_WIDTH: usize = 0;
    const HAS_WITNESS_POLYNOMIALS: bool = false;
    const HAS_CUSTOM_GATES: bool = true;
    const CAN_ACCESS_NEXT_TRACE_STEP: bool = true;
}

use crate::plonk::polynomials::*;

#[derive(Clone)]
pub struct PolynomialStorage<E: Engine> {
    pub state_map: std::collections::HashMap<PolyIdentifier, Vec<Variable>>,
    pub witness_map: std::collections::HashMap<PolyIdentifier, Vec<E::Fr>>,
    pub setup_map: std::collections::HashMap<PolyIdentifier, Vec<E::Fr>>,
}

impl<E: Engine> PolynomialStorage<E> {
    pub fn new() -> Self {
        Self {
            state_map: std::collections::HashMap::new(),
            witness_map: std::collections::HashMap::new(),
            setup_map: std::collections::HashMap::new(),
        }
    }

    pub fn get_value(&self, poly: &PolynomialInConstraint, n: usize) -> Result<E::Fr, SynthesisError> {
        match poly {
            PolynomialInConstraint(PolyIdentifier::VariablesPolynomial(_), TimeDilation(_)) => {
                unreachable!("should not try to get value of the state polynomial, get variable first instead");
            },
            PolynomialInConstraint(PolyIdentifier::GateSetupPolynomial(gate_descr, idx), TimeDilation(dilation)) => {
                let final_index = n + dilation;
                let identifier = PolyIdentifier::GateSetupPolynomial(gate_descr, *idx);
                let value = *self.setup_map
                    .get(&identifier)
                    .ok_or(SynthesisError::AssignmentMissing)?
                    .get(final_index)
                    .ok_or(SynthesisError::AssignmentMissing)?;

                Ok(value)
            },
            PolynomialInConstraint(PolyIdentifier::WitnessPolynomial(_), TimeDilation(_)) => {
                unimplemented!()
            },
            _ => {
                unreachable!();
            }
        }
    }

    pub fn get_variable(&self, poly: &PolynomialInConstraint, n: usize) -> Result<Variable, SynthesisError> {
        match poly {
            PolynomialInConstraint(PolyIdentifier::VariablesPolynomial(idx), TimeDilation(dilation)) => {
                let final_index = n + dilation;
                let identifier = PolyIdentifier::VariablesPolynomial(*idx);
                let value = *self.state_map
                    .get(&identifier)
                    .ok_or(SynthesisError::AssignmentMissing)?
                    .get(final_index)
                    .ok_or(SynthesisError::AssignmentMissing)?;

                Ok(value)
            },
            _ => {
                unreachable!("should not try to get variable of setup or witness polynomial");
            }
        }
    }
}

#[derive(Clone)]
pub struct GateDensityStorage<E: Engine>(pub std::collections::HashMap<Box<dyn GateInternal<E>>, BitVec>);

impl<E: Engine> GateDensityStorage<E> {
    pub fn new() -> Self {
        Self(std::collections::HashMap::new())
    }
}

pub struct GateConstantCoefficientsStorage<E: Engine>(pub std::collections::HashMap<Box<dyn GateInternal<E>>, Vec<E::Fr>>);

impl<E: Engine> GateConstantCoefficientsStorage<E> {
    pub fn new() -> Self {
        Self(std::collections::HashMap::new())
    }
}

pub type TrivialAssembly<E, P, MG> = Assembly<E, P, MG, SynthesisModeTesting>;
pub type ProvingAssembly<E, P, MG> = Assembly<E, P, MG, SynthesisModeProve>;
pub type SetupAssembly<E, P, MG> = Assembly<E, P, MG, SynthesisModeGenerateSetup>;

#[derive(Clone)]
pub struct Assembly<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode> {
    pub inputs_storage: PolynomialStorage<E>,
    pub aux_storage: PolynomialStorage<E>,
    pub num_input_gates: usize,
    pub num_aux_gates: usize,
    pub max_constraint_degree: usize,
    pub main_gate: MG,
    pub input_assingments: Vec<E::Fr>,
    pub aux_assingments: Vec<E::Fr>,
    pub num_inputs: usize,
    pub num_aux: usize,
    pub trace_step_for_batch: Option<usize>,
    pub is_finalized: bool,

    pub gates: std::collections::HashSet<Box<dyn GateInternal<E>>>,
    pub all_queried_polys_in_constraints: std::collections::HashSet<PolynomialInConstraint>,
    // pub sorted_setup_polynomial_ids: Vec<PolyIdentifier>,
    pub sorted_gates: Vec<Box<dyn GateInternal<E>>>,
    pub aux_gate_density: GateDensityStorage<E>,
    pub explicit_zero_variable: Option<Variable>,
    pub explicit_one_variable: Option<Variable>,

    pub tables: Vec<Arc<LookupTableApplication<E>>>,
    pub multitables: Vec<Arc<MultiTableApplication<E>>>,
    pub table_selectors: std::collections::HashMap<String, BitVec>,
    pub multitable_selectors: std::collections::HashMap<String, BitVec>,
    pub table_ids_poly: Vec<E::Fr>,
    pub total_length_of_all_tables: usize,

    pub individual_table_entries: std::collections::HashMap<String, Vec<Vec<E::Fr>>>,
    pub individual_multitable_entries: std::collections::HashMap<String, Vec<Vec<E::Fr>>>,
    pub known_table_ids: Vec<E::Fr>,
    pub num_table_lookups: usize,
    pub num_multitable_lookups: usize,

    _marker_p: std::marker::PhantomData<P>,
    _marker_s: std::marker::PhantomData<S>,
}

impl<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode> ConstraintSystem<E> for Assembly<E, P, MG, S> {
    type Params = P;
    type MainGate = MG;

    // allocate a variable
    fn alloc<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError> 
    {

        self.num_aux += 1;
        let index = self.num_aux;
        if S::PRODUCE_WITNESS {
            let value = value()?;
            self.aux_assingments.push(value);
        }

        Ok(Variable(Index::Aux(index)))
    }

    // allocate an input variable
    fn alloc_input<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError> 
    {
        self.num_inputs += 1;
        let index = self.num_inputs;
        if S::PRODUCE_WITNESS {
            let value = value()?;
            self.input_assingments.push(value);
        }
        
        let input_var = Variable(Index::Input(index));

        let mut main_gate = MainGateTerm::<E>::new();
        main_gate.sub_assign(ArithmeticTerm::from_variable(input_var));

        let dummy = Self::get_dummy_variable();
        let (variables_assignments, coefficients_assignments) = MG::format_term(main_gate, dummy).expect("must make empty padding gate");

        let n = self.num_input_gates;
        Self::allocate_into_storage(
            &MG::default(), 
            &mut self.inputs_storage, 
            n, 
            &coefficients_assignments, 
            &variables_assignments, 
            &[]
        )?;

        self.num_input_gates += 1;

        Ok(input_var)
    }

    fn get_main_gate(&self) -> &MG {
        &self.main_gate
    }

    fn begin_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
        debug_assert!(self.trace_step_for_batch.is_none());
        self.num_aux_gates += 1;
        let n = self.num_aux_gates;
        self.trace_step_for_batch = Some(n);

        Ok(())
    }

    fn new_gate_in_batch<G: Gate<E>>(&mut self, 
        gate: &G,
        coefficients_assignments: &[E::Fr],
        variables_assignments: &[Variable],
        witness_assignments: &[E::Fr]
    ) -> Result<(), SynthesisError> {
        // check that gate is ok for config
        // debug_assert!(check_gate_is_allowed_for_params::<E, P, G>(&gate), format!("supplied params do not work with gate {:?}", gate));

        let n = self.trace_step_for_batch.unwrap();
        // make zero-enumerated index
        let n = n - 1;
        
        Self::allocate_into_storage(
            gate,
            &mut self.aux_storage,
            n,
            coefficients_assignments,
            variables_assignments,
            witness_assignments,
        )?;

        self.add_gate_into_list(gate);

        if S::PRODUCE_SETUP {
            if let Some(tracker) = self.aux_gate_density.0.get_mut(gate.as_internal() as &dyn GateInternal<E>) {
                if tracker.len() != n {
                    let padding = n - tracker.len();
                    tracker.grow(padding, false);
                }
                tracker.push(true);
                debug_assert_eq!(n+1, tracker.len());
            } else {
                self.aux_gate_density.0.insert(gate.clone().into_internal(), BitVec::new());
                let tracker = self.aux_gate_density.0.get_mut(gate.as_internal() as &dyn GateInternal<E>).unwrap();
                tracker.grow(n, false);
                tracker.push(true);
                debug_assert_eq!(n+1, tracker.len());
            }
        }

        Ok(())
    }

    fn allocate_variables_without_gate(&mut self, 
        variables_assignments: &[Variable],
        witness_assignments: &[E::Fr]
    ) -> Result<(), SynthesisError> {
        let n = self.trace_step_for_batch.expect("may only be called in a batch");
        // make zero-enumerated index
        let n = n - 1;

        let empty_coefficients = Self::MainGate::empty_coefficients();

        let gate = Self::MainGate::default();

        Self::allocate_into_storage(
            &gate,
            &mut self.aux_storage,
            n,
            &empty_coefficients,
            variables_assignments,
            witness_assignments,
        )?;

        if S::PRODUCE_SETUP {
            let apply_gate = false;

            let tracker = self.aux_gate_density.0.get_mut(gate.as_internal() as &dyn GateInternal<E>).unwrap();
            if tracker.len() != n {
                let padding = n - tracker.len();
                tracker.grow(padding, false);
            }
            tracker.push(apply_gate);
            debug_assert_eq!(n+1, tracker.len());
        }

        Ok(())
    }

    fn end_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
        debug_assert!(self.trace_step_for_batch.is_some());
        let n = self.trace_step_for_batch.take().unwrap();
        debug_assert_eq!(n, self.num_aux_gates, "invalid batch id");

        Ok(())
    }

    fn get_value(&self, var: Variable) -> Result<E::Fr, SynthesisError> {
        if !S::PRODUCE_WITNESS {
            return Err(SynthesisError::AssignmentMissing);
        }
        let value = match var {
            Variable(Index::Aux(0)) => {
                // use crate::rand::Rng;

                // let mut rng = crate::rand::thread_rng();
                // let value: E::Fr = rng.gen();

                // value
                E::Fr::zero()
                // return Err(SynthesisError::AssignmentMissing);
            }
            Variable(Index::Input(0)) => {
                return Err(SynthesisError::AssignmentMissing);
            }
            Variable(Index::Input(input)) => {
                self.input_assingments[input - 1]
            },
            Variable(Index::Aux(aux)) => {
                self.aux_assingments[aux - 1]
            }
        };

        Ok(value)
    }

    fn get_dummy_variable() -> Variable {
        Self::dummy_variable()
    }

    fn get_explicit_zero(&mut self) -> Result<Variable, SynthesisError> {
        if let Some(var) = self.explicit_zero_variable {
            return Ok(var);
        }

        let value = E::Fr::zero();
        let zero = self.alloc(|| Ok(value))?;

        let self_term = ArithmeticTerm::from_variable(zero);
        let other_term = ArithmeticTerm::constant(value);
        let mut term = MainGateTerm::new();
        term.add_assign(self_term);
        term.sub_assign(other_term);

        self.allocate_main_gate(term)?;

        self.explicit_zero_variable = Some(zero);

        Ok(zero)
    }

    fn get_explicit_one(&mut self) -> Result<Variable, SynthesisError> {
        if let Some(var) = self.explicit_one_variable {
            return Ok(var);
        }

        let value = E::Fr::one();
        let one = self.alloc(|| Ok(value))?;

        let self_term = ArithmeticTerm::from_variable(one);
        let other_term = ArithmeticTerm::constant(value);
        let mut term = MainGateTerm::new();
        term.add_assign(self_term);
        term.sub_assign(other_term);

        self.allocate_main_gate(term)?;

        self.explicit_one_variable = Some(one);

        Ok(one)
    }

    fn add_table(&mut self, table: LookupTableApplication<E>) -> Result<Arc<LookupTableApplication<E>>, SynthesisError> {
        assert!(table.applies_over().len() == 3, "only support tables of width 3");
        assert!(table.can_be_combined(), "can only add tables that are combinable");
        assert!(!self.known_table_ids.contains(&table.table_id()));
        let table_name = table.functional_name();
        let table_id = table.table_id();
        let number_of_entries = table.size();

        let shared = Arc::from(table);
        let res = shared.clone();

        self.tables.push(shared);
        self.individual_table_entries.insert(table_name.clone(), vec![]);
        self.table_selectors.insert(table_name, BitVec::new());
        self.known_table_ids.push(table_id);

        self.total_length_of_all_tables += number_of_entries;

        Ok(res)
    }

    fn get_table(&self, name: &str) -> Result<Arc<LookupTableApplication<E>>, SynthesisError> {
        for t in self.tables.iter() {
            if t.functional_name() == name {
                return Ok(Arc::clone(t));
            }
        }

        Err(SynthesisError::AssignmentMissing)
    }

    fn add_multitable(&mut self, table: MultiTableApplication<E>) -> Result<(), SynthesisError> {
        let table_name = table.functional_name();
        let mut exists = false;
        for t in self.multitables.iter() {
            if t.functional_name() == table_name {
                exists = true;
            }
        }
        assert!(exists == false);
        self.multitables.push(Arc::from(table));
        self.multitable_selectors.insert(table_name.clone(), BitVec::new());
        self.individual_table_entries.insert(table_name.clone(), vec![]);

        Ok(())
    }

    fn get_multitable(&self, functional_name: &str) -> Result<Arc<MultiTableApplication<E>>, SynthesisError> {
        for t in self.multitables.iter() {
            if t.functional_name() == functional_name {
                return Ok(Arc::clone(t));
            }
        }

        Err(SynthesisError::AssignmentMissing)
    }

    #[track_caller]
    fn apply_single_lookup_gate(&mut self, variables: &[Variable], table: Arc<LookupTableApplication<E>>) -> Result<(), SynthesisError> {
        let n = self.trace_step_for_batch.expect("may only add table constraint in a transaction");
        // make zero-enumerated index
        let n = n - 1;

        if S::PRODUCE_SETUP {
            debug_assert!(self.tables.contains(&table));
            assert!(table.can_be_combined() == true);
            assert!(table.applies_over().len() == 3);

            let table_name = table.functional_name();

            // we need to:
            // - mark that this table applies at this row
            // - add values into the list to later on make a sorted polynomial

            let tracker = self.table_selectors.get_mut(&table_name).unwrap();
            if tracker.len() != n {
                let padding = n - tracker.len();
                tracker.grow(padding, false);
            }
            tracker.push(true);
            debug_assert_eq!(n+1, tracker.len());

            // keep track of what table is applied at what row
            self.table_ids_poly.resize(n, E::Fr::zero());
            self.table_ids_poly.push(table.table_id());
        }

        if S::PRODUCE_WITNESS {
            let table_name = table.functional_name();

            // add values for lookup table sorting later
            let keys_and_values_len = table.applies_over().len();
            let mut table_entries = Vec::with_capacity(keys_and_values_len + 1);
            for v in variables.iter() {
                let value = self.get_value(*v).unwrap();
                table_entries.push(value);
            }
            table_entries.push(table.table_id());        

            let entries = self.individual_table_entries.get_mut(&table_name).unwrap();
            assert_eq!(variables.len(), table.applies_over().len());

            let valid_entries = table.is_valid_entry(&table_entries[..keys_and_values_len]);
            assert!(valid_entries);

            if !valid_entries {
                return Err(SynthesisError::Unsatisfiable);
            }

            entries.push(table_entries);
        }

        self.num_table_lookups += 1;

        Ok(())
    }
    
    #[track_caller]
    fn apply_multi_lookup_gate(&mut self, variables: &[Variable], table: Arc<MultiTableApplication<E>>) -> Result<(), SynthesisError> {
        let n = self.trace_step_for_batch.expect("may only add table constraint in a transaction");
        // make zero-enumerated index
        let n = n - 1;

        if S::PRODUCE_SETUP {
            let table_name = table.functional_name();
            let tracker = self.multitable_selectors.get_mut(&table_name).unwrap();
            if tracker.len() != n {
                let padding = n - tracker.len();
                tracker.grow(padding, false);
            }
            tracker.push(true);
            debug_assert_eq!(n+1, tracker.len());

            self.table_ids_poly.resize(n, E::Fr::zero());
            self.table_ids_poly.push(table.table_id());
        }

        if S::PRODUCE_WITNESS {
            let table_name = table.functional_name();
            let mut table_entries = Vec::with_capacity(table.applies_over().len());
            for v in variables.iter() {
                let value = self.get_value(*v).unwrap();
                table_entries.push(value);
            }
    
            let entries = self.individual_table_entries.get_mut(&table_name).unwrap();
            assert_eq!(variables.len(), table.applies_over().len());
    
            assert!(table.is_valid_entry(&table_entries));

            entries.push(table_entries);
        }

        self.num_multitable_lookups += 1;

        Ok(())
    }

    fn get_current_step_number(&self) -> usize {
        self.n()
    }

    fn get_current_aux_gate_number(&self) -> usize {
        self.num_aux_gates
    }
}

impl<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode> Assembly<E, P, MG, S> {
    fn allocate_into_storage<G: Gate<E>>(
        gate: &G,
        storage: &mut PolynomialStorage<E>, 
        n: usize,
        coefficients_assignments: &[E::Fr],
        variables_assignments: &[Variable],
        witness_assignments: &[E::Fr]
    ) -> Result<(), SynthesisError> {
        let dummy = Self::get_dummy_variable();
        let zero = E::Fr::zero();

        if S::PRODUCE_SETUP {
            let mut coeffs_it = coefficients_assignments.iter();

            for setup_poly in gate.setup_polynomials().into_iter() {
                let poly_ref = storage.setup_map.entry(setup_poly).or_insert(vec![]);
                if poly_ref.len() < n {
                    poly_ref.resize(n, E::Fr::zero());
                }
                poly_ref.push(*coeffs_it.next().unwrap_or(&zero));
            }

            debug_assert!(coeffs_it.next().is_none(), "must consume all the coefficients for gate");

        }

        let mut variable_it = variables_assignments.iter();

        for var_poly in gate.variable_polynomials().into_iter() {
            let poly_ref = storage.state_map.entry(var_poly).or_insert(vec![]);
            if poly_ref.len() < n {
                poly_ref.resize(n, dummy);
            }
            if poly_ref.len() == n {
                // we consume variable only ONCE
                let var = *variable_it.next().unwrap_or(&dummy);
                poly_ref.push(var);
            }
        }

        debug_assert!(variable_it.next().is_none(), "must consume all variables for gate");

        let mut witness_it = witness_assignments.iter();

        for key in gate.witness_polynomials().into_iter() {
            let poly_ref = storage.witness_map.entry(key).or_insert(vec![]);
            if poly_ref.len() < n {
                poly_ref.resize(n, E::Fr::zero());
            } 
            poly_ref.push(*witness_it.next().unwrap_or(&zero));
        }

        Ok(())
    }

    pub fn n(&self) -> usize {
        self.num_input_gates + self.num_aux_gates
    }

    fn add_gate_into_list<G: Gate<E>>(&mut self, gate: &G) {
        if !self.gates.contains(gate.as_internal() as &dyn GateInternal<E>) {
            self.gates.insert(gate.clone().into_internal());

            // self.add_gate_setup_polys_into_list(gate);
            for p in gate.all_queried_polynomials().into_iter() {
                self.all_queried_polys_in_constraints.insert(p);
            }
            
            self.sorted_gates.push(gate.clone().into_internal());

            let degree = gate.degree();
            if self.max_constraint_degree < degree {
                self.max_constraint_degree = degree;
            }
        }        
    }

    pub fn new() -> Self {
        let mut tmp = Self {
            inputs_storage: PolynomialStorage::new(),
            aux_storage: PolynomialStorage::new(),

            max_constraint_degree: 0,

            num_input_gates: 0,
            num_aux_gates: 0,

            num_inputs: 0,
            num_aux: 0,

            input_assingments: vec![],
            aux_assingments: vec![],

            main_gate: MG::default(),

            trace_step_for_batch: None,

            gates: std::collections::HashSet::new(),
            all_queried_polys_in_constraints: std::collections::HashSet::new(),

            aux_gate_density: GateDensityStorage::new(),

            // sorted_setup_polynomial_ids: vec![],
            sorted_gates: vec![],

            is_finalized: false,

            explicit_zero_variable: None,
            explicit_one_variable: None,

            tables: vec![],
            multitables: vec![],
            table_selectors: std::collections::HashMap::new(),
            multitable_selectors: std::collections::HashMap::new(),
            table_ids_poly: vec![],
            total_length_of_all_tables: 0,

            individual_table_entries: std::collections::HashMap::new(),
            individual_multitable_entries: std::collections::HashMap::new(),

            known_table_ids: vec![],

            num_table_lookups: 0,
            num_multitable_lookups: 0,

            _marker_p: std::marker::PhantomData,
            _marker_s: std::marker::PhantomData,
        };

        tmp.add_gate_into_list(&MG::default());

        tmp
    }

    // return variable that is not in a constraint formally, but has some value
    fn dummy_variable() -> Variable {
        Variable(Index::Aux(0))
    }

    pub fn finalize(&mut self) {
        if self.is_finalized {
            return;
        }

        let mut min_space_for_lookups = self.num_table_lookups;
        for table in self.tables.iter() {
            let table_num_rows = table.size();
            min_space_for_lookups += table_num_rows;
        }
        min_space_for_lookups += self.n();

        let new_size_candidates = [(self.n() + 1).next_power_of_two() - 1, (min_space_for_lookups + 1).next_power_of_two() - 1];

        let new_size = *new_size_candidates.iter().max().unwrap();

        let dummy = Self::get_dummy_variable();

        // let empty_gate = MainGateTerm::<E>::new();
        let empty_vars = vec![dummy; <Self as ConstraintSystem<E>>::Params::STATE_WIDTH];
        let empty_witness = vec![E::Fr::zero(); <Self as ConstraintSystem<E>>::Params::WITNESS_WIDTH];

        // let mg = MG::default();

        for _ in self.n()..new_size {

            self.begin_gates_batch_for_step().unwrap();

            self.allocate_variables_without_gate(
                &empty_vars,
                &empty_witness
            ).unwrap();

            self.end_gates_batch_for_step().unwrap();

            // self.new_single_gate_for_trace_step(
            //     &mg, 
            //     &coeffs, 
            //     &vars,
            //     &[]
            // ).expect("must add padding gate");
        }

        let new_size_for_aux = new_size - self.num_input_gates;

        if S::PRODUCE_SETUP {
            // pad gate selectors
            for (_, tracker) in self.aux_gate_density.0.iter_mut() {
                tracker.grow(new_size_for_aux, false);
            }

            // pad lookup selectors
            for (_, selector) in self.table_selectors.iter_mut() {
                selector.grow(new_size_for_aux, false);
            }

            // pad special purpose table selector poly
            self.table_ids_poly.resize(new_size_for_aux, E::Fr::zero());
        }

        assert!((self.n()+1).is_power_of_two());
        self.is_finalized = true;
    }

    fn get_storage_for_trace_step(&self, step: usize) -> &PolynomialStorage<E> {
        if step < self.num_input_gates {
            &self.inputs_storage
        } else {
            &self.aux_storage
        }
    }

    pub fn is_satisfied(&self) -> bool {
        if !S::PRODUCE_SETUP || !S::PRODUCE_WITNESS {
            // only testing mode can run this check for now
            return true;
        }
        // expect a small number of inputs

        if self.n() == 0 {
            return true;
        }

        // TODO: handle public inputs

        // for i in 0..self.num_input_gates {
        //     let gate = self.input_assingments
        // }

        // let one = E::Fr::one();
        // let mut minus_one = E::Fr::one();
        // minus_one.negate();

        let n = self.n() - 1;

        let worker = Worker::new();

        let storage = self.make_assembled_poly_storage(&worker, false).unwrap();

        for (gate_type, density) in self.aux_gate_density.0.iter() {
            for (gate_index, is_applicable) in density.iter().enumerate() {
                if is_applicable == false {
                    continue;
                }

                let trace_index = self.num_input_gates + gate_index;

                let last = trace_index == n;

                let value = gate_type.verify_on_row(trace_index, &storage, last);

                if value.is_zero() == false {
                    println!("Unsatisfied at aux gate {} (zero enumerated)", gate_index);
                    println!("Constraint value = {}", value);
                    println!("Gate {:?}", gate_type.name());
                    return false;
                }
            }
        }

        true
    }

    pub(crate) fn make_permutations(&self, worker: &Worker) -> Result<Vec<Polynomial::<E::Fr, Values>>, SynthesisError> {
        assert!(self.is_finalized);

        if !S::PRODUCE_SETUP {
            return Err(SynthesisError::AssignmentMissing);
        }

        let num_gates = self.n();
        let num_partitions = self.num_inputs + self.num_aux;
        let num_inputs = self.num_inputs;
        // in the partition number i there is a set of indexes in V = (a, b, c) such that V_j = i
        let mut partitions = vec![vec![]; num_partitions + 1];

        let mut poly_ids = vec![];
        for i in 0..P::STATE_WIDTH {
            let id = PolyIdentifier::VariablesPolynomial(i);
            poly_ids.push(id);
        }

        // gate_idx is zero-enumerated here
        for gate_idx in 0..num_gates
        {
            let storage = self.get_storage_for_trace_step(gate_idx);
            for (state_poly_index, poly_id) in poly_ids.iter().enumerate() {
                let variables_vec_ref = storage.state_map.get(&poly_id).expect("must get a variables polynomial");
                let storage_idx = if gate_idx < self.num_input_gates {
                    gate_idx
                } else {
                    gate_idx - self.num_input_gates
                };

                let v = variables_vec_ref[storage_idx];
                match v {
                    Variable(Index::Aux(0)) => {
                        // Dummy variables do not participate in the permutation
                    },
                    Variable(Index::Input(0)) => {
                        unreachable!("There must be no input with index 0");
                    },
                    Variable(Index::Input(index)) => {
                        let i = index; // inputs are [1, num_inputs]
                        partitions[i].push((state_poly_index, gate_idx+1));
                    },
                    Variable(Index::Aux(index)) => {
                        let i = index + num_inputs; // aux are [num_inputs + 1, ..]
                        partitions[i].push((state_poly_index, gate_idx+1));
                    },
                }
            }
        }

        // sanity check
        assert_eq!(partitions[0].len(), 0);

        let domain = Domain::new_for_size(num_gates as u64).expect("must have enough roots of unity to fit the circuit");

        // now we need to make root at it's cosets
        let domain_elements = materialize_domain_elements_with_natural_enumeration(
            &domain, &worker
        );

        // domain_elements.pop().unwrap();

        use crate::ff::SqrtField;
        use crate::ff::LegendreSymbol;

        let mut non_residues = vec![];
        non_residues.push(E::Fr::one());
        non_residues.extend(make_non_residues::<E::Fr>(P::STATE_WIDTH - 1));

        assert_eq!(non_residues.len(), 4);

        let mut sigmas = vec![];
        for i in 0..P::STATE_WIDTH {
            let mut sigma_i = Polynomial::from_values_unpadded(domain_elements.clone()).unwrap();
            sigma_i.scale(&worker, non_residues[i]);
            sigmas.push(sigma_i);
        }

        let mut permutations = vec![vec![]; num_partitions + 1];

        fn rotate<T: Sized>(mut vec: Vec<T>) -> Vec<T> {
            if vec.len() > 1 {
                let mut els: Vec<_> = vec.drain(0..1).collect();
                debug_assert_eq!(els.len(), 1);
                // els.reverse();
                vec.push(els.pop().unwrap());
            }

            vec
        }

        for (i, partition) in partitions.into_iter().enumerate().skip(1) {
            // copy-permutation should have a cycle around the partition

            // we do not need to permute over partitions of length 1,
            // as this variable only happends in one place
            if partition.len() == 1 {
                continue;
            }

            let permutation = rotate(partition.clone());
            permutations[i] = permutation.clone();

            for (original, new) in partition.into_iter()
                                    .zip(permutation.into_iter()) 
            {
                // (column_idx, trace_step_idx)
                let new_zero_enumerated = new.1 - 1;
                let mut new_value = domain_elements[new_zero_enumerated];

                // we have shuffled the values, so we need to determine FROM
                // which of k_i * {1, omega, ...} cosets we take a value
                // for a permutation polynomial
                new_value.mul_assign(&non_residues[new.0]);

                // check to what witness polynomial the variable belongs
                let place_into = &mut sigmas[original.0].as_mut();

                let original_zero_enumerated = original.1 - 1;
                place_into[original_zero_enumerated] = new_value;
            }
        }

        Ok(sigmas)
    }

    fn make_setup_polynomials(
        &self,
        with_finalization: bool
    ) -> Result<std::collections::HashMap<PolyIdentifier, Polynomial<E::Fr, Values>>, SynthesisError> {
        if with_finalization {
            assert!(self.is_finalized);
        }

        if !S::PRODUCE_SETUP {
            return Err(SynthesisError::AssignmentMissing);
        }

        let total_num_gates = self.n();
        let num_input_gates = self.num_input_gates;

        let mut map = std::collections::HashMap::new();

        let setup_poly_ids: Vec<_> = self.aux_storage.setup_map.keys().collect();

        for &id in setup_poly_ids.into_iter() {
            let mut assembled_poly = vec![E::Fr::zero(); total_num_gates];
            if num_input_gates != 0 {
                let input_gates_coeffs = &mut assembled_poly[..num_input_gates];
                input_gates_coeffs.copy_from_slice(&self.inputs_storage.setup_map.get(&id).unwrap()[..]);
            }

            {
                let src = &self.aux_storage.setup_map.get(&id).unwrap()[..];
                let src_len = src.len();
                let aux_gates_coeffs = &mut assembled_poly[num_input_gates..(num_input_gates+src_len)];
                aux_gates_coeffs.copy_from_slice(src);
            }

            let as_poly = Polynomial::from_values_unpadded(assembled_poly)?;

            map.insert(id, as_poly);
        }

        Ok(map)
    }

    #[track_caller]
    pub fn create_setup<C: Circuit<E>>(
        &self,
        worker: &Worker
    ) -> Result<Setup<E, C>, SynthesisError> {
        assert!(self.is_finalized);

        assert!(S::PRODUCE_SETUP);

        let claimed_gates_list = C::declare_used_gates()?;
        let known_gates_list = &self.sorted_gates;

        assert_eq!(&claimed_gates_list, known_gates_list, "trying to perform setup for a circuit that has different gates set from synthesized one: circuit claims {:?}, in synthesis we had {:?}", &claimed_gates_list, &known_gates_list);

        // check for consistency
        {
            assert!(&<Self as ConstraintSystem<E>>::MainGate::default().into_internal() == &claimed_gates_list[0]);
            assert!(&C::MainGate::default().into_internal() == &claimed_gates_list[0]);
            // dbg!(&claimed_gates_list[0]);
            // let as_any = (&claimed_gates_list[0]) as &dyn std::any::Any;
            // match as_any.downcast_ref::<<Self as ConstraintSystem<E>>::MainGate>() {
            //     Some(..) => {

            //     },
            //     None => {
            //         println!("Type mismatch: first gate among used gates must be the main gate of CS");
            //         // panic!("first gate among used gates must be the main gate of CS");
            //     }
            // }
        }

        let mut setup = Setup::<E, C>::empty();

        setup.n = self.n();
        setup.num_inputs = self.num_inputs;
        setup.state_width = <Self as ConstraintSystem<E>>::Params::STATE_WIDTH;
        setup.num_witness_polys = <Self as ConstraintSystem<E>>::Params::WITNESS_WIDTH;

        setup.non_residues = make_non_residues::<E::Fr>(setup.state_width - 1);

        let (mut setup_polys_values_map, permutation_polys) = self.perform_setup(&worker)?;
        for gate in known_gates_list.iter() {
            let setup_polys = gate.setup_polynomials();
            for id in setup_polys.into_iter() {
                let values = setup_polys_values_map.remove(&id).expect("must contain setup poly").clone_padded_to_domain()?;
                let mon = values.icoset_fft_for_generator(&worker, &E::Fr::one());

                setup.gate_setup_monomials.push(mon);
            }
        }

        for perm in permutation_polys.into_iter() {
            let mon = perm.icoset_fft_for_generator(&worker, &E::Fr::one());

            setup.permutation_monomials.push(mon);
        }

        let gate_selector_values = self.output_gate_selectors(&worker)?;
        
        if known_gates_list.len() > 1 {
            assert_eq!(gate_selector_values.len(), known_gates_list.len(), "numbers of selectors and known gates mismatch");
        }

        for values in gate_selector_values.into_iter() {
            let poly = Polynomial::from_values(values)?;
            let mon = poly.icoset_fft_for_generator(&worker, &E::Fr::one());

            setup.gate_selectors_monomials.push(mon);
        }

        if self.tables.len() > 0 && self.num_table_lookups > 0 {
            // we have lookup tables, so add them to setup

            let num_lookups = self.num_table_lookups;
            setup.total_lookup_entries_length = num_lookups;

            let table_tails = self.calculate_t_polynomial_values_for_single_application_tables()?;
            assert_eq!(table_tails.len(), 4);

            let tails_len = table_tails[0].len();

            // total number of gates, Input + Aux
            let size = self.n();

            let copy_start = size - tails_len;

            for tail in table_tails.into_iter() {
                let mut values = vec![E::Fr::zero(); size];
                values[copy_start..].copy_from_slice(&tail[..]);

                let poly = Polynomial::from_values(values)?;
                let mon = poly.icoset_fft_for_generator(&worker, &E::Fr::one());

                setup.lookup_tables_monomials.push(mon);
            }

            let selector_for_lookup_values = self.calculate_lookup_selector_values()?;
            let poly = Polynomial::from_values(selector_for_lookup_values)?;
            let mon = poly.icoset_fft_for_generator(&worker, &E::Fr::one());
            setup.lookup_selector_monomial = Some(mon);

            let table_type_values = self.calculate_table_type_values()?;
            let poly = Polynomial::from_values(table_type_values)?;
            let mon = poly.icoset_fft_for_generator(&worker, &E::Fr::one());
            setup.lookup_table_type_monomial = Some(mon);
        }

        Ok(setup)
    }

    pub fn perform_setup(
        &self, 
        worker: &Worker
    ) -> Result<
    (std::collections::HashMap<PolyIdentifier, Polynomial<E::Fr, Values>>, Vec<Polynomial<E::Fr, Values>>), 
    SynthesisError
    > {
        let map = self.make_setup_polynomials(true)?;
        let permutation_polys = self.make_permutations(&worker)?;

        Ok((map, permutation_polys))
    }

    pub fn output_gate_selectors(&self, worker: &Worker) -> Result<Vec<Vec<E::Fr>>, SynthesisError> {
        // if self.sorted_gates.len() == 1 {
        //     return Ok(vec![]);
        // }

        let num_gate_selectors = self.sorted_gates.len();

        let one = E::Fr::one();
        let empty_poly_values = vec![E::Fr::zero(); self.n()];
        let mut poly_values = vec![empty_poly_values.clone(); num_gate_selectors];
        let num_input_gates = self.num_input_gates;

        // first gate in sorted in a main gate and applies on public inputs
        for p in poly_values[0][..num_input_gates].iter_mut() {
            *p = one;
        }

        worker.scope(poly_values.len(), |scope, chunk| {
            for (i, lh) in poly_values.chunks_mut(chunk)
                            .enumerate() {
                scope.spawn(move |_| {
                    // we take `values_per_leaf` values from each of the polynomial
                    // and push them into the conbinations
                    let base_idx = i*chunk;
                    for (j, lh) in lh.iter_mut().enumerate() {
                        let idx = base_idx + j;
                        let id = &self.sorted_gates[idx];
                        let density = self.aux_gate_density.0.get(id).unwrap();
                        let poly_mut_slice: &mut [E::Fr] = &mut lh[num_input_gates..];
                        for (i, d) in density.iter().enumerate() {
                            if d {
                                poly_mut_slice[i] = one;
                            }
                        }
                    }
                });
            }
        });

        Ok(poly_values)
    }

    pub fn calculate_t_polynomial_values_for_single_application_tables(&self) -> 
        Result<Vec<Vec<E::Fr>>, SynthesisError> {

        if !S::PRODUCE_SETUP {
            return Err(SynthesisError::AssignmentMissing);
        }
        
        if self.tables.len() == 0 {
            return Ok(vec![])
        }
        
        // we should pass over every table and append it

        let mut width = 0;
        for table in self.tables.iter() {
            if width == 0 {
                width = table.width();
            } else {
                assert_eq!(width, table.width());
            }
        }

        assert_eq!(width, 3, "only support tables that span over 3 polynomials for now");

        let mut column_contributions = vec![vec![]; width + 1];

        for table in self.tables.iter() {
            let entries = table.get_table_values_for_polys();
            // these are individual column vectors, so just copy
            for (idx, e) in entries.into_iter().enumerate() {
                column_contributions[idx].extend(e);
            }

            let table_id = table.table_id();
            let pad_to_len = column_contributions[0].len();

            column_contributions.last_mut().unwrap().resize(pad_to_len, table_id);
        }


        Ok(column_contributions)
    }

    pub fn calculate_interleaved_t_polys(&self) -> 
        Result<Vec<Vec<E::Fr>>, SynthesisError> {
        assert!(self.is_finalized);
        
        if self.tables.len() == 0 {
            return Ok(vec![])
        }
        
        // we should pass over every table and append it

        let mut width = 0;
        for table in self.tables.iter() {
            if width == 0 {
                width = table.width();
            } else {
                assert_eq!(width, table.width());
            }
        }

        assert_eq!(width, 3, "only support tables that span over 3 polynomials for now");

        assert!(self.is_finalized);
        if self.tables.len() == 0 {
            return Ok(vec![]);
        }

        // total number of gates, Input + Aux
        let size = self.n();

        let aux_gates_start = self.num_input_gates;

        let mut contributions = vec![vec![E::Fr::zero(); size]; 4];

        // make it shifted for ease of rotations later one
        let mut place_into_idx = aux_gates_start + 1;

        let lookup_selector = self.calculate_lookup_selector_values()?;

        for single_application in self.tables.iter() {
            let entries = single_application.get_table_values_for_polys();
            assert_eq!(entries.len(), 3);
            let table_id = single_application.table_id();
            let num_entries = single_application.size();

            assert_eq!(entries[0].len(), num_entries);

            for entry_idx in 0..num_entries {
                'inner: loop {
                    if lookup_selector[place_into_idx].is_zero() {
                        // we can place a table row into the poly
                        for column in 0..3 {
                            contributions[column][place_into_idx] = entries[column][entry_idx];
                        }
                        contributions[3][place_into_idx] = table_id;
                        place_into_idx += 1;

                        break 'inner;
                    } else {
                        // go for a next one 
                        place_into_idx += 1;
                    }
                }
            }
        }

        Ok(contributions)
    }

    pub fn calculate_s_poly_contributions_from_witness(&self) ->
        Result<Vec<Vec<E::Fr>>, SynthesisError> 
    {
        if self.tables.len() == 0 {
            return Ok(vec![]);
        }
        // we first form a set of all occured witness values,
        // then include table entries to the set
        // and then sort this set

        let mut kv_set_entries = vec![];
        let mut contributions_per_column = vec![vec![]; 4];
        for single_application in self.tables.iter() {
            // copy all queries from witness
            let table_name = single_application.functional_name();
            for kv_values in self.individual_table_entries.get(&table_name).unwrap() {
                let entry = KeyValueSet::<E>::from_slice(&kv_values[..3]);
                kv_set_entries.push(entry);
            }

            // copy table elements themselves

            let entries = single_application.get_table_values_for_polys();
            // those are full values of polynomials, so we have to virtually transpose

            let size = entries[0].len();
            for i in 0..size {
                let entry = KeyValueSet::new([entries[0][i], entries[1][i], entries[2][i]]);
                kv_set_entries.push(entry)
            }

            kv_set_entries.sort();

            // now copy backward with addition of the table id

            for kv in kv_set_entries.iter() {
                for i in 0..3 {
                    contributions_per_column[i].push(kv.inner[i]);
                }
            }

            let table_id = single_application.table_id();
            let pad_to_len = contributions_per_column[0].len();
            contributions_per_column.last_mut().unwrap().resize(pad_to_len, table_id);

            kv_set_entries.truncate(0);
        }

        Ok(contributions_per_column)
    }

    pub fn calculate_table_type_values(
        &self
    ) -> 
        Result<Vec<E::Fr>, SynthesisError>
    {
        assert!(self.is_finalized);

        if !S::PRODUCE_SETUP {
            return Err(SynthesisError::AssignmentMissing);
        }

        if self.tables.len() == 0 {
            return Ok(vec![]);
        }

        let table_ids_vector_on_aux_gates = &self.table_ids_poly;

        let num_aux_gates = self.num_aux_gates;

        // total number of gates, Input + Aux
        let size = self.n();

        let aux_gates_start = self.num_input_gates;
        let aux_gates_end = aux_gates_start + num_aux_gates;

        let mut values = vec![E::Fr::zero(); size];
        assert_eq!(num_aux_gates, table_ids_vector_on_aux_gates.len());

        values[aux_gates_start..aux_gates_end].copy_from_slice(table_ids_vector_on_aux_gates);

        Ok(values)
    }

    pub fn calculate_lookup_selector_values(
        &self
    ) -> Result<Vec<E::Fr>, SynthesisError> {
        assert!(self.is_finalized);

        if !S::PRODUCE_SETUP {
            return Err(SynthesisError::AssignmentMissing);
        }

        if self.tables.len() == 0 {
            return Ok(vec![]);
        }

        // total number of gates, Input + Aux
        let size = self.n();

        let aux_gates_start = self.num_input_gates;

        let num_aux_gates = self.num_aux_gates;
        // input + aux gates without t-polys

        let mut lookup_selector_values = vec![E::Fr::zero(); size];

        for single_application in self.tables.iter() {
            let table_name = single_application.functional_name();
            let selector_bitvec = self.table_selectors.get(&table_name).unwrap();

            for aux_gate_idx in 0..num_aux_gates {
                if selector_bitvec[aux_gate_idx] {
                    let global_gate_idx = aux_gate_idx + aux_gates_start;
                    // place 1 into selector
                    lookup_selector_values[global_gate_idx] = E::Fr::one();
                }
            }
        }

        Ok(lookup_selector_values)
    }

    pub fn calculate_masked_lookup_entries(
        &self,
        storage: &AssembledPolynomialStorage<E>
    ) -> 
        Result<Vec<Vec<E::Fr>>, SynthesisError>
    {
        assert!(self.is_finalized);
        if self.tables.len() == 0 {
            return Ok(vec![]);
        }

        // total number of gates, Input + Aux
        let size = self.n();

        let aux_gates_start = self.num_input_gates;

        let num_aux_gates = self.num_aux_gates;
        // input + aux gates without t-polys

        let mut contributions_per_column = vec![vec![E::Fr::zero(); size]; 3];
        for single_application in self.tables.iter() {
            let table_name = single_application.functional_name();
            let keys_and_values = single_application.applies_over();
            let selector_bitvec = self.table_selectors.get(&table_name).unwrap();

            for aux_gate_idx in 0..num_aux_gates {
                if selector_bitvec[aux_gate_idx] {
                    let global_gate_idx = aux_gate_idx + aux_gates_start;

                    // place value into f poly
                    for (idx, &poly_id) in keys_and_values.iter().enumerate() {
                        let value = storage.get_poly_at_step(poly_id, global_gate_idx);
                        contributions_per_column[idx][global_gate_idx] = value;
                    }
                }
            }
        }

        Ok(contributions_per_column)
    }

    pub fn calculate_masked_lookup_entries_using_selector<'a>(
        &self,
        storage: &AssembledPolynomialStorage<E>,
        selector: &PolynomialProxy<'a, E::Fr, Values>
    ) -> 
        Result<Vec<Vec<E::Fr>>, SynthesisError>
    {
        assert!(self.is_finalized);
        if self.tables.len() == 0 {
            return Ok(vec![]);
        }

        // total number of gates, Input + Aux
        let size = self.n();

        let aux_gates_start = self.num_input_gates;

        let num_aux_gates = self.num_aux_gates;
        // input + aux gates without t-polys

        let selector_ref = selector.as_ref().as_ref();

        let one = E::Fr::one();

        let mut contributions_per_column = vec![vec![E::Fr::zero(); size]; 3];
        for single_application in self.tables.iter() {
            let keys_and_values = single_application.applies_over();

            for aux_gate_idx in 0..num_aux_gates {
                let global_gate_idx = aux_gate_idx + aux_gates_start;

                if selector_ref[global_gate_idx] == one {
                    // place value into f poly
                    for (idx, &poly_id) in keys_and_values.iter().enumerate() {
                        let value = storage.get_poly_at_step(poly_id, global_gate_idx);
                        contributions_per_column[idx][global_gate_idx] = value;
                    }
                }
            }
        }

        Ok(contributions_per_column)
    }

    fn sort_by_t(
        witness_entries: &Vec<Vec<E::Fr>>,
        table_entries: &Vec<Vec<E::Fr>>,
    ) -> Result< Vec<Vec<E::Fr>>, SynthesisError> {
        assert_eq!(witness_entries.len(), table_entries.len());

        if witness_entries.len() == 0 {
            return Ok(vec![]);
        }

        // make s = f sorted by t (elements in s appear in the same order as elements in t)

        let entries_len = table_entries[0].len();
        let witnesses_len = witness_entries[0].len();

        let mut index_lookups_for_sorting = vec![std::collections::HashMap::with_capacity(entries_len); witness_entries.len()];

        for (idx, table) in table_entries.iter().enumerate() {
            for (entry_index, &entry_value) in table.iter().enumerate() {
                // make a reverse lookup field element -> index
                index_lookups_for_sorting[idx].insert(entry_value, entry_index);
            }
        }

        let mut column_contributions = vec![];

        for (idx, witness_column) in witness_entries.iter().enumerate() {
            let mut indexes = vec![usize::max_value(); witnesses_len];
            for (witness_index, witness_value) in witness_column.iter().enumerate() {
                let reverse_lookup_index = index_lookups_for_sorting[idx].get(witness_value).unwrap();
                indexes[witness_index] = *reverse_lookup_index;
            }

            indexes.sort();

            println!("sorted_index = {:?}", indexes);

            let mut s_for_column = Vec::with_capacity(witnesses_len);

            for sorted_index in indexes.into_iter() {
                let table = &table_entries[idx];
                s_for_column.push(table[sorted_index]);
            }

            column_contributions.push(s_for_column);
        }


        Ok(column_contributions)
    }

    pub fn make_state_and_witness_polynomials(
        &self,
        worker: &Worker,
        with_finalization: bool
    ) -> Result<(Vec<Vec<E::Fr>>, Vec<Vec<E::Fr>>), SynthesisError>
    {
        if with_finalization {
            assert!(self.is_finalized);
        }

        if !S::PRODUCE_WITNESS {
            return Err(SynthesisError::AssignmentMissing);
        }

        let mut full_assignments = if with_finalization {
            vec![Vec::with_capacity((self.n()+1).next_power_of_two()); P::STATE_WIDTH]
        } else {
            vec![Vec::with_capacity(self.n()+1); P::STATE_WIDTH]
        };

        let pad_to = if with_finalization {
            (self.n()+1).next_power_of_two()
        } else {
            self.n()+1
        };

        let num_input_gates = self.num_input_gates;
        let num_aux_gates = self.num_aux_gates;

        full_assignments[0].extend_from_slice(&self.input_assingments);
        assert!(full_assignments[0].len() == num_input_gates);
        for i in 1..P::STATE_WIDTH {
            full_assignments[i].resize(num_input_gates, E::Fr::zero());
        }

        let dummy = Self::get_dummy_variable();

        worker.scope(full_assignments.len(), |scope, chunk| {
            for (i, lh) in full_assignments.chunks_mut(chunk)
                            .enumerate() {
                scope.spawn(move |_| {
                    // we take `values_per_leaf` values from each of the polynomial
                    // and push them into the conbinations
                    let base_idx = i*chunk;
                    for (j, lh) in lh.iter_mut().enumerate() {
                        let idx = base_idx + j;
                        let id = PolyIdentifier::VariablesPolynomial(idx);
                        let poly_ref = self.aux_storage.state_map.get(&id).unwrap();
                        for i in 0..num_aux_gates {
                            let var = poly_ref.get(i).unwrap_or(&dummy);
                            let value = self.get_value(*var).unwrap();
                            lh.push(value);
                        }
                    }
                });
            }
        });

        for p in full_assignments.iter_mut() {
            p.resize(pad_to - 1, E::Fr::zero());
        }

        for a in full_assignments.iter() {
            assert_eq!(a.len(), pad_to - 1);
        }

        Ok((full_assignments, vec![]))
    }

    pub fn make_assembled_poly_storage<'a>(
        &self, 
        worker: &Worker, 
        with_finalization: bool
    ) -> Result<AssembledPolynomialStorage<'a, E>, SynthesisError> {
        if with_finalization {
            assert!(self.is_finalized);
        }
        
        let (state_polys, witness_polys) = self.make_state_and_witness_polynomials(&worker, with_finalization)?;

        let mut state_polys_map = std::collections::HashMap::new();
        for (idx, poly) in state_polys.into_iter().enumerate() {
            let key = PolyIdentifier::VariablesPolynomial(idx);
            let p = Polynomial::from_values_unpadded(poly)?;
            let p = PolynomialProxy::from_owned(p);
            state_polys_map.insert(key, p);
        }

        let mut witness_polys_map = std::collections::HashMap::new();
        for (idx, poly) in witness_polys.into_iter().enumerate() {
            let key = PolyIdentifier::WitnessPolynomial(idx);
            let p = Polynomial::from_values_unpadded(poly)?;
            let p = PolynomialProxy::from_owned(p);
            witness_polys_map.insert(key, p);
        }

        let mut setup_map = std::collections::HashMap::new();
        let mut gate_selectors_map = std::collections::HashMap::new();

        if S::PRODUCE_SETUP {
            let setup_polys_map = self.make_setup_polynomials(with_finalization)?;
            let gate_selectors = self.output_gate_selectors(&worker)?;

            for (gate, poly) in self.sorted_gates.iter().zip(gate_selectors.into_iter()) {
                // let key = gate.clone();
                let key = PolyIdentifier::GateSelector(gate.name());
                let p = Polynomial::from_values_unpadded(poly)?;
                let p = PolynomialProxy::from_owned(p);
                gate_selectors_map.insert(key, p);
            }

            for (key, p) in setup_polys_map.into_iter() {
                let p = PolynomialProxy::from_owned(p);
                setup_map.insert(key, p);
            }
        }

        let assembled = AssembledPolynomialStorage::<E> {
            state_map: state_polys_map,
            witness_map: witness_polys_map,
            setup_map: setup_map,
            scratch_space: std::collections::HashMap::new(),
            gate_selectors: gate_selectors_map,
            named_polys: std::collections::HashMap::new(),
            is_bitreversed: false,
            lde_factor: 1
        };

        Ok(assembled)
    }

    pub fn create_monomial_storage<'a, 'b>(
        worker: &Worker,
        omegas_inv: &OmegasInvBitreversed<E::Fr>,
        value_form_storage: &'a AssembledPolynomialStorage<E>,
        include_setup: bool,
    ) -> Result<AssembledPolynomialStorageForMonomialForms<'b, E>, SynthesisError> {
        assert_eq!(value_form_storage.lde_factor, 1);
        assert!(value_form_storage.is_bitreversed == false);

        let mut monomial_storage = AssembledPolynomialStorageForMonomialForms::<E>::new();

        for (&k, v) in value_form_storage.state_map.iter() {
            let mon_form = v.as_ref().clone_padded_to_domain()?.ifft_using_bitreversed_ntt(
                &worker, 
                omegas_inv, 
                &E::Fr::one()
            )?;
            let mon_form = PolynomialProxy::from_owned(mon_form);
            monomial_storage.state_map.insert(k, mon_form);
        }

        for (&k, v) in value_form_storage.witness_map.iter() {
            let mon_form = v.as_ref().clone_padded_to_domain()?.ifft_using_bitreversed_ntt(
                &worker, 
                omegas_inv, 
                &E::Fr::one()
            )?;
            let mon_form = PolynomialProxy::from_owned(mon_form);
            monomial_storage.witness_map.insert(k, mon_form);
        }

        if include_setup {
            for (&k, v) in value_form_storage.gate_selectors.iter() {
                let mon_form = v.as_ref().clone_padded_to_domain()?.ifft_using_bitreversed_ntt(
                    &worker, 
                    omegas_inv, 
                    &E::Fr::one()
                )?;
                let mon_form = PolynomialProxy::from_owned(mon_form);
                monomial_storage.gate_selectors.insert(k, mon_form);
            }

            for (&k, v) in value_form_storage.setup_map.iter() {
                let mon_form = v.as_ref().clone_padded_to_domain()?.ifft_using_bitreversed_ntt(
                    &worker, 
                    omegas_inv, 
                    &E::Fr::one()
                )?;
                let mon_form = PolynomialProxy::from_owned(mon_form);
                monomial_storage.setup_map.insert(k, mon_form);
            }
        }

        Ok(monomial_storage)
    }

    pub fn prover_stub(self, worker: &Worker) -> Result<(), SynthesisError> {
        use crate::pairing::CurveAffine;

        assert!(self.is_finalized);

        let num_state_polys = <Self as ConstraintSystem<E>>::Params::STATE_WIDTH;
        let num_witness_polys = <Self as ConstraintSystem<E>>::Params::WITNESS_WIDTH;

        let mut values_storage = self.make_assembled_poly_storage(worker, true)?;
        let permutation_polys = self.make_permutations(&worker)?;
        assert_eq!(permutation_polys.len(), num_state_polys);

        for (idx, poly) in permutation_polys.into_iter().enumerate() {
            let key = PolyIdentifier::PermutationPolynomial(idx);
            let poly = PolynomialProxy::from_owned(poly);
            values_storage.setup_map.insert(key, poly);
        }

        let required_domain_size = self.n() + 1;
        assert!(required_domain_size.is_power_of_two());

        let omegas_bitreversed = BitReversedOmegas::<E::Fr>::new_for_domain_size(required_domain_size);
        let omegas_inv_bitreversed = <OmegasInvBitreversed::<E::Fr> as CTPrecomputations::<E::Fr>>::new_for_domain_size(required_domain_size);

        let mut ldes_storage = AssembledPolynomialStorage::<E>::new(
            true, 
            self.max_constraint_degree.next_power_of_two()
        );
        let monomials_storage = Self::create_monomial_storage(
            &worker, 
            &omegas_inv_bitreversed, 
            &values_storage, 
            true
        )?;

        // step 1 - commit state and witness, enumerated. Also commit sorted polynomials for table arguments
        let mut commitments = vec![];

        for i in 0..num_state_polys {
            // TODO
            let id = PolyIdentifier::VariablesPolynomial(i);
            let commitment = E::G1Affine::zero();
            commitments.push(commitment);
        }

        for i in 0..num_witness_polys {
            // TODO
            let id = PolyIdentifier::WitnessPolynomial(i);
            let commitment = E::G1Affine::zero();
            commitments.push(commitment);
        }

        // step 1.5 - if there are lookup tables then draw random "eta" to linearlize over tables
        let mut lookup_data: Option<lookup_tables::LookupDataHolder<E>> = if self.tables.len() > 0 {
            // TODO: Some of those values are actually part of the setup, but deal with it later

            let eta = E::Fr::from_str("987").unwrap();

            // these are selected rows from witness (where lookup applies)
            // TODO: split and output selector without calculating witness

            let selector_for_lookup_values = self.calculate_lookup_selector_values()?;
            let table_type_values = self.calculate_table_type_values()?;

            let mut table_contributions_values = self.calculate_masked_lookup_entries(&values_storage)?;

            assert_eq!(table_contributions_values.len(), 3);

            let witness_len = table_contributions_values[0].len();
            assert_eq!(witness_len, required_domain_size - 1);

            let mut f_poly_values_aggregated = table_contributions_values.drain(0..1).collect::<Vec<_>>().pop().unwrap();

            let mut current = eta;
            for t in table_contributions_values.into_iter() {
                let op = BinopAddAssignScaled::new(current);
                binop_over_slices(&worker, &op, &mut f_poly_values_aggregated, &t);
                
                current.mul_assign(&eta);
            }

            // add table type marker
            let op = BinopAddAssignScaled::new(current);
            binop_over_slices(&worker, &op, &mut f_poly_values_aggregated, &table_type_values);

            let table_type_poly_monomial = {
                let mon = Polynomial::from_values(table_type_values)?;
                let mon = mon.ifft_using_bitreversed_ntt(
                    &worker,
                    &omegas_inv_bitreversed,
                    &E::Fr::one()
                )?;

                mon
            };

            // these are unsorted rows of lookup tables
            let mut t_poly_ends = self.calculate_t_polynomial_values_for_single_application_tables()?;

            assert_eq!(t_poly_ends.len(), 4);

            let mut t_poly_values_aggregated = t_poly_ends.drain(0..1).collect::<Vec<_>>().pop().unwrap();
            let mut current = eta;
            for t in t_poly_ends.into_iter() {
                let op = BinopAddAssignScaled::new(current);
                binop_over_slices(&worker, &op, &mut t_poly_values_aggregated, &t);
                
                current.mul_assign(&eta);
            }

            let copy_start = witness_len - t_poly_values_aggregated.len();
            let mut full_t_poly_values = vec![E::Fr::zero(); witness_len];
            let mut full_t_poly_values_shifted = full_t_poly_values.clone();

            full_t_poly_values[copy_start..].copy_from_slice(&t_poly_values_aggregated);
            full_t_poly_values_shifted[(copy_start - 1)..(witness_len-1)].copy_from_slice(&t_poly_values_aggregated);

            assert!(full_t_poly_values[0].is_zero());

            dbg!(&full_t_poly_values);

            // // these are unsorted rows of lookup tables
            // let mut interleaved_t_polys = self.calculate_interleaved_t_polys()?;

            // assert_eq!(interleaved_t_polys.len(), 4);
            // assert_eq!(interleaved_t_polys[0].len(), witness_len);

            // let mut t_poly_values_aggregated = interleaved_t_polys.drain(0..1).collect::<Vec<_>>().pop().unwrap();
            // let mut current = eta;
            // for t in interleaved_t_polys.into_iter() {
            //     let op = BinopAddAssignScaled::new(current);
            //     binop_over_slices(&worker, &op, &mut t_poly_values_aggregated, &t);
                
            //     current.mul_assign(&eta);
            // }

            // let mut full_t_poly_values = vec![E::Fr::zero(); witness_len];
            // let mut full_t_poly_values_shifted = full_t_poly_values.clone();

            // full_t_poly_values[..].copy_from_slice(&t_poly_values_aggregated);
            // full_t_poly_values_shifted[..(witness_len-1)].copy_from_slice(&t_poly_values_aggregated[1..]);

            // assert!(full_t_poly_values[0].is_zero());

            let mut s_poly_ends = self.calculate_s_poly_contributions_from_witness()?;
            assert_eq!(s_poly_ends.len(), 4);

            let mut s_poly_values_aggregated = s_poly_ends.drain(0..1).collect::<Vec<_>>().pop().unwrap();
            let mut current = eta;
            for t in s_poly_ends.into_iter() {
                let op = BinopAddAssignScaled::new(current);
                binop_over_slices(&worker, &op, &mut s_poly_values_aggregated, &t);
                
                current.mul_assign(&eta);
            }

            let sorted_copy_start = witness_len - s_poly_values_aggregated.len();

            let mut full_s_poly_values = vec![E::Fr::zero(); witness_len];
            let mut full_s_poly_values_shifted = full_s_poly_values.clone();

            full_s_poly_values[sorted_copy_start..].copy_from_slice(&s_poly_values_aggregated);
            full_s_poly_values_shifted[(sorted_copy_start - 1)..(witness_len-1)].copy_from_slice(&s_poly_values_aggregated);

            assert!(full_s_poly_values[0].is_zero());

            let selector_poly = Polynomial::<E::Fr, Values>::from_values(selector_for_lookup_values)?.ifft_using_bitreversed_ntt(
                &worker,
                &omegas_inv_bitreversed,
                &E::Fr::one()
            )?;

            let data = lookup_tables::LookupDataHolder::<E> {
                eta,
                f_poly_unpadded_values: Some(Polynomial::from_values_unpadded(f_poly_values_aggregated)?),
                t_poly_unpadded_values: Some(Polynomial::from_values_unpadded(full_t_poly_values)?),
                t_shifted_unpadded_values: Some(Polynomial::from_values_unpadded(full_t_poly_values_shifted)?),
                s_poly_unpadded_values: Some(Polynomial::from_values_unpadded(full_s_poly_values)?),
                s_shifted_unpadded_values: Some(Polynomial::from_values_unpadded(full_s_poly_values_shifted)?),
                // f_poly_monomial: None,
                t_poly_monomial: None,
                s_poly_monomial: None,
                selector_poly_monomial: Some(selector_poly),
                table_type_poly_monomial: Some(table_type_poly_monomial),
            };

            Some(data)
        } else {
            None
        };

        if self.multitables.len() > 0 {
            unimplemented!("do not support multitables yet")
        }

        // step 2 - grand product arguments

        let beta_for_copy_permutation = E::Fr::from_str("123").unwrap();
        let gamma_for_copy_permutation = E::Fr::from_str("456").unwrap();

        // copy permutation grand product argument

        let mut grand_products_protos_with_gamma = vec![];

        for i in 0..num_state_polys {
            let id = PolyIdentifier::VariablesPolynomial(i);

            let mut p = values_storage.state_map.get(&id).unwrap().as_ref().clone();
            p.add_constant(&worker, &gamma_for_copy_permutation);

            grand_products_protos_with_gamma.push(p);
        }
        
        let required_domain_size = required_domain_size;

        let domain = Domain::new_for_size(required_domain_size as u64)?;

        let mut domain_elements = materialize_domain_elements_with_natural_enumeration(
            &domain, 
            &worker
        );

        domain_elements.pop().expect("must pop last element for omega^i");

        let non_residues = make_non_residues::<E::Fr>(num_state_polys - 1);

        let mut domain_elements_poly_by_beta = Polynomial::from_values_unpadded(domain_elements)?;
        domain_elements_poly_by_beta.scale(&worker, beta_for_copy_permutation);

        // we take A, B, C, ... values and form (A + beta * X * non_residue + gamma), etc and calculate their grand product

        let mut z_num = {
            let mut grand_products_proto_it = grand_products_protos_with_gamma.iter().cloned();

            let mut z_1 = grand_products_proto_it.next().unwrap();
            z_1.add_assign(&worker, &domain_elements_poly_by_beta);

            for (mut p, non_res) in grand_products_proto_it.zip(non_residues.iter()) {
                p.add_assign_scaled(&worker, &domain_elements_poly_by_beta, non_res);
                z_1.mul_assign(&worker, &p);
            }

            z_1
        };

        // we take A, B, C, ... values and form (A + beta * perm_a + gamma), etc and calculate their grand product

        let mut permutation_polynomials_values_of_size_n_minus_one = vec![];
        
        for idx in 0..num_state_polys {
            let key = PolyIdentifier::PermutationPolynomial(idx);

            let mut coeffs = values_storage.get_poly(key).clone().into_coeffs();
            coeffs.pop().unwrap();

            let p = Polynomial::from_values_unpadded(coeffs)?;
            permutation_polynomials_values_of_size_n_minus_one.push(p);
        }

        let z_den = {
            assert_eq!(
                permutation_polynomials_values_of_size_n_minus_one.len(), 
                grand_products_protos_with_gamma.len()
            );
            let mut grand_products_proto_it = grand_products_protos_with_gamma.into_iter();
            let mut permutation_polys_it = permutation_polynomials_values_of_size_n_minus_one.iter();

            let mut z_2 = grand_products_proto_it.next().unwrap();
            z_2.add_assign_scaled(&worker, permutation_polys_it.next().unwrap(), &beta_for_copy_permutation);

            for (mut p, perm) in grand_products_proto_it
                                            .zip(permutation_polys_it) {
                // permutation polynomials 
                p.add_assign_scaled(&worker, &perm, &beta_for_copy_permutation);
                z_2.mul_assign(&worker, &p);
            }

            z_2.batch_inversion(&worker)?;

            z_2
        };

        z_num.mul_assign(&worker, &z_den);
        drop(z_den);

        let z = z_num.calculate_shifted_grand_product(&worker)?;
        drop(z_num);

        assert!(z.size().is_power_of_two());

        assert!(z.as_ref()[0] == E::Fr::one());

        let copy_permutation_z_in_monomial_form = z.ifft_using_bitreversed_ntt(
            &worker, 
            &omegas_inv_bitreversed, 
            &E::Fr::one()
        )?;

        // // TODO: move to setup
        // let mut permutation_polys_in_monomial_form = vec![];
        // for idx in 0.. {
        //     let key = PolyIdentifier::PermutationPolynomial(idx);
        //     let p = p.ifft_using_bitreversed_ntt(
        //         &worker, 
        //         &omegas_inv_bitreversed, 
        //         &E::Fr::one()
        //     )?;

        //     permutation_polys_in_monomial_form.push(p);
        // }

        let mut beta_for_lookup = None; 
        let mut gamma_for_lookup = None;

        let lookup_z_poly_in_monomial_form = if let Some(data) = lookup_data.as_mut() {
            let beta_for_lookup_permutation = E::Fr::from_str("789").unwrap();
            let gamma_for_lookup_permutation = E::Fr::from_str("1230").unwrap();

            beta_for_lookup = Some(beta_for_lookup_permutation);
            gamma_for_lookup = Some(gamma_for_lookup_permutation);
    
            let mut beta_plus_one = beta_for_lookup_permutation;
            beta_plus_one.add_assign(&E::Fr::one());
            let mut gamma_beta = gamma_for_lookup_permutation;
            gamma_beta.mul_assign(&beta_plus_one);
    
            let expected = gamma_beta.pow([(required_domain_size-1) as u64]);
            
            let f_poly_unpadded_values = data.f_poly_unpadded_values.take().unwrap();
            let t_poly_unpadded_values = data.t_poly_unpadded_values.take().unwrap();
            let t_shifted_unpadded_values = data.t_shifted_unpadded_values.take().unwrap();
            let s_poly_unpadded_values = data.s_poly_unpadded_values.take().unwrap();
            let s_shifted_unpadded_values = data.s_shifted_unpadded_values.take().unwrap();

            // Z(x*omega) = Z(x) * 
            // (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) / 
            // (\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) 

            let mut z_num = {
                // (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega))

                let mut t = t_poly_unpadded_values.clone();
                t.add_assign_scaled(&worker, &t_shifted_unpadded_values, &beta_for_lookup_permutation);
                t.add_constant(&worker, &gamma_beta);

                let mut tmp = f_poly_unpadded_values.clone();
                tmp.add_constant(&worker, &gamma_for_lookup_permutation);
                tmp.scale(&worker, beta_plus_one);

                t.mul_assign(&worker, &tmp);
                drop(tmp);

                t
            };

            let z_den = {
                // (\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))

                let mut t = s_poly_unpadded_values.clone();
                t.add_assign_scaled(&worker, &s_shifted_unpadded_values, &beta_for_lookup_permutation);
                t.add_constant(&worker, &gamma_beta);

                t.batch_inversion(&worker)?;

                t
            };

            z_num.mul_assign(&worker, &z_den);
            drop(z_den);

            let z = z_num.calculate_shifted_grand_product(&worker)?;
            drop(z_num);

            assert!(z.size().is_power_of_two());

            assert_eq!(z.as_ref()[0], E::Fr::one());
            assert_eq!(*z.as_ref().last().unwrap(), expected);

            let t_poly_monomial = t_poly_unpadded_values.clone_padded_to_domain()?.ifft_using_bitreversed_ntt(
                &worker,
                &omegas_inv_bitreversed,
                &E::Fr::one()
            )?;

            let s_poly_monomial = s_poly_unpadded_values.clone_padded_to_domain()?.ifft_using_bitreversed_ntt(
                &worker,
                &omegas_inv_bitreversed,
                &E::Fr::one()
            )?;

            // data.f_poly_monomial = Some(f_poly_monomial);
            data.t_poly_monomial = Some(t_poly_monomial);
            data.s_poly_monomial = Some(s_poly_monomial);

            let z = z.ifft_using_bitreversed_ntt(
                &worker, 
                &omegas_inv_bitreversed, 
                &E::Fr::one()
            )?;

            Some(z)
        } else {
            None
        };

        // TODO: commit to z polynomials

        // now draw alpha and add all the contributions to the quotient polynomial

        let alpha = E::Fr::from_str("1234567890").unwrap();
        
        let mut total_powers_of_alpha_for_gates = 0;
        for g in self.sorted_gates.iter() {
            total_powers_of_alpha_for_gates += g.num_quotient_terms();
        }

        println!("Have {} terms from {} gates", total_powers_of_alpha_for_gates, self.sorted_gates.len());

        let mut current_alpha = E::Fr::one();
        let mut powers_of_alpha_for_gates = Vec::with_capacity(total_powers_of_alpha_for_gates);
        powers_of_alpha_for_gates.push(current_alpha);
        for _ in 1..total_powers_of_alpha_for_gates {
            current_alpha.mul_assign(&alpha);
            powers_of_alpha_for_gates.push(current_alpha);
        }

        assert_eq!(powers_of_alpha_for_gates.len(), total_powers_of_alpha_for_gates);

        let mut all_gates = self.sorted_gates.clone();
        let num_different_gates = self.sorted_gates.len();

        let mut challenges_slice = &powers_of_alpha_for_gates[..];

        let mut lde_factor = num_state_polys;
        for g in self.sorted_gates.iter() {
            let degree = g.degree();
            if degree > lde_factor {
                lde_factor = degree;
            }
        }

        assert!(lde_factor <= 4);

        let coset_factor = E::Fr::multiplicative_generator();

        let mut t_poly = {
            let gate = all_gates.drain(0..1).into_iter().next().unwrap();
            assert!(<Self as ConstraintSystem<E>>::MainGate::default().into_internal() == gate);
            let gate = <Self as ConstraintSystem<E>>::MainGate::default();
            let num_challenges = gate.num_quotient_terms();
            let (for_gate, rest) = challenges_slice.split_at(num_challenges);
            challenges_slice = rest;

            let input_values = self.input_assingments.clone();

            let mut t = gate.contribute_into_quotient_for_public_inputs(
                required_domain_size,
                &input_values,
                &mut ldes_storage,
                &monomials_storage,
                for_gate,
                &omegas_bitreversed,
                &omegas_inv_bitreversed,
                &worker
            )?;

            if num_different_gates > 1 {
                // we have to multiply by the masking poly (selector)
                let key = PolyIdentifier::GateSelector(gate.name());
                let monomial_selector = monomials_storage.gate_selectors.get(&key).unwrap().as_ref();
                let lde = monomial_selector.clone_padded_to_domain()?.bitreversed_lde_using_bitreversed_ntt(
                    &worker, 
                    lde_factor, 
                    &omegas_bitreversed, 
                    &coset_factor
                )?;

                t.mul_assign(&worker, &lde);
                drop(lde);
            }

            t
        };

        let non_main_gates = all_gates;

        for gate in non_main_gates.into_iter() {
            let num_challenges = gate.num_quotient_terms();
            let (for_gate, rest) = challenges_slice.split_at(num_challenges);
            challenges_slice = rest;

            let mut contribution = gate.contribute_into_quotient(
                required_domain_size,
                &mut ldes_storage,
                &monomials_storage,
                for_gate,
                &omegas_bitreversed,
                &omegas_inv_bitreversed,
                &worker
            )?;

            {
                // we have to multiply by the masking poly (selector)
                let key = PolyIdentifier::GateSelector(gate.name());
                let monomial_selector = monomials_storage.gate_selectors.get(&key).unwrap().as_ref();
                let lde = monomial_selector.clone_padded_to_domain()?.bitreversed_lde_using_bitreversed_ntt(
                    &worker, 
                    lde_factor, 
                    &omegas_bitreversed, 
                    &coset_factor
                )?;

                contribution.mul_assign(&worker, &lde);
                drop(lde);
            }

            t_poly.add_assign(&worker, &contribution);
        }

        assert_eq!(challenges_slice.len(), 0);

        println!("Power of alpha for a start of normal permutation argument = {}", total_powers_of_alpha_for_gates);

        // perform copy-permutation argument

        // we precompute L_{0} here cause it's necessary for both copy-permutation and lookup permutation

        // z(omega^0) - 1 == 0
        let l_0 = calculate_lagrange_poly::<E::Fr>(&worker, required_domain_size.next_power_of_two(), 0)?;

        let l_0_coset_lde_bitreversed = l_0.bitreversed_lde_using_bitreversed_ntt(
            &worker, 
            lde_factor, 
            &omegas_bitreversed, 
            &coset_factor
        )?;

        let mut copy_grand_product_alphas = None;
        let x_poly_lde_bitreversed = {
            // now compute the permutation argument

            let z_coset_lde_bitreversed = copy_permutation_z_in_monomial_form.clone().bitreversed_lde_using_bitreversed_ntt(
                &worker, 
                lde_factor, 
                &omegas_bitreversed, 
                &coset_factor
            )?;

            assert!(z_coset_lde_bitreversed.size() == required_domain_size*lde_factor);

            let z_shifted_coset_lde_bitreversed = z_coset_lde_bitreversed.clone_shifted_assuming_bitreversed(
                lde_factor,
                &worker,
            )?;

            assert!(z_shifted_coset_lde_bitreversed.size() == required_domain_size*lde_factor);

            // For both Z_1 and Z_2 we first check for grand products
            // z*(X)(A + beta*X + gamma)(B + beta*k_1*X + gamma)(C + beta*K_2*X + gamma) - 
            // - (A + beta*perm_a(X) + gamma)(B + beta*perm_b(X) + gamma)(C + beta*perm_c(X) + gamma)*Z(X*Omega)== 0

            // we use evaluations of the polynomial X and K_i * X on a large domain's coset
            let mut contrib_z = z_coset_lde_bitreversed.clone();

            // precompute x poly
            let mut x_poly = Polynomial::from_values(vec![
                coset_factor;
                required_domain_size*lde_factor
            ])?;
            x_poly.distribute_powers(&worker, z_shifted_coset_lde_bitreversed.omega);
            x_poly.bitreverse_enumeration(&worker);
    
            assert_eq!(x_poly.size(), required_domain_size * lde_factor);

            // A + beta*X + gamma

            let mut tmp = ldes_storage.state_map.get(&PolyIdentifier::VariablesPolynomial(0)).unwrap().as_ref().clone();
            tmp.add_constant(&worker, &gamma_for_copy_permutation);
            tmp.add_assign_scaled(&worker, &x_poly, &beta_for_copy_permutation);
            contrib_z.mul_assign(&worker, &tmp);

            assert_eq!(non_residues.len() + 1, num_state_polys);

            for (poly_idx, non_res) in (1..num_state_polys).zip(non_residues.iter()) {
                let mut factor = beta_for_copy_permutation;
                factor.mul_assign(&non_res);

                let key = PolyIdentifier::VariablesPolynomial(poly_idx);
                tmp.reuse_allocation(&ldes_storage.state_map.get(&key).unwrap().as_ref());
                tmp.add_constant(&worker, &gamma_for_copy_permutation);
                tmp.add_assign_scaled(&worker, &x_poly, &factor);
                contrib_z.mul_assign(&worker, &tmp);
            }

            t_poly.add_assign_scaled(&worker, &contrib_z, &current_alpha);

            let alpha_0 = current_alpha;

            drop(contrib_z);

            let mut contrib_z = z_shifted_coset_lde_bitreversed;

            // A + beta*perm_a + gamma

            for idx in 0..num_state_polys {
                let key = PolyIdentifier::VariablesPolynomial(idx);

                tmp.reuse_allocation(&ldes_storage.state_map.get(&key).unwrap().as_ref());
                tmp.add_constant(&worker, &gamma_for_copy_permutation);

                let key = PolyIdentifier::PermutationPolynomial(idx);
                let perm = monomials_storage.get_poly(key).clone().bitreversed_lde_using_bitreversed_ntt(
                    &worker, 
                    lde_factor, 
                    &omegas_bitreversed, 
                    &coset_factor
                )?;
                tmp.add_assign_scaled(&worker, &perm, &beta_for_copy_permutation);
                contrib_z.mul_assign(&worker, &tmp);
                drop(perm);
            }

            t_poly.sub_assign_scaled(&worker, &contrib_z, &current_alpha);

            drop(contrib_z);

            drop(tmp);

            // Z(x) * L_{0}(x) - 1 == 0
            current_alpha.mul_assign(&alpha);

            let alpha_1 = current_alpha;

            {
                let mut z_minus_one_by_l_0 = z_coset_lde_bitreversed;
                z_minus_one_by_l_0.sub_constant(&worker, &E::Fr::one());

                z_minus_one_by_l_0.mul_assign(&worker, &l_0_coset_lde_bitreversed);

                t_poly.add_assign_scaled(&worker, &z_minus_one_by_l_0, &current_alpha);
            }

            copy_grand_product_alphas = Some([alpha_0, alpha_1]);

            x_poly
        };

        // let inverse_divisor_on_coset_lde_natural_ordering = {
        //     let mut vanishing_poly_inverse_bitreversed =
        //         evaluate_vanishing_polynomial_of_degree_on_domain_size::<E::Fr>(
        //             required_domain_size as u64,
        //             &E::Fr::multiplicative_generator(),
        //             (required_domain_size * lde_factor) as u64,
        //             &worker,
        //         )?;
        //     vanishing_poly_inverse_bitreversed.batch_inversion(&worker)?;
        //     // vanishing_poly_inverse_bitreversed.bitreverse_enumeration(&worker)?;

        //     vanishing_poly_inverse_bitreversed
        // };

        // t_poly.bitreverse_enumeration(&worker);

        // t_poly.mul_assign(&worker, &inverse_divisor_on_coset_lde_natural_ordering);

        // drop(inverse_divisor_on_coset_lde_natural_ordering);

        // add contribution from grand product for loopup polys if there is one

        let mut lookup_grand_product_alphas = None;
        if let Some(z_poly_in_monomial_form) = lookup_z_poly_in_monomial_form.as_ref() {
            let beta_for_lookup_permutation = beta_for_lookup.unwrap();
            let gamma_for_lookup_permutation = gamma_for_lookup.unwrap();

            let mut beta_plus_one = beta_for_lookup_permutation;
            beta_plus_one.add_assign(&E::Fr::one());
            let mut gamma_beta = gamma_for_lookup_permutation;
            gamma_beta.mul_assign(&beta_plus_one);
    
            let expected = gamma_beta.pow([(required_domain_size-1) as u64]);

            current_alpha.mul_assign(&alpha);

            let alpha_0 = current_alpha;

            // same grand product argument for lookup permutation except divisor is now with one point cut

            let z_lde = z_poly_in_monomial_form.clone().bitreversed_lde_using_bitreversed_ntt(
                &worker, 
                lde_factor, 
                &omegas_bitreversed, 
                &coset_factor
            )?;

            let z_lde_shifted = z_lde.clone_shifted_assuming_bitreversed(
                lde_factor,
                &worker
            )?;

            // We make an small ad-hoc modification here and instead of dividing some contributions by
            // (X^n - 1)/(X - omega^{n-1}) we move (X - omega^{n-1}) to the numerator and join the divisions

            // Numerator degree is at max 4n, so it's < 4n after division

            // ( Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -  
            // - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) )*(X - omega^{n-1})

            let data = lookup_data.as_ref().unwrap();

            let s_lde = data.s_poly_monomial.as_ref().unwrap().clone().bitreversed_lde_using_bitreversed_ntt(
                &worker, 
                lde_factor, 
                &omegas_bitreversed, 
                &coset_factor
            )?;

            let s_lde_shifted = s_lde.clone_shifted_assuming_bitreversed(
                lde_factor,
                &worker
            )?;

            // Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))

            let mut contribution = s_lde;
            contribution.add_assign_scaled(&worker, &s_lde_shifted, &beta_for_lookup_permutation);
            contribution.add_constant(&worker, &gamma_beta);
            contribution.mul_assign(&worker, &z_lde_shifted);

            drop(s_lde_shifted);
            drop(z_lde_shifted);

            let t_lde = data.t_poly_monomial.as_ref().unwrap().clone().bitreversed_lde_using_bitreversed_ntt(
                &worker, 
                lde_factor, 
                &omegas_bitreversed, 
                &coset_factor
            )?;

            let t_lde_shifted = t_lde.clone_shifted_assuming_bitreversed(
                lde_factor,
                &worker
            )?;

            let f_lde = {
                // add up ldes of a,b,c and table_type poly and multiply by selector
                
                let a_ref = get_from_map_unchecked(
                    PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
                    &ldes_storage
                );
                let mut tmp = a_ref.clone();
                drop(a_ref);

                let eta = lookup_data.as_ref().unwrap().eta;

                let mut current = eta;

                let b_ref = get_from_map_unchecked(
                    PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)),
                    &ldes_storage
                );

                tmp.add_assign_scaled(&worker, b_ref, &current);

                drop(b_ref);
                current.mul_assign(&eta);

                let c_ref = get_from_map_unchecked(
                    PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)),
                    &ldes_storage
                );

                tmp.add_assign_scaled(&worker, c_ref, &current);

                drop(c_ref);
                current.mul_assign(&eta);

                let table_type_lde = lookup_data.as_ref().unwrap().table_type_poly_monomial.as_ref().unwrap().clone()
                    .bitreversed_lde_using_bitreversed_ntt(
                        &worker, 
                        lde_factor, 
                        &omegas_bitreversed, 
                        &coset_factor
                )?;

                tmp.add_assign_scaled(&worker, &table_type_lde, &current);

                drop(table_type_lde);

                let lookup_selector_lde = lookup_data.as_ref().unwrap().selector_poly_monomial.as_ref().unwrap().clone()
                    .bitreversed_lde_using_bitreversed_ntt(
                        &worker, 
                        lde_factor, 
                        &omegas_bitreversed, 
                        &coset_factor
                )?;

                tmp.mul_assign(&worker, &lookup_selector_lde);

                drop(lookup_selector_lde);

                tmp

            };

            //  - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega))

            let mut tmp = f_lde;
            tmp.add_constant(&worker, &gamma_for_lookup_permutation);
            tmp.mul_assign(&worker, &z_lde);
            tmp.scale(&worker, beta_plus_one);

            let mut t = t_lde;
            t.add_assign_scaled(&worker, &t_lde_shifted, &beta_for_lookup_permutation);
            t.add_constant(&worker, &gamma_beta);

            tmp.mul_assign(&worker, &t);
            
            drop(t);
            drop(t_lde_shifted);

            contribution.sub_assign(&worker, &tmp);

            contribution.scale(&worker, current_alpha);

            // multiply by (X - omega^{n-1})

            let last_omega = domain.generator.pow(&[(required_domain_size - 1) as u64]);
            let mut x_minus_last_omega = x_poly_lde_bitreversed;
            x_minus_last_omega.sub_constant(&worker, &last_omega);

            contribution.mul_assign(&worker, &x_minus_last_omega);
            drop(x_minus_last_omega);

            // we do not need to do addition multiplications for terms below cause multiplication by lagrange poly
            // does everything for us

            // check that (Z(x) - 1) * L_{0} == 0
            current_alpha.mul_assign(&alpha);

            let alpha_1 = current_alpha;
            
            tmp.reuse_allocation(&z_lde);
            tmp.sub_constant(&worker, &E::Fr::one());
            tmp.mul_assign(&worker, &l_0_coset_lde_bitreversed);

            drop(l_0_coset_lde_bitreversed);

            contribution.add_assign_scaled(&worker, &tmp, &current_alpha);

            // check that (Z(x) - expected) * L_{n-1}  == 0

            current_alpha.mul_assign(&alpha);

            let alpha_2 = current_alpha;

            let l_last = calculate_lagrange_poly::<E::Fr>(&worker, required_domain_size.next_power_of_two(), required_domain_size - 1)?;

            let l_last_coset_lde_bitreversed = l_last.bitreversed_lde_using_bitreversed_ntt(
                &worker, 
                lde_factor, 
                &omegas_bitreversed, 
                &coset_factor
            )?;

            tmp.reuse_allocation(&z_lde);
            tmp.sub_constant(&worker, &expected);
            tmp.mul_assign(&worker, &l_last_coset_lde_bitreversed);

            drop(l_last_coset_lde_bitreversed);

            contribution.add_assign_scaled(&worker, &tmp, &current_alpha);

            drop(tmp);
            drop(z_lde);

            // contribution.bitreverse_enumeration(&worker);


            // let vanishing_for_lookup = {
            //     let normally_enumerated = calculate_inverse_vanishing_polynomial_with_last_point_cut(
            //         &worker,
            //         required_domain_size * lde_factor,
            //         required_domain_size,
            //         coset_factor,
            //     )?;

            //     normally_enumerated
            // };

            // contribution.mul_assign(&worker, &vanishing_for_lookup);

            t_poly.add_assign(&worker, &contribution);

            drop(contribution);

            lookup_grand_product_alphas = Some([alpha_0, alpha_1, alpha_2]);
        } else {
            drop(x_poly_lde_bitreversed);
            drop(l_0_coset_lde_bitreversed);
        }

        // perform the division

        let inverse_divisor_on_coset_lde_natural_ordering = {
            let mut vanishing_poly_inverse_bitreversed =
                evaluate_vanishing_polynomial_of_degree_on_domain_size::<E::Fr>(
                    required_domain_size as u64,
                    &E::Fr::multiplicative_generator(),
                    (required_domain_size * lde_factor) as u64,
                    &worker,
                )?;
            vanishing_poly_inverse_bitreversed.batch_inversion(&worker)?;
            // vanishing_poly_inverse_bitreversed.bitreverse_enumeration(&worker)?;

            vanishing_poly_inverse_bitreversed
        };

        // don't forget to bitreverse

        t_poly.bitreverse_enumeration(&worker);

        t_poly.mul_assign(&worker, &inverse_divisor_on_coset_lde_natural_ordering);

        drop(inverse_divisor_on_coset_lde_natural_ordering);
        
        let t_poly = t_poly.icoset_fft_for_generator(&worker, &coset_factor);

        println!("Lde factor = {}", lde_factor);

        // println!("Quotient poly = {:?}", t_poly.as_ref());

        {
            // degree is 4n-4
            let l = t_poly.as_ref().len();
            assert_eq!(&t_poly.as_ref()[(l-4)..], &[E::Fr::zero(); 4][..]);
        }

        println!("Quotient poly degree = {}", get_degree::<E::Fr>(&t_poly));

        let t_poly_parts = t_poly.break_into_multiples(required_domain_size)?;
        // TODO: commit every of them

        // draw opening point
        let z = E::Fr::from_str("333444555").unwrap();
        let omega = domain.generator;

        // Now perform the linearization.
        // First collect and evalute all the polynomials that are necessary for linearization

        let mut opening_requests_before_linearization = std::collections::HashSet::new();
        let mut all_queries = std::collections::HashSet::new();
        let mut sorted_opening_requests = vec![];
        let mut sorted_selector_for_opening = vec![];
        let mut polys_in_linearization = std::collections::HashSet::new();

        let num_gate_types = self.sorted_gates.len();

        for gate in self.sorted_gates.iter() {
            for q in gate.all_queried_polynomials().into_iter() {
                all_queries.insert(q);
            }
            if gate.benefits_from_linearization() {
                // it's better to linearize the gate
                if num_gate_types > 1 {
                    // there are various gates, so we need to query the selector
                    sorted_selector_for_opening.push(gate.box_clone());
                }

                for q in gate.needs_opened_for_linearization().into_iter() {
                    if !opening_requests_before_linearization.contains(&q) {
                        opening_requests_before_linearization.insert(q.clone());
                        sorted_opening_requests.push(q);
                    }
                }

                for q in gate.linearizes_over().into_iter() {
                    polys_in_linearization.insert(q);
                }
            } else {
                // we will linearize over the selector, so we do not need to query it
            }
        }

        // println!("Will need to open for linearization: {:?}", sorted_opening_requests);
        // println!("Will open selectors: {:?}", sorted_selector_for_opening);

        // println!("All queries = {:?}", all_queries);

        // Sanity check: we open everything either in linearization or in plain text! 
        {
            let must_open_without_linearization: Vec<_> = all_queries.difference(&polys_in_linearization).collect();

            for p in must_open_without_linearization.into_iter() {
                assert!(opening_requests_before_linearization.contains(&p));
            }
        }

        let mut query_values = vec![];
        let mut query_values_map_over_dilation = std::collections::HashMap::new();
        let mut query_values_map = std::collections::HashMap::new();

        for p in sorted_opening_requests.into_iter() {
            let (id, dilation_value) = p.into_id_and_raw_dilation();
            let poly_ref = match id {
                p @ PolyIdentifier::VariablesPolynomial(..) => {
                    monomials_storage.state_map.get(&p).unwrap().as_ref()
                },
                p @ PolyIdentifier::WitnessPolynomial(..) => {
                    monomials_storage.witness_map.get(&p).unwrap().as_ref()
                },
                p @ PolyIdentifier::GateSetupPolynomial(..) => {
                    monomials_storage.setup_map.get(&p).unwrap().as_ref()
                },
                _ => {
                    unreachable!();
                }
            };

            let mut opening_point = z;
            for _ in 0..dilation_value {
                opening_point.mul_assign(&omega);
            }

            let value = poly_ref.evaluate_at(&worker, opening_point);

            let entry = query_values_map_over_dilation.entry(dilation_value).or_insert(vec![]);
            entry.push(value);

            query_values_map.insert(p, value);
            query_values.push(value);
        }

        // TODO: commit queried values

        // println!("Queried values = {:?}", query_values_map_over_dilation);

        // also open selectors

        let mut selector_values = vec![];
        for s in sorted_selector_for_opening.into_iter() {
            let key = PolyIdentifier::GateSelector(s.name());
            let poly_ref = monomials_storage.gate_selectors.get(&key).unwrap().as_ref();
            let value = poly_ref.evaluate_at(&worker, z);

            selector_values.push(value);
        }

        // TODO: commit queried selector values

        // copy-permutation polynomials queries

        let mut copy_permutation_queries = vec![];

        for idx in 0..(num_state_polys-1) {
            let key = PolyIdentifier::PermutationPolynomial(idx);
            let value = monomials_storage.get_poly(key).evaluate_at(&worker, z);

            copy_permutation_queries.push(value);
        }

        // TODO: commit queried permutation poly values

        // copy-permutation grand product query

        let mut z_omega = z;
        z_omega.mul_assign(&domain.generator);
        let copy_permutation_z_at_z_omega = copy_permutation_z_in_monomial_form.evaluate_at(&worker, z_omega);

        // TODO: commit copy-permutation grand product query

        // we've computed everything, so perform linearization

        let mut challenges_slice = &powers_of_alpha_for_gates[..];

        let mut all_gates = self.sorted_gates.clone();

        let mut r_poly = {
            let gate = all_gates.drain(0..1).into_iter().next().unwrap();
            assert!(gate.benefits_from_linearization(), "main gate is expected to benefit from linearization!");
            assert!(<Self as ConstraintSystem<E>>::MainGate::default().into_internal() == gate);
            let gate = <Self as ConstraintSystem<E>>::MainGate::default();
            let num_challenges = gate.num_quotient_terms();
            let (for_gate, rest) = challenges_slice.split_at(num_challenges);
            challenges_slice = rest;

            let input_values = self.input_assingments.clone();

            let mut r = gate.contribute_into_linearization_for_public_inputs(
                required_domain_size,
                &input_values,
                z,
                &query_values_map,
                &monomials_storage,
                for_gate,
                &worker
            )?;

            let mut selectors_it = selector_values.clone().into_iter();

            if num_different_gates > 1 {
                // first multiply r by the selector value at z
                r.scale(&worker, selectors_it.next().unwrap());
            }

            // now proceed per gate
            for gate in all_gates.into_iter() {
                let num_challenges = gate.num_quotient_terms();
                let (for_gate, rest) = challenges_slice.split_at(num_challenges);
                challenges_slice = rest;

                if gate.benefits_from_linearization() {
                    // gate benefits from linearization, so make temporary value
                    let tmp = gate.contribute_into_linearization(
                        required_domain_size,
                        z,
                        &query_values_map,
                        &monomials_storage,
                        for_gate,
                        &worker
                    )?;

                    let selector_value = selectors_it.next().unwrap();

                    r.add_assign_scaled(&worker, &tmp, &selector_value);
                } else {
                    // we linearize over the selector, so take a selector and scale it
                    let gate_value_at_z = gate.contribute_into_verification_equation(
                        required_domain_size,
                        z,
                        &query_values_map,
                        for_gate
                    )?;

                    let key = PolyIdentifier::GateSelector(gate.name());
                    let gate_selector_ref = monomials_storage.gate_selectors.get(&key).expect("must get monomial form of gate selector").as_ref();

                    r.add_assign_scaled(&worker, gate_selector_ref, &gate_value_at_z);
                }
            }

            assert!(selectors_it.next().is_none());
            assert_eq!(challenges_slice.len(), 0);

            r
        };

        // add contributions from copy-permutation and lookup-permutation

        // copy-permutation linearization comtribution
        {
            // + (a(z) + beta*z + gamma)*()*()*()*Z(x)

            let [alpha_0, alpha_1] = copy_grand_product_alphas.expect("there must be powers of alpha for copy permutation");

            let some_one = Some(E::Fr::one());
            let mut non_residues_iterator = some_one.iter().chain(&non_residues);

            let mut factor = alpha_0;

            for idx in 0..num_state_polys {
                let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(idx));
                let wire_value = query_values_map.get(&key)
                    .ok_or(SynthesisError::AssignmentMissing)?;
                let mut t = z;
                let non_res = non_residues_iterator.next().unwrap();
                t.mul_assign(&non_res);
                t.mul_assign(&beta_for_copy_permutation);
                t.add_assign(&wire_value);
                t.add_assign(&gamma_for_copy_permutation);

                factor.mul_assign(&t);
            }

            assert!(non_residues_iterator.next().is_none());

            r_poly.add_assign_scaled(&worker, &copy_permutation_z_in_monomial_form, &factor);

            // - (a(z) + beta*perm_a + gamma)*()*()*z(z*omega) * beta * perm_d(X)

            let mut factor = alpha_0;
            factor.mul_assign(&beta_for_copy_permutation);
            factor.mul_assign(&copy_permutation_z_at_z_omega);

            for idx in 0..(num_state_polys-1) {
                let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(idx));
                let wire_value = query_values_map.get(&key)
                    .ok_or(SynthesisError::AssignmentMissing)?;
                let permutation_at_z = copy_permutation_queries[idx];
                let mut t = permutation_at_z;
                
                t.mul_assign(&beta_for_copy_permutation);
                t.add_assign(&wire_value);
                t.add_assign(&gamma_for_copy_permutation);

                factor.mul_assign(&t);
            }

            let key = PolyIdentifier::PermutationPolynomial(num_state_polys - 1);
            let last_permutation_poly_ref = monomials_storage.get_poly(key);

            r_poly.sub_assign_scaled(&worker, last_permutation_poly_ref, &factor);

            // + L_0(z) * Z(x)

            let mut factor = evaluate_l0_at_point(required_domain_size as u64, z)?;
            factor.mul_assign(&alpha_1);

            r_poly.add_assign_scaled(&worker, &copy_permutation_z_in_monomial_form, &factor);
        }

        // lookup grand product linearization

        // due to separate divisor it's not obvious if this is beneficial without some tricks
        // like multiplication by (1 - L_{n-1}) or by (x - omega^{n-1})

        // Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -  
        // Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) == 0
        // check that (Z(x) - 1) * L_{0} == 0
        // check that (Z(x) - expected) * L_{n-1} == 0, or (Z(x*omega) - expected)* L_{n-2} == 0

        // f(x) does not need to be opened as it's made of table selector and witnesses
        // if we pursue the strategy from the linearization of a copy-permutation argument
        // then we leave something like s(x) from the Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) term,
        // and Z(x) from Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) term,
        // with terms with lagrange polys as multipliers left intact

        let lookup_queries = if let Some(lookup_z_poly) = lookup_z_poly_in_monomial_form.as_ref() {
            let [alpha_0, alpha_1, alpha_2] = lookup_grand_product_alphas.expect("there must be powers of alpha for lookup permutation");

            let s_at_z_omega = lookup_data.as_ref().unwrap().s_poly_monomial.as_ref().unwrap().evaluate_at(&worker, z_omega);
            let grand_product_at_z_omega = lookup_z_poly.evaluate_at(&worker, z_omega);
            let t_at_z = lookup_data.as_ref().unwrap().t_poly_monomial.as_ref().unwrap().evaluate_at(&worker, z);
            let t_at_z_omega = lookup_data.as_ref().unwrap().t_poly_monomial.as_ref().unwrap().evaluate_at(&worker, z_omega);
            let selector_at_z = lookup_data.as_ref().unwrap().selector_poly_monomial.as_ref().unwrap().evaluate_at(&worker, z);
            let table_type_at_z = lookup_data.as_ref().unwrap().table_type_poly_monomial.as_ref().unwrap().evaluate_at(&worker, z);

            let l_0_at_z = evaluate_lagrange_poly_at_point(0, &domain, z)?;
            let l_n_minus_one_at_z = evaluate_lagrange_poly_at_point(required_domain_size - 1, &domain, z)?;

            let beta_for_lookup_permutation = beta_for_lookup.unwrap();
            let gamma_for_lookup_permutation = gamma_for_lookup.unwrap();

            let mut beta_plus_one = beta_for_lookup_permutation;
            beta_plus_one.add_assign(&E::Fr::one());
            let mut gamma_beta = gamma_for_lookup_permutation;
            gamma_beta.mul_assign(&beta_plus_one);

            // (Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -  
            // Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)))*(X - omega^{n-1})

            let last_omega = domain.generator.pow(&[(required_domain_size - 1) as u64]);
            let mut z_minus_last_omega = z;
            z_minus_last_omega.sub_assign(&last_omega);

            // s(x) from the Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))
            let mut factor = grand_product_at_z_omega; // we do not need to account for additive terms
            factor.mul_assign(&alpha_0);
            factor.mul_assign(&z_minus_last_omega);

            r_poly.add_assign_scaled(&worker, lookup_data.as_ref().unwrap().s_poly_monomial.as_ref().unwrap(), &factor);

            // Z(x) from - alpha_0 * Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) 
            // + alpha_1 * Z(x) * L_{0}(z) + alpha_2 * Z(x) * L_{n-1}(z)

            // accumulate coefficient
            let mut factor = t_at_z_omega;
            factor.mul_assign(&beta_for_lookup_permutation);
            factor.add_assign(&t_at_z);
            factor.add_assign(&gamma_beta);

            // (\gamma + f(x))

            let mut f_reconstructed = E::Fr::zero();
            let mut current = E::Fr::one();
            let eta = lookup_data.as_ref().unwrap().eta;
            // a,b,c
            for idx in 0..(num_state_polys-1) {
                let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(idx));
                let mut value = *query_values_map.get(&key)
                    .ok_or(SynthesisError::AssignmentMissing)?;

                value.mul_assign(&current);
                f_reconstructed.add_assign(&value);

                current.mul_assign(&eta);
            }

            // and table type
            let mut t = table_type_at_z;
            t.mul_assign(&current);
            f_reconstructed.add_assign(&t);

            f_reconstructed.mul_assign(&selector_at_z);
            f_reconstructed.add_assign(&gamma_for_lookup_permutation);

            // end of (\gamma + f(x)) part

            factor.mul_assign(&f_reconstructed);
            factor.mul_assign(&beta_plus_one);
            factor.negate(); // don't forget minus sign
            factor.mul_assign(&alpha_0);

            // Multiply by (z - omega^{n-1})

            factor.mul_assign(&z_minus_last_omega);

            // L_{0}(z) in front of Z(x)

            let mut tmp = l_0_at_z;
            tmp.mul_assign(&alpha_1);
            factor.add_assign(&tmp);

            // L_{n-1}(z) in front of Z(x)

            let mut tmp = l_n_minus_one_at_z;
            tmp.mul_assign(&alpha_2);
            factor.add_assign(&tmp);

            r_poly.add_assign_scaled(&worker, lookup_z_poly, &factor);

            let query = LookupQuery::<E> {
                s_at_z_omega,
                grand_product_at_z_omega,
                t_at_z,
                t_at_z_omega,
                selector_at_z,
                table_type_at_z,
            };

            Some(query)
        } else {
            None
        };

        let linearization_at_z = r_poly.evaluate_at(&worker, z);

        // don't forget to evaluate quotient too

        let t_at_z = {
            let mut result = E::Fr::zero();
            let mut current = E::Fr::one();
            let z_in_domain_size = z.pow(&[required_domain_size as u64]);
            for p in t_poly_parts.iter() {
                let mut subvalue_at_z = p.evaluate_at(&worker, z);
                subvalue_at_z.mul_assign(&current);
                result.add_assign(&subvalue_at_z);
                current.mul_assign(&z_in_domain_size);
            }

            result
        };

        // linearization is done, now perform sanity check
        // this is effectively a verification procedure

        {
            let vanishing_at_z = evaluate_vanishing_for_size(&z, required_domain_size as u64);

            // first let's aggregate gates

            let mut t_num_on_full_domain = E::Fr::zero();

            let challenges_slice = &powers_of_alpha_for_gates[..];

            let mut all_gates = self.sorted_gates.clone();

            // we've suffered and linearization polynomial captures all the gates except the public input!

            {
                let mut tmp = linearization_at_z;
                // add input values

                let gate = all_gates.drain(0..1).into_iter().next().unwrap();
                assert!(gate.benefits_from_linearization(), "main gate is expected to benefit from linearization!");
                assert!(<Self as ConstraintSystem<E>>::MainGate::default().into_internal() == gate);
                let gate = <Self as ConstraintSystem<E>>::MainGate::default();
                let num_challenges = gate.num_quotient_terms();
                let (for_gate, _) = challenges_slice.split_at(num_challenges);

                let input_values = self.input_assingments.clone();

                let mut inputs_term = gate.add_inputs_into_quotient(
                    required_domain_size,
                    &input_values,
                    z,
                    for_gate,
                )?;

                if num_different_gates > 1 {
                    let selector_value = selector_values[0];
                    inputs_term.mul_assign(&selector_value);
                }

                tmp.add_assign(&inputs_term);

                t_num_on_full_domain.add_assign(&tmp);
            } 

            // now aggregate leftovers from grand product for copy permutation
            {
                // - alpha_0 * (a + perm(z) * beta + gamma)*()*(d + gamma) * z(z*omega)
                let [alpha_0, alpha_1] = copy_grand_product_alphas.expect("there must be powers of alpha for copy permutation");

                let mut factor = alpha_0;
                factor.mul_assign(&copy_permutation_z_at_z_omega);

                for idx in 0..(num_state_polys-1) {
                    let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(idx));
                    let wire_value = query_values_map.get(&key)
                        .ok_or(SynthesisError::AssignmentMissing)?;
                    let permutation_at_z = copy_permutation_queries[idx];
                    let mut t = permutation_at_z;
                    
                    t.mul_assign(&beta_for_copy_permutation);
                    t.add_assign(&wire_value);
                    t.add_assign(&gamma_for_copy_permutation);

                    factor.mul_assign(&t);
                }

                let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(num_state_polys-1));
                let mut tmp = *query_values_map.get(&key)
                    .ok_or(SynthesisError::AssignmentMissing)?;
                tmp.add_assign(&gamma_for_copy_permutation);

                factor.mul_assign(&tmp);

                t_num_on_full_domain.sub_assign(&factor);

                // - L_0(z) * alpha_1

                let mut l_0_at_z = evaluate_l0_at_point(required_domain_size as u64, z)?;
                l_0_at_z.mul_assign(&alpha_1);
    
                t_num_on_full_domain.sub_assign(&l_0_at_z);
            }

            // and if exists - grand product for lookup permutation

            {
                if lookup_queries.is_some() {

                    let [alpha_0, alpha_1, alpha_2] = lookup_grand_product_alphas.expect("there must be powers of alpha for lookup permutation");

                    let lookup_queries = lookup_queries.clone().expect("lookup queries must be made");

                    let beta_for_lookup_permutation = beta_for_lookup.unwrap();
                    let gamma_for_lookup_permutation = gamma_for_lookup.unwrap();
                    let mut beta_plus_one = beta_for_lookup_permutation;
                    beta_plus_one.add_assign(&E::Fr::one());
                    let mut gamma_beta = gamma_for_lookup_permutation;
                    gamma_beta.mul_assign(&beta_plus_one);
            
                    let expected = gamma_beta.pow([(required_domain_size-1) as u64]);

                    // in a linearization we've taken terms:
                    // - s(x) from the alpha_0 * Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))
                    // - and Z(x) from - alpha_0 * Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) (term in full) + 
                    // + alpha_1 * (Z(x) - 1) * L_{0}(z) + alpha_2 * (Z(x) - expected) * L_{n-1}(z)

                    // first make alpha_0 * Z(x*omega)*(\gamma*(1 + \beta) + \beta * s(x*omega)))

                    let mut tmp = lookup_queries.s_at_z_omega;
                    tmp.mul_assign(&beta_for_lookup_permutation);
                    tmp.add_assign(&gamma_beta);
                    tmp.mul_assign(&lookup_queries.grand_product_at_z_omega);
                    tmp.mul_assign(&alpha_0);

                    // (z - omega^{n-1}) for this part
                    let last_omega = domain.generator.pow(&[(required_domain_size - 1) as u64]);
                    let mut z_minus_last_omega = z;
                    z_minus_last_omega.sub_assign(&last_omega);

                    tmp.mul_assign(&z_minus_last_omega);

                    t_num_on_full_domain.add_assign(&tmp);

                    // // - alpha_1 * L_{0}(z)

                    let mut l_0_at_z = evaluate_l0_at_point(required_domain_size as u64, z)?;
                    l_0_at_z.mul_assign(&alpha_1);
        
                    t_num_on_full_domain.sub_assign(&l_0_at_z);

                    // // - alpha_2 * expected L_{n-1}(z)

                    let mut l_n_minus_one_at_z = evaluate_lagrange_poly_at_point(required_domain_size - 1, &domain, z)?;
                    l_n_minus_one_at_z.mul_assign(&expected);
                    l_n_minus_one_at_z.mul_assign(&alpha_2);

                    t_num_on_full_domain.sub_assign(&l_n_minus_one_at_z);
                }
            }

            let mut lhs = t_at_z;
            lhs.mul_assign(&vanishing_at_z);

            let rhs = t_num_on_full_domain;

            if lhs != rhs {
                dbg!("Circuit is not satisfied");
                return Err(SynthesisError::Unsatisfiable);
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::pairing::Engine;
    use crate::pairing::ff::PrimeField;

    struct TestCircuit4<E:Engine>{
        _marker: PhantomData<E>
    }

    impl<E: Engine> Circuit<E> for TestCircuit4<E> {
        type MainGate = Width4MainGateWithDNext;

        fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
            let a = cs.alloc(|| {
                Ok(E::Fr::from_str("10").unwrap())
            })?;

            println!("A = {:?}", a);

            let b = cs.alloc(|| {
                Ok(E::Fr::from_str("20").unwrap())
            })?;

            println!("B = {:?}", b);

            let c = cs.alloc(|| {
                Ok(E::Fr::from_str("200").unwrap())
            })?;

            println!("C = {:?}", c);

            let d = cs.alloc(|| {
                Ok(E::Fr::from_str("100").unwrap())
            })?;

            println!("D = {:?}", d);

            let one = E::Fr::one();

            let mut two = one;
            two.double();

            let mut negative_one = one;
            negative_one.negate();

            // 2a - b = 0

            let two_a = ArithmeticTerm::from_variable_and_coeff(a, two);
            let minus_b = ArithmeticTerm::from_variable_and_coeff(b, negative_one);
            let mut term = MainGateTerm::new();
            term.add_assign(two_a);
            term.add_assign(minus_b);

            cs.allocate_main_gate(term)?;

            // c - a*b == 0 

            let mut ab_term = ArithmeticTerm::from_variable(a).mul_by_variable(b);
            ab_term.scale(&negative_one);
            let c_term = ArithmeticTerm::from_variable(c);
            let mut term = MainGateTerm::new();
            term.add_assign(c_term);
            term.add_assign(ab_term);

            cs.allocate_main_gate(term)?;

            // d - 100 == 0 

            let hundred = ArithmeticTerm::constant(E::Fr::from_str("100").unwrap());
            let d_term = ArithmeticTerm::from_variable(d);
            let mut term = MainGateTerm::new();
            term.add_assign(d_term);
            term.sub_assign(hundred);

            cs.allocate_main_gate(term)?;

            // let gamma = cs.alloc_input(|| {
            //     Ok(E::Fr::from_str("20").unwrap())
            // })?;

            let gamma = cs.alloc(|| {
                Ok(E::Fr::from_str("20").unwrap())
            })?;

            // gamma - b == 0 

            let gamma_term = ArithmeticTerm::from_variable(gamma);
            let b_term = ArithmeticTerm::from_variable(b);
            let mut term = MainGateTerm::new();
            term.add_assign(gamma_term);
            term.sub_assign(b_term);

            cs.allocate_main_gate(term)?;

            // 2a
            let mut term = MainGateTerm::<E>::new();
            term.add_assign(ArithmeticTerm::from_variable_and_coeff(a, two));

            let dummy = CS::get_dummy_variable();

            // 2a - d_next = 0

            let (vars, mut coeffs) = CS::MainGate::format_term(term, dummy)?;
            *coeffs.last_mut().unwrap() = negative_one;

            // here d is equal = 2a, so we need to place b there
            // and compensate it with -b somewhere before

            cs.new_single_gate_for_trace_step(&CS::MainGate::default(), 
                &coeffs, 
                &vars, 
                &[]
            )?;

            let mut term = MainGateTerm::<E>::new();
            term.add_assign(ArithmeticTerm::from_variable(b));

            // b + 0 + 0 - b = 0
            let (mut vars, mut coeffs) = CS::MainGate::format_term(term, dummy)?;
            coeffs[3] = negative_one;
            vars[3] = b;

            cs.new_single_gate_for_trace_step(&CS::MainGate::default(), 
                &coeffs, 
                &vars, 
                &[]
            )?;

            Ok(())
        }
    }

    struct TestCircuit4WithLookups<E:Engine>{
        _marker: PhantomData<E>
    }

    impl<E: Engine> Circuit<E> for TestCircuit4WithLookups<E> {
        type MainGate = Width4MainGateWithDNext;

        fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
            Ok(
                vec![
                    Width4MainGateWithDNext::default().into_internal(),
                    TestBitGate::default().into_internal()
                ]
            )
        }

        fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
            let columns = vec![PolyIdentifier::VariablesPolynomial(0), PolyIdentifier::VariablesPolynomial(1), PolyIdentifier::VariablesPolynomial(2)];
            let range_table = LookupTableApplication::new_range_table_of_width_3(2, columns.clone())?;
            let range_table_name = range_table.functional_name();

            let xor_table = LookupTableApplication::new_xor_table(2, columns.clone())?;
            let xor_table_name = xor_table.functional_name();

            let and_table = LookupTableApplication::new_and_table(2, columns)?;
            let and_table_name = and_table.functional_name();

            cs.add_table(range_table)?;
            cs.add_table(xor_table)?;
            cs.add_table(and_table)?;

            let a = cs.alloc(|| {
                Ok(E::Fr::from_str("10").unwrap())
            })?;

            println!("A = {:?}", a);

            let b = cs.alloc(|| {
                Ok(E::Fr::from_str("20").unwrap())
            })?;

            println!("B = {:?}", b);

            let c = cs.alloc(|| {
                Ok(E::Fr::from_str("200").unwrap())
            })?;

            println!("C = {:?}", c);

            let d = cs.alloc(|| {
                Ok(E::Fr::from_str("100").unwrap())
            })?;

            println!("D = {:?}", d);

            let e = cs.alloc(|| {
                Ok(E::Fr::from_str("2").unwrap())
            })?;

            let binary_x_value = E::Fr::from_str("3").unwrap();
            let binary_y_value = E::Fr::from_str("1").unwrap();

            let binary_x = cs.alloc(|| {
                Ok(binary_x_value)
            })?;

            let binary_y = cs.alloc(|| {
                Ok(binary_y_value)
            })?;

            let one = E::Fr::one();

            let mut two = one;
            two.double();

            let mut negative_one = one;
            negative_one.negate();

            // 2a - b = 0

            let two_a = ArithmeticTerm::from_variable_and_coeff(a, two);
            let minus_b = ArithmeticTerm::from_variable_and_coeff(b, negative_one);
            let mut term = MainGateTerm::new();
            term.add_assign(two_a);
            term.add_assign(minus_b);

            cs.allocate_main_gate(term)?;

            // c - a*b == 0 

            let mut ab_term = ArithmeticTerm::from_variable(a).mul_by_variable(b);
            ab_term.scale(&negative_one);
            let c_term = ArithmeticTerm::from_variable(c);
            let mut term = MainGateTerm::new();
            term.add_assign(c_term);
            term.add_assign(ab_term);

            cs.allocate_main_gate(term)?;

            let dummy = CS::get_dummy_variable();

            // and table (gate #2 zero enumerated)
            {
                let table = cs.get_table(&and_table_name)?;
                let num_keys_and_values = table.width();

                let and_result_value = table.query(&[binary_x_value, binary_y_value])?[0];

                let binary_z = cs.alloc(|| {
                    Ok(and_result_value)
                })?;

                cs.begin_gates_batch_for_step()?;

                let vars = [binary_x, binary_y, binary_z, dummy];
                cs.allocate_variables_without_gate(
                    &vars,
                    &[]
                )?;

                cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

                cs.end_gates_batch_for_step()?;
            }

            // d - 100 == 0 

            let hundred = ArithmeticTerm::constant(E::Fr::from_str("100").unwrap());
            let d_term = ArithmeticTerm::from_variable(d);
            let mut term = MainGateTerm::new();
            term.add_assign(d_term);
            term.sub_assign(hundred);

            cs.allocate_main_gate(term)?;

            let var_zero = cs.get_explicit_zero()?;

            // range table (gate #4 zero enumerated)
            {
                let table = cs.get_table(&range_table_name)?;
                let num_keys_and_values = table.width();

                cs.begin_gates_batch_for_step()?;

                let mut term = MainGateTerm::<E>::new();
                term.add_assign(ArithmeticTerm::from_variable_and_coeff(e, E::Fr::zero()));
                term.add_assign(ArithmeticTerm::from_variable_and_coeff(var_zero, E::Fr::zero()));
                term.add_assign(ArithmeticTerm::from_variable_and_coeff(var_zero, E::Fr::zero()));

                let (vars, coeffs) = CS::MainGate::format_linear_term_with_duplicates(term, dummy)?;

                cs.new_gate_in_batch(
                    &CS::MainGate::default(),
                    &coeffs,
                    &vars,
                    &[]
                )?;

                cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

                cs.end_gates_batch_for_step()?;
            }

            // xor table (gate #5 zero enumerated)
            {
                let table = cs.get_table(&xor_table_name)?;
                let num_keys_and_values = table.width();

                let xor_result_value = table.query(&[binary_x_value, binary_y_value])?[0];

                let binary_z = cs.alloc(|| {
                    Ok(xor_result_value)
                })?;

                cs.begin_gates_batch_for_step()?;

                let vars = [binary_x, binary_y, binary_z, dummy];
                cs.allocate_variables_without_gate(
                    &vars,
                    &[]
                )?;

                cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

                cs.end_gates_batch_for_step()?;
            }

            let one = cs.get_explicit_one()?;

            cs.new_single_gate_for_trace_step(
                &TestBitGate::default(),
                &[],
                &[one],
                &[],
            )?;

            Ok(())
        }
    }


    struct TestCircuit4WithLookupsManyGatesSmallTable<E:Engine>{
        _marker: PhantomData<E>
    }

    impl<E: Engine> Circuit<E> for TestCircuit4WithLookupsManyGatesSmallTable<E> {
        type MainGate = Width4MainGateWithDNext;

        fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
            Ok(
                vec![
                    Width4MainGateWithDNext::default().into_internal(),
                ]
            )
        }

        fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
            let columns = vec![PolyIdentifier::VariablesPolynomial(0), PolyIdentifier::VariablesPolynomial(1), PolyIdentifier::VariablesPolynomial(2)];
            let range_table = LookupTableApplication::new_range_table_of_width_3(2, columns.clone())?;
            let range_table_name = range_table.functional_name();

            let xor_table = LookupTableApplication::new_xor_table(2, columns.clone())?;
            let xor_table_name = xor_table.functional_name();

            let and_table = LookupTableApplication::new_and_table(2, columns)?;
            let and_table_name = and_table.functional_name();

            cs.add_table(range_table)?;
            cs.add_table(xor_table)?;
            cs.add_table(and_table)?;

            let a = cs.alloc(|| {
                Ok(E::Fr::from_str("10").unwrap())
            })?;

            let b = cs.alloc(|| {
                Ok(E::Fr::from_str("20").unwrap())
            })?;

            let c = cs.alloc(|| {
                Ok(E::Fr::from_str("200").unwrap())
            })?;

            let binary_x_value = E::Fr::from_str("3").unwrap();
            let binary_y_value = E::Fr::from_str("1").unwrap();

            let binary_x = cs.alloc(|| {
                Ok(binary_x_value)
            })?;

            let binary_y = cs.alloc(|| {
                Ok(binary_y_value)
            })?;

            let mut negative_one = E::Fr::one();
            negative_one.negate();

            for _ in 0..((1 << 11) - 100) {
                // c - a*b == 0 

                let mut ab_term = ArithmeticTerm::from_variable(a).mul_by_variable(b);
                ab_term.scale(&negative_one);
                let c_term = ArithmeticTerm::from_variable(c);
                let mut term = MainGateTerm::new();
                term.add_assign(c_term);
                term.add_assign(ab_term);

                cs.allocate_main_gate(term)?;
            }

            let dummy = CS::get_dummy_variable();

            // and table
            {
                let table = cs.get_table(&and_table_name)?;
                let num_keys_and_values = table.width();

                let and_result_value = table.query(&[binary_x_value, binary_y_value])?[0];

                let binary_z = cs.alloc(|| {
                    Ok(and_result_value)
                })?;

                cs.begin_gates_batch_for_step()?;

                let vars = [binary_x, binary_y, binary_z, dummy];
                cs.allocate_variables_without_gate(
                    &vars,
                    &[]
                )?;

                cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

                cs.end_gates_batch_for_step()?;
            }

            let var_zero = cs.get_explicit_zero()?;

            // range table
            {
                let table = cs.get_table(&range_table_name)?;
                let num_keys_and_values = table.width();

                cs.begin_gates_batch_for_step()?;

                let mut term = MainGateTerm::<E>::new();
                term.add_assign(ArithmeticTerm::from_variable_and_coeff(binary_y, E::Fr::zero()));
                term.add_assign(ArithmeticTerm::from_variable_and_coeff(var_zero, E::Fr::zero()));
                term.add_assign(ArithmeticTerm::from_variable_and_coeff(var_zero, E::Fr::zero()));

                let (vars, coeffs) = CS::MainGate::format_linear_term_with_duplicates(term, dummy)?;

                cs.new_gate_in_batch(
                    &CS::MainGate::default(),
                    &coeffs,
                    &vars,
                    &[]
                )?;

                cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

                cs.end_gates_batch_for_step()?;
            }

            // xor table
            {
                let table = cs.get_table(&xor_table_name)?;
                let num_keys_and_values = table.width();

                let xor_result_value = table.query(&[binary_x_value, binary_y_value])?[0];

                let binary_z = cs.alloc(|| {
                    Ok(xor_result_value)
                })?;

                cs.begin_gates_batch_for_step()?;

                let vars = [binary_x, binary_y, binary_z, dummy];
                cs.allocate_variables_without_gate(
                    &vars,
                    &[]
                )?;

                cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

                cs.end_gates_batch_for_step()?;
            }

            Ok(())
        }
    }

    #[test]
    fn test_trivial_circuit_with_gate_agnostic_cs() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;

        let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        assert!(assembly.gates.len() == 1);

        // println!("Assembly state polys = {:?}", assembly.storage.state_map);

        // println!("Assembly setup polys = {:?}", assembly.storage.setup_map);    

        println!("Assembly contains {} gates", assembly.n());
        assembly.finalize();
        assert!(assembly.is_satisfied());

        assembly.finalize();

        let worker = Worker::new();

        let (_storage, _permutation_polys) = assembly.perform_setup(&worker).unwrap();
    }

    #[test]
    fn test_circuit_with_combined_lookups() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;

        let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4WithLookups::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        println!("Assembly contains {} gates", assembly.n());
        assert!(assembly.is_satisfied());

        assert!(assembly.gates.len() == 2);

        assembly.finalize();

        println!("Finalized assembly contains {} gates", assembly.n());

        let worker = Worker::new();

        let _ = assembly.prover_stub(&worker).unwrap();
    }

    #[test]
    fn test_circuit_with_small_tables() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;

        let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4WithLookupsManyGatesSmallTable::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        println!("Assembly contains {} gates", assembly.n());
        assert!(assembly.is_satisfied());

        assert_eq!(assembly.gates.len(), 1);
        assert_eq!(assembly.total_length_of_all_tables, 4*4 + 4*4 + 4);

        assembly.finalize();

        println!("Finalized assembly contains {} gates", assembly.n());

        let worker = Worker::new();

        let _ = assembly.prover_stub(&worker).unwrap();
    }

    #[test]
    fn test_circuit_setup_and_prove() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;
        use crate::plonk::better_better_cs::verifier::*;
        use crate::plonk::better_better_cs::setup::VerificationKey;

        let mut assembly
            = SetupAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        println!("Assembly contains {} gates", assembly.n());
        assert!(assembly.is_satisfied());

        assembly.finalize();

        println!("Finalized assembly contains {} gates", assembly.n());

        let worker = Worker::new();

        let setup = assembly.create_setup::<TestCircuit4<Bn256>>(&worker).unwrap();

        let mut assembly = ProvingAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();
        circuit.synthesize(&mut assembly).expect("must work");
        assembly.finalize();

        let size = assembly.n().next_power_of_two();

        use crate::plonk::commitments::transcript::keccak_transcript::RollingKeccakTranscript;
        use crate::kate_commitment::*;

        let crs_mons = Crs::<Bn256, CrsForMonomialForm>::crs_42(size, &worker);

        let proof = assembly.create_proof::<TestCircuit4<Bn256>, RollingKeccakTranscript<Fr>>(
            &worker, 
            &setup, 
            &crs_mons, 
            None
        ).unwrap();

        let vk = VerificationKey::from_setup(&setup, &worker, &crs_mons).unwrap();

        let valid = verify::<Bn256, TestCircuit4<Bn256>, RollingKeccakTranscript<Fr>>(
            &vk,
            &proof,
            None,
        ).unwrap();

        dbg!(valid);
    }

    #[derive(Clone, Debug, Hash, Default)]
    pub struct TestBitGate;

    impl<E: Engine> GateInternal<E> for TestBitGate {
        fn name(&self) -> &'static str {
            "Test bit gate on A"
        }

        fn degree(&self) -> usize {
            2
        }

        fn can_include_public_inputs(&self) -> bool {
            false
        }

        fn all_queried_polynomials(&self) -> Vec<PolynomialInConstraint> {
            vec![
                PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            ]
        }

        fn setup_polynomials(&self) -> Vec<PolyIdentifier> {
            vec![]
        }

        fn variable_polynomials(&self) -> Vec<PolyIdentifier> {
            vec![
                PolyIdentifier::VariablesPolynomial(0),
            ]
        }

        fn benefits_from_linearization(&self) -> bool {
            false
        }

        fn linearizes_over(&self) -> Vec<PolynomialInConstraint> {
            vec![]
        }

        fn needs_opened_for_linearization(&self) -> Vec<PolynomialInConstraint> {
            vec![]
        }

        fn num_quotient_terms(&self) -> usize {
            1
        }

        fn verify_on_row<'a>(&self, row: usize, poly_storage: &AssembledPolynomialStorage<'a, E>, _last_row: bool) -> E::Fr {
            let q_a = poly_storage.get_poly_at_step(PolyIdentifier::VariablesPolynomial(0), row);
            
            // (A - 1) * A 
            let mut tmp = q_a;
            tmp.sub_assign(&E::Fr::one());
            tmp.mul_assign(&q_a);

            tmp
        }

        fn contribute_into_quotient<'a, 'b>(
            &self, 
            domain_size: usize,
            poly_storage: &mut AssembledPolynomialStorage<'a, E>,
            monomials_storage: & AssembledPolynomialStorageForMonomialForms<'b, E>,
            challenges: &[E::Fr],
            omegas_bitreversed: &BitReversedOmegas<E::Fr>,
            _omegas_inv_bitreversed: &OmegasInvBitreversed<E::Fr>,
            worker: &Worker
        ) -> Result<Polynomial<E::Fr, Values>, SynthesisError> {
            assert!(domain_size.is_power_of_two());
            assert_eq!(challenges.len(), <Self as GateInternal<E>>::num_quotient_terms(&self));

            let lde_factor = poly_storage.lde_factor;
            assert!(lde_factor.is_power_of_two());

            assert!(poly_storage.is_bitreversed);

            let coset_factor = E::Fr::multiplicative_generator();
           
            for p in <Self as GateInternal<E>>::all_queried_polynomials(&self).into_iter() {
                ensure_in_map_or_create(&worker, 
                    p, 
                    domain_size, 
                    omegas_bitreversed, 
                    lde_factor, 
                    coset_factor, 
                    monomials_storage, 
                    poly_storage
                )?;
            }

            let ldes_storage = &*poly_storage;

            // (A - 1) * A 
            let a_ref = get_from_map_unchecked(
                PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
                ldes_storage
            );

            let mut tmp = a_ref.clone();
            drop(a_ref);

            let one = E::Fr::one();

            tmp.map(&worker,
                |el| {
                    let mut tmp = *el;
                    tmp.sub_assign(&one);
                    tmp.mul_assign(&*el);

                    *el = tmp;
                }, 
            );

            tmp.scale(&worker, challenges[0]);

            Ok(tmp)
        }

        fn contribute_into_linearization<'a>(
            &self, 
            _domain_size: usize,
            _at: E::Fr,
            _queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
            _monomials_storage: & AssembledPolynomialStorageForMonomialForms<'a, E>,
            _challenges: &[E::Fr],
            _worker: &Worker
        ) -> Result<Polynomial<E::Fr, Coefficients>, SynthesisError> {
            unreachable!("this gate does not contribute into linearization");
        }
        fn contribute_into_verification_equation(
            &self, 
            _domain_size: usize,
            _at: E::Fr,
            queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
            challenges: &[E::Fr],
        ) -> Result<E::Fr, SynthesisError> {
            assert_eq!(challenges.len(), 1);
            // (A-1) * A 
            let a_value = *queried_values.get(&PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)))
                .ok_or(SynthesisError::AssignmentMissing)?;
            let mut result = a_value;
            result.sub_assign(&E::Fr::one());
            result.mul_assign(&a_value);
            result.mul_assign(&challenges[0]);
            
            Ok(result)
        }

        fn put_public_inputs_into_selector_id(&self) -> Option<usize> {
            None
        }

        fn box_clone(&self) -> Box<dyn GateInternal<E>> {
            Box::from(self.clone())
        }
        fn contribute_into_linearization_commitment(
            &self, 
            _domain_size: usize,
            _at: E::Fr,
            _queried_values: &std::collections::HashMap<PolynomialInConstraint, E::Fr>,
            _commitments_storage: &std::collections::HashMap<PolyIdentifier, E::G1Affine>,
            _challenges: &[E::Fr],
        ) -> Result<E::G1, SynthesisError> {
            unreachable!("this gate does not contribute into linearization");
        }
    }

    impl<E: Engine> Gate<E> for TestBitGate {}
}