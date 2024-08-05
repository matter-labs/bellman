use crate::pairing::ff::{Field, PrimeField, PrimeFieldRepr};
use crate::pairing::{Engine, CurveAffine, CurveProjective};
use crate::bit_vec::BitVec;

use crate::{SynthesisError};
#[cfg(feature = "allocator")]
use std::alloc::{Allocator, Global};
use std::collections::HashMap;
use std::marker::PhantomData;

use crate::worker::Worker;
use crate::plonk::domains::*;
use crate::plonk::polynomials::*;

pub use crate::plonk::cs::variable::*;
use crate::plonk::better_cs::utils::*;
pub use super::lookup_tables::{*};

use crate::plonk::fft::cooley_tukey_ntt::*;

use super::utils::*;
pub use super::data_structures::*;
pub use super::setup::*;

pub use super::gates::main_gate_with_d_next::*;

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
    fn synthesize<CS: ConstraintSystem<E> + 'static>(&self, cs: &mut CS) -> Result<(), SynthesisError>;
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
    fn all_queried_polynomials(&self) -> &'static [PolynomialInConstraint];
    fn setup_polynomials(&self) -> &'static [PolyIdentifier];
    fn variable_polynomials(&self) -> &'static [PolyIdentifier];
    #[inline]
    fn witness_polynomials(&self) -> &'static [PolyIdentifier] {
        &[]
    }
    fn benefits_from_linearization(&self) -> bool;
    fn linearizes_over(&self) -> &'static [PolynomialInConstraint];
    fn needs_opened_for_linearization(&self) -> &'static [PolynomialInConstraint];
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

use serde::{Serialize, Deserialize};
use smallvec::SmallVec;

pub const DEFAULT_SMALLVEC_CAPACITY: usize = 8;

pub trait MainGate<E: Engine>: Gate<E> {
    const NUM_LINEAR_TERMS: usize;
    const NUM_VARIABLES: usize;
    const NUM_VARIABLES_ON_NEXT_STEP: usize;

    fn range_of_multiplicative_term() -> std::ops::Range<usize>;
    fn range_of_linear_terms() -> std::ops::Range<usize>;
    fn index_for_constant_term() -> usize;
    fn range_of_next_step_linear_terms() -> std::ops::Range<usize>;
    fn format_term(instance: MainGateTerm<E>, padding: Variable) -> Result<(SmallVec<[Variable; DEFAULT_SMALLVEC_CAPACITY]>, SmallVec<[E::Fr; DEFAULT_SMALLVEC_CAPACITY]>), SynthesisError>;
    fn format_linear_term_with_duplicates(instance: MainGateTerm<E>, padding: Variable) -> Result<(SmallVec<[Variable; DEFAULT_SMALLVEC_CAPACITY]>, SmallVec<[E::Fr; DEFAULT_SMALLVEC_CAPACITY]>), SynthesisError>;
    fn dummy_vars_to_inscribe(dummy: Variable) -> SmallVec<[Variable; DEFAULT_SMALLVEC_CAPACITY]>;
    fn empty_coefficients() -> SmallVec<[E::Fr; DEFAULT_SMALLVEC_CAPACITY]>;
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
    Product(smallvec::SmallVec<[Variable; 2]>, E::Fr),
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
                let terms = smallvec::smallvec![this, other];

                ArithmeticTerm::Product(terms, coeff)
            },
            ArithmeticTerm::Constant(coeff) => {
                let terms = smallvec::smallvec![other];

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

const DEFAULT_SMALLVEC_CAPACITY_FOR_TERM: usize = 8;

#[derive(Clone, Debug)]
pub struct MainGateTerm<E: Engine>{
    pub(crate) terms: smallvec::SmallVec<[ArithmeticTerm<E>; DEFAULT_SMALLVEC_CAPACITY_FOR_TERM]>,
    pub(crate) vars_scratch: std::collections::HashMap<Variable, usize>,
    pub(crate) num_multiplicative_terms: usize,
    pub(crate) num_constant_terms: usize
}

impl<E: Engine> MainGateTerm<E> {
    pub fn new() -> Self {
        Self {
            terms: smallvec::smallvec![],
            vars_scratch: std::collections::HashMap::with_capacity(DEFAULT_SMALLVEC_CAPACITY_FOR_TERM),
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

pub(crate) struct SimpleBitmap(u64, usize);

impl SimpleBitmap {
    pub(crate) fn new() -> Self {
        Self(0u64, 0)
    }

    pub(crate) fn get_next_unused(&mut self) -> usize {
        for i in 0..64 {
            if self.get(i) == false {
                return i;
            }
        }

        unreachable!()
    }

    pub(crate) fn get(&self, idx: usize) -> bool {
        1u64 << idx & self.0 > 0
    }

    pub(crate) fn set(&mut self, idx: usize) {
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

#[cfg(not(feature="allocator"))]
macro_rules! new_vec_with_allocator {
    ($capacity:expr) => {
        Vec::with_capacity($capacity)
    }
}

#[cfg(feature="allocator")]
macro_rules! new_vec_with_allocator {
    ($capacity:expr) => {
        Vec::with_capacity_in($capacity, A::default())
    }
}


use crate::plonk::polynomials::*;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[cfg_attr(feature = "allocator",serde(bound(serialize = "A: serde::Serialize", deserialize = "'de: 'static, A: serde::Deserialize<'de>")))]
#[cfg_attr(not(feature = "allocator"), serde(bound(deserialize = "'de: 'static")))]
pub struct PolynomialStorage<E: Engine, #[cfg(feature = "allocator")]A: Allocator + Default = Global> {
    #[cfg(feature = "allocator")]
    #[cfg_attr(feature = "allocator", serde(serialize_with = "serialize_hashmap_with_allocator"))]
    #[cfg_attr(feature = "allocator", serde(deserialize_with = "deserialize_hashmap_with_allocator"))]
    pub state_map: std::collections::HashMap<PolyIdentifier, Vec<Variable, A>>,
    #[cfg(not(feature = "allocator"))]
    pub state_map: std::collections::HashMap<PolyIdentifier, Vec<Variable>>,
    pub witness_map: std::collections::HashMap<PolyIdentifier, Vec<E::Fr>>,
    #[cfg(feature = "allocator")]
    #[cfg_attr(feature = "allocator", serde(serialize_with = "serialize_hashmap_with_allocator"))]
    #[cfg_attr(feature = "allocator", serde(deserialize_with = "deserialize_hashmap_with_allocator"))]
    pub setup_map: std::collections::HashMap<PolyIdentifier, Vec<E::Fr, A>>,
    #[cfg(not(feature = "allocator"))]
    pub setup_map: std::collections::HashMap<PolyIdentifier, Vec<E::Fr>>,
}

macro_rules! impl_poly_storage {
    (impl PolynomialStorage $inherent:tt) => {
        #[cfg(feature = "allocator")]
        impl<E: Engine, A: Allocator + Default> PolynomialStorage<E, A> $inherent

        #[cfg(not(feature = "allocator"))]
        impl<E: Engine> PolynomialStorage<E> $inherent
    };
}

impl_poly_storage! {
    impl PolynomialStorage {
        pub fn new() -> Self {
            Self {
                state_map: std::collections::HashMap::new(),
                witness_map: std::collections::HashMap::new(),
                setup_map: std::collections::HashMap::new(),
            }
        }


        pub fn new_specialized_for_proving_assembly_and_state_4(size: usize) -> Self {
            assert!(size <= 1 << <E::Fr as PrimeField>::S);
            let mut state_map = std::collections::HashMap::new();
            for idx in 0..4{
                state_map.insert(PolyIdentifier::VariablesPolynomial(idx), new_vec_with_allocator!(size));
            }
            Self {
                state_map,
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
}
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(bound(serialize = "dyn GateInternal<E>: serde::Serialize", deserialize = "'de: 'static, dyn GateInternal<E>: serde::Deserialize<'de>"))]
pub struct GateDensityStorage<E: Engine>(pub std::collections::HashMap<Box<dyn GateInternal<E>>, BitVec>);

impl<E: Engine> Default for GateDensityStorage<E>{
    fn default() -> Self {
        Self(std::collections::HashMap::new())
    }
}

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

#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[cfg_attr(feature = "allocator", serde(bound(serialize = "MG: serde::Serialize, A: serde::Serialize", deserialize = "'de: 'static, MG: serde::Deserialize<'de>, A: serde::Deserialize<'de>")))]
#[cfg_attr(not(feature = "allocator"), serde(bound(serialize = "MG: serde::Serialize", deserialize = "'de: 'static, MG: serde::Deserialize<'de>")))]
pub struct Assembly<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode, #[cfg(feature = "allocator")]A: Allocator + Default = Global> {
    #[cfg(feature = "allocator")]
    pub inputs_storage: PolynomialStorage<E, A>,
    #[cfg(not(feature = "allocator"))]
    pub inputs_storage: PolynomialStorage<E>,
    #[cfg(feature = "allocator")]
    pub aux_storage: PolynomialStorage<E, A>,
    #[cfg(not(feature = "allocator"))]
    pub aux_storage: PolynomialStorage<E>,
    pub num_input_gates: usize,
    pub num_aux_gates: usize,
    pub max_constraint_degree: usize,
    pub main_gate: MG,
    pub input_assingments: Vec<E::Fr>,
    #[cfg(feature = "allocator")]
    #[cfg_attr(feature = "allocator", serde(serialize_with = "serialize_vec_with_allocator"))]
    #[cfg_attr(feature = "allocator", serde(deserialize_with = "deserialize_vec_with_allocator"))]
    pub aux_assingments: Vec<E::Fr, A>,
    #[cfg(not(feature = "allocator"))]
    pub aux_assingments: Vec<E::Fr>,
    pub num_inputs: usize,
    pub num_aux: usize,
    pub trace_step_for_batch: Option<usize>,
    pub is_finalized: bool,
    #[serde(skip)]
    pub gates: std::collections::HashSet<Box<dyn GateInternal<E>>>,
    pub all_queried_polys_in_constraints: std::collections::HashSet<PolynomialInConstraint>,
    // pub sorted_setup_polynomial_ids: Vec<PolyIdentifier>,
    #[serde(skip)]
    pub sorted_gates: Vec<Box<dyn GateInternal<E>>>,
    #[serde(skip)]
    pub aux_gate_density: GateDensityStorage<E>,
    pub explicit_zero_variable: Option<Variable>,
    pub explicit_one_variable: Option<Variable>,
    #[serde(skip)]
    pub tables: Vec<Arc<LookupTableApplication<E>>>,
    #[serde(skip)]
    pub multitables: Vec<Arc<MultiTableApplication<E>>>,
    pub table_selectors: std::collections::HashMap<String, BitVec>,
    pub multitable_selectors: std::collections::HashMap<String, BitVec>,
    pub table_ids_poly: Vec<E::Fr>,
    pub total_length_of_all_tables: usize,

    pub individual_table_canonical_sorted_entries: std::collections::HashMap<String, Vec<[E::Fr; 3]>>,
    pub individual_table_entries_lookups: std::collections::HashMap<String, std::collections::HashMap<[E::Fr; 3], usize>>,
    #[cfg(feature = "allocator")]
    #[cfg_attr(feature = "allocator", serde(serialize_with = "serialize_hashmap_with_allocator"))]
    #[cfg_attr(feature = "allocator", serde(deserialize_with = "deserialize_hashmap_with_allocator"))]
    pub individual_table_entries: std::collections::HashMap<String, Vec<u32, A>>,
    #[cfg(not(feature = "allocator"))]
    pub individual_table_entries: std::collections::HashMap<String, Vec<u32>>,
    #[cfg(feature = "allocator")]
    #[cfg_attr(feature = "allocator", serde(serialize_with = "serialize_2d_vec_with_allocator"))]
    #[cfg_attr(feature = "allocator", serde(deserialize_with = "deserialize_2d_vec_with_allocator"))]
    pub reusable_buffer_for_lookup_entries: Vec<Vec<u32, A>>,
    #[cfg(not(feature = "allocator"))]
    pub reusable_buffer_for_lookup_entries: Vec<Vec<u32>>,
    pub individual_multitable_entries: std::collections::HashMap<String, Vec<Vec<E::Fr>>>,
    pub known_table_ids: HashMap<String, E::Fr>,
    pub known_table_names: Vec<String>,
    pub num_table_lookups: usize,
    pub num_multitable_lookups: usize,

    _marker_p: std::marker::PhantomData<P>,
    _marker_s: std::marker::PhantomData<S>,
    #[cfg(feature = "allocator")]
    _marker_a: std::marker::PhantomData<A>,
}

cfg_if!{
    if #[cfg(feature = "allocator")]{
        use serde::de::{Visitor, SeqAccess, MapAccess};

        struct VecVisitor<T, B: Allocator> {
            m1: PhantomData<T>,
            m2: PhantomData<B>,
        }
        
        impl<T, B: Allocator> VecVisitor<T, B>{
            pub fn new() -> Self{
                Self{
                    m1: PhantomData,
                    m2: PhantomData,
                }
            }
        }
        
        impl<'de, T, B> Visitor<'de> for VecVisitor<T, B>
        where
            T: Deserialize<'de>,
            B: Allocator + Default,
        {
            type Value = Vec<T, B>;
        
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("a sequence")
            }
        
            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let size_hint = seq.size_hint();
                let size_hint = std::cmp::min(size_hint.unwrap_or(0), 4096);
                let mut values = Vec::with_capacity_in(size_hint, B::default());
        
                while let Ok(result) = seq.next_element() {
                    match result{
                        Some(value) => values.push(value),
                        None => (),
                    }
                }
        
        
                Ok(values)
            }
        }
        
        struct TwoDVecVisitor<T, B: Allocator> {
            m1: PhantomData<T>,
            m2: PhantomData<B>,
        }
        
        
        impl<T, B: Allocator> TwoDVecVisitor<T, B>{
            pub fn new() -> Self{
                Self{
                    m1: PhantomData,
                    m2: PhantomData,
                }
            }
        }
        
        impl<'de, T, B> Visitor<'de> for TwoDVecVisitor<T, B>
        where
            T: Deserialize<'de>,
            B: Allocator + Default,
        {
            type Value = Vec<Vec<T, B>>;
        
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("a sequence")
            }
        
            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let size_hint = seq.size_hint();
                let size_hint = std::cmp::min(size_hint.unwrap_or(0), 4096);
                let mut final_result = Vec::with_capacity(size_hint);
                
        
                while let Ok(result)  = seq.next_element::<Vec<T>>() {
                    match result{                
                        Some(sub_vec) => {
                            let size_hint = seq.size_hint();
                            let size_hint = std::cmp::min(size_hint.unwrap_or(0), 4096);
                            let mut values = Vec::with_capacity_in(size_hint, B::default());
                            for el in sub_vec{
                                values.push(el)
                            }
                            final_result.push(values);                    
                    },
                        None => (),
                    }
                }
        
                Ok(final_result)
            }
        }
        
        struct MapVisitor<K, T, B: Allocator> {
            m0: PhantomData<K>,
            m1: PhantomData<T>,
            m2: PhantomData<B>,
        }
        
        impl<K, T, B: Allocator> MapVisitor<K, T, B>{
            pub fn new() -> Self{
                Self{
                    m0: PhantomData,
                    m1: PhantomData,
                    m2: PhantomData,
                }
            }
        }
        
        impl<'de, K, T, B> Visitor<'de> for MapVisitor<K, T, B>
        where
            T: Deserialize<'de>,
            K: Deserialize<'de>,
            B: Allocator + Default
        {
            type Value = HashMap<K, Vec<T, B>>;
        
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("a 2d sequence")
            }
        
            fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
            where
                A: MapAccess<'de>,
            {
                let size_hint = map.size_hint();
                let size_hint = std::cmp::min(size_hint.unwrap_or(0), 4096);
                let mut final_map = HashMap::with_capacity(size_hint);
        
                while let Ok(entry) = map.next_entry::<K, Vec<T>>() {
                    let mut values = vec![];
                    match entry{
                        Some((key, sub_vec)) => {
                            for el in sub_vec{
                                values.push(el)
                            }
                        },
                        None => (),
                    }
                }
        
        
                Ok(final_map)
            }
        }
        
        fn serialize_vec_with_allocator<T: serde::Serialize, S, A: Allocator + serde::Serialize>(data: &Vec<T, A>, serializer: S) -> Result<S::Ok, S::Error> where S: serde::Serializer {
            data.serialize(serializer)
        }
        
        fn deserialize_vec_with_allocator<'de, D, T: serde::Deserialize<'de>, A: Allocator + Default + serde::Deserialize<'de>>(deserializer: D) -> Result<Vec<T, A>, D::Error> where D: serde::Deserializer<'de> {
            deserializer.deserialize_seq(VecVisitor::new())
        }
        
        fn serialize_2d_vec_with_allocator<T: serde::Serialize, S, A: Allocator + serde::Serialize>(data: &Vec<Vec<T, A>>, serializer: S) -> Result<S::Ok, S::Error> where S: serde::Serializer {
            use serde::ser::SerializeSeq;
            let mut seq = serializer.serialize_seq(Some(data.len()))?;
            for sub_vec in data {
                for el in sub_vec{
                    seq.serialize_element(el)?;
                }
            }
            seq.end()
        }
        
        fn deserialize_2d_vec_with_allocator<'de, D, T: serde::Deserialize<'de>, A: Allocator + Default + serde::Deserialize<'de>>(deserializer: D) -> Result<Vec<Vec<T, A>>, D::Error> where D: serde::Deserializer<'de> {    
            deserializer.deserialize_seq(TwoDVecVisitor::new())
        }
        
        fn serialize_hashmap_with_allocator<K: serde::Serialize, T: serde::Serialize, S, A: Allocator + Default + serde::Serialize>(data: &HashMap<K, Vec<T, A>>, serializer: S) -> Result<S::Ok, S::Error> where S: serde::Serializer {
            use serde::ser::{SerializeMap, SerializeSeq};
        
            let mut s = serializer.serialize_map(Some(data.len()))?;
            for (k, v) in data{
                s.serialize_key(k)?;
                for el in v.iter(){
                    s.serialize_value(el)?;
                }
            }
            s.end()
        }
        
        fn deserialize_hashmap_with_allocator<'de, D,K: serde::Deserialize<'de>, T: serde::Deserialize<'de>, A: Allocator + Default>(deserializer: D) -> Result<HashMap<K, Vec<T, A>>, D::Error> where D: serde::Deserializer<'de> {
            deserializer.deserialize_map(MapVisitor::new())
        }
    }
}


macro_rules! impl_assembly {
    {impl Assembly $inherent:tt} => {
        #[cfg(feature = "allocator")]
        impl<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode, A: Allocator + Default + 'static + Send + Sync> Assembly<E, P, MG, S, A> $inherent
        
        #[cfg(not(feature = "allocator"))]
        impl<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode> Assembly<E, P, MG, S> $inherent
    };
    {impl ConstraintSystem $inherent:tt} =>{
        #[cfg(feature = "allocator")]
        impl<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode, A: Allocator + Default + 'static + Send + Sync> ConstraintSystem<E> for Assembly<E, P, MG, S, A>  $inherent

        #[cfg(not(feature = "allocator"))]
        impl<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode> ConstraintSystem<E> for Assembly<E, P, MG, S> $inherent
    }
}

impl_assembly!{
    impl ConstraintSystem {
        type Params = P;
        type MainGate = MG;
    
        // allocate a variable
        #[inline]
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
    
        #[inline]
        fn get_main_gate(&self) -> &MG {
            &self.main_gate
        }
    
        #[inline]
        fn begin_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
            debug_assert!(self.trace_step_for_batch.is_none());
            let n = self.num_aux_gates;
            self.num_aux_gates += 1;
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
            debug_assert_eq!(n+1, self.num_aux_gates, "invalid batch id");
    
            Ok(())
        }
    
        #[inline]
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
    
        #[inline]
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
            assert!(!self.known_table_ids.contains_key(&table.functional_name()), "can not add a duplicate table for name {}", table.functional_name());
            let table_name = table.functional_name();
            let table_id = table.table_id();
            let number_of_entries = table.size();
    
            // ensure sorted format when we add table
            let mut entries = Self::ensure_sorted_table(&table);
            assert_eq!(entries.len(), 3);
            let mut entries_as_arrays = Vec::with_capacity(entries[0].len());
            let mut entries_into_table_row = std::collections::HashMap::with_capacity(entries[0].len());
            let column_2 = entries.pop().unwrap();
            let column_1 = entries.pop().unwrap();
            let column_0 = entries.pop().unwrap();
            for (idx, ((a, b), c)) in column_0.into_iter().zip(column_1.into_iter()).zip(column_2.into_iter()).enumerate() {
                entries_as_arrays.push([a, b, c]);
                entries_into_table_row.insert([a, b, c], idx);
            }
    
            let shared = Arc::from(table);
            let res = shared.clone();
    
            self.tables.push(shared);
            self.individual_table_canonical_sorted_entries.insert(table_name.clone(), entries_as_arrays);
            self.individual_table_entries_lookups.insert(table_name.clone(), entries_into_table_row);            
            let buffer_for_current_table = if let Some(mut buffer) =  self.reusable_buffer_for_lookup_entries.pop(){
                buffer.clear();

                buffer
            }else{
                // println!("allocating new buffer for table {}", table_name);
                
                new_vec_with_allocator!(0)
            };
            
            self.individual_table_entries.insert(table_name.clone(), buffer_for_current_table);
            self.known_table_names.push(table_name.clone());
            self.table_selectors.insert(table_name.clone(), BitVec::new());
            self.known_table_ids.insert(table_name, table_id);
    
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
            self.individual_table_entries.insert(table_name.clone(), new_vec_with_allocator!(0));
    
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
    
            if S::PRODUCE_SETUP {
                debug_assert!(self.tables.contains(&table));
                assert!(table.can_be_combined() == true);
                assert!(table.applies_over().len() == 3);
    
                let table_name = table.functional_name();
                let table_id = table.table_id();
            
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
                self.table_ids_poly.push(table_id);
            }
    
            if S::PRODUCE_WITNESS {
                let table_name = table.functional_name();
    
                // add values for lookup table sorting later
                let keys_and_values_len = table.applies_over().len();
                let mut table_entries = arrayvec::ArrayVec::<_, 3>::new();
                for v in variables.iter() {
                    let value = self.get_value(*v).unwrap();
                    table_entries.push(value);
                }
                use std::convert::TryInto;
                let table_entries_as_array: [_; 3] = table_entries.into_inner().unwrap();
    
                let entries = self.individual_table_entries.get_mut(&table_name).unwrap();
                assert_eq!(variables.len(), table.applies_over().len());

                // // This check is substituted by the lookup from values into index below
                // let valid_entries = table.is_valid_entry(&table_entries_as_array[..keys_and_values_len]);
                // assert!(valid_entries);
    
                // if !valid_entries {
                //     return Err(SynthesisError::Unsatisfiable);
                // }

                let row_idx = self.individual_table_entries_lookups.get(&table_name).unwrap().get(&table_entries_as_array);
                assert!(row_idx.is_some(), "table most likely doesn't contain a row for {:?}", table_entries_as_array);
    
                entries.push(*row_idx.unwrap() as u32);
            }
    
            self.num_table_lookups += 1;
    
            Ok(())
        }
        
        #[track_caller]
        fn apply_multi_lookup_gate(&mut self, variables: &[Variable], table: Arc<MultiTableApplication<E>>) -> Result<(), SynthesisError> {
            unimplemented!("not implementing multitable for now");
        }
    
        fn get_current_step_number(&self) -> usize {
            self.n()
        }
    
        fn get_current_aux_gate_number(&self) -> usize {
            self.num_aux_gates
        }
    }
    
}

impl_assembly!{
    impl Assembly{
        fn allocate_into_storage<G: Gate<E>>(
            gate: &G,
            #[cfg(feature = "allocator")]
            storage: &mut PolynomialStorage<E, A>, 
            #[cfg(not(feature = "allocator"))]
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

                for &setup_poly in gate.setup_polynomials().into_iter() {
                    let poly_ref = storage.setup_map.entry(setup_poly).or_insert(new_vec_with_allocator!(0));
                    if poly_ref.len() < n {
                        poly_ref.resize(n, E::Fr::zero());
                    }
                    poly_ref.push(*coeffs_it.next().unwrap_or(&zero));
                }

                debug_assert!(coeffs_it.next().is_none(), "must consume all the coefficients for gate");

            }

            let mut variable_it = variables_assignments.iter();

            for &var_poly in gate.variable_polynomials().into_iter() {
                let poly_ref = storage.state_map.entry(var_poly).or_insert(new_vec_with_allocator!(0));
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

            for &key in gate.witness_polynomials().into_iter() {
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
                for &p in gate.all_queried_polynomials().into_iter() {
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
                aux_assingments: new_vec_with_allocator!(0),

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

                individual_table_canonical_sorted_entries: std::collections::HashMap::new(),
                individual_table_entries_lookups: std::collections::HashMap::new(),
                individual_table_entries: std::collections::HashMap::new(),
                reusable_buffer_for_lookup_entries: vec![],
                individual_multitable_entries: std::collections::HashMap::new(),

                known_table_ids: HashMap::new(),
                known_table_names: vec![],

                num_table_lookups: 0,
                num_multitable_lookups: 0,

                _marker_p: std::marker::PhantomData,
                _marker_s: std::marker::PhantomData,
                #[cfg(feature = "allocator")]
                _marker_a: std::marker::PhantomData,
            };

            tmp.add_gate_into_list(&MG::default());

            tmp
        }        

        pub fn new_specialized_for_proving_assembly_and_state_4(domain_size: usize, aux_size: usize, num_lookup_tables: usize, max_num_lookup_entries: usize) -> Self {
            assert!(domain_size <= 1 << <E::Fr as PrimeField>::S);
            let reusable_buffer_for_lookup_entries = (0..num_lookup_tables).map(|_| new_vec_with_allocator!(max_num_lookup_entries)).collect();
            let mut tmp = Self {
                inputs_storage: PolynomialStorage::new(),
                aux_storage: PolynomialStorage::new_specialized_for_proving_assembly_and_state_4(domain_size),

                max_constraint_degree: 0,

                num_input_gates: 0,
                num_aux_gates: 0,

                num_inputs: 0,
                num_aux: 0,

                input_assingments: vec![],
                aux_assingments: new_vec_with_allocator!(aux_size),

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

                individual_table_canonical_sorted_entries: std::collections::HashMap::new(),
                individual_table_entries_lookups: std::collections::HashMap::new(),
                individual_table_entries: std::collections::HashMap::new(),
                reusable_buffer_for_lookup_entries: reusable_buffer_for_lookup_entries,
                individual_multitable_entries: std::collections::HashMap::new(),

                known_table_ids: HashMap::new(),
                known_table_names: vec![],

                num_table_lookups: 0,
                num_multitable_lookups: 0,

                _marker_p: std::marker::PhantomData,
                _marker_s: std::marker::PhantomData,
                #[cfg(feature = "allocator")]
                _marker_a: std::marker::PhantomData,
            };

            tmp.add_gate_into_list(&MG::default());

            tmp
        }


        // return variable that is not in a constraint formally, but has some value
        const fn dummy_variable() -> Variable {
            Variable(Index::Aux(0))
        }

        pub fn finalize(&mut self) {
            if self.is_finalized {
                return;
            }

            // the lookup argument (as in the paper) will make two polynomials to fit jointly sorted set
            // but in practice we fit it into one. For this purpose we only need to add as many empty rows
            // as there is all the tables (as the worst case)

            // let mut min_space_for_lookups = self.num_table_lookups;
            // for table in self.tables.iter() {
            //     let table_num_rows = table.size();
            //     min_space_for_lookups += table_num_rows;
            // }
            // min_space_for_lookups += self.n();
            // let new_size_candidates = [(self.n() + 1).next_power_of_two() - 1, (min_space_for_lookups + 1).next_power_of_two() - 1];

            // In better case it's enough for us to have num_lookups + total length of tables to be smaller
            // than problem size, so joinly sorted set fits into 1 polynomial, and we use zeroes as padding values

            let total_number_of_table_entries = self.num_table_lookups + self.total_length_of_all_tables;
            let new_size_candidates = [(self.n() + 1).next_power_of_two() - 1, (total_number_of_table_entries + 1).next_power_of_two() - 1];

            let new_size = *new_size_candidates.iter().max().unwrap();
            assert!(
                new_size <= 1usize << E::Fr::S, 
                "Padded circuit size is {}, that is larget than number of roots of unity 2^{}. Padded from {} gates and {} lookup table accesses", 
                new_size, 
                E::Fr::S,
                self.n(),
                total_number_of_table_entries,
            );
            assert!(
                new_size <= (1usize << E::Fr::S) / <Self as ConstraintSystem<E>>::Params::STATE_WIDTH, 
                "Circuit size is {}, that is larget than number of roots of unity 2^{} for copy-permutation over {} polys. Padded from {} gates and {} lookup table accesses", 
                new_size, 
                E::Fr::S, 
                <Self as ConstraintSystem<E>>::Params::STATE_WIDTH,
                self.n(),
                total_number_of_table_entries,
            );

            let dummy = Self::get_dummy_variable();

            let empty_vars = vec![dummy; <Self as ConstraintSystem<E>>::Params::STATE_WIDTH];
            let empty_witness = vec![E::Fr::zero(); <Self as ConstraintSystem<E>>::Params::WITNESS_WIDTH];

            for _ in self.n()..new_size {
                self.begin_gates_batch_for_step().unwrap();

                self.allocate_variables_without_gate(
                    &empty_vars,
                    &empty_witness
                ).expect("must add padding gate");

                self.end_gates_batch_for_step().unwrap();
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

        // Caller can specify how large circuit should be artificially inflated
        // if possible. Will panic if size is already too large
        pub fn finalize_to_size_log_2(&mut self, size_log_2: usize) {
            if self.is_finalized {
                return;
            }

            assert!(size_log_2 <= E::Fr::S as usize);

            // the lookup argument (as in the paper) will make two polynomials to fit jointly sorted set
            // but in practice we fit it into one. For this purpose we only need to add as many empty rows
            // as there is all the tables (as the worst case)

            // In better case it's enough for us to have num_lookups + total length of tables to be smaller
            // than problem size, so joinly sorted set fits into 1 polynomial, and we use zeroes as padding values

            let total_number_of_table_entries = self.num_table_lookups + self.total_length_of_all_tables;
            let new_size_candidates = [(self.n() + 1).next_power_of_two() - 1, (total_number_of_table_entries + 1).next_power_of_two() - 1];

            let new_size = *new_size_candidates.iter().max().unwrap();
            assert!(
                new_size <= 1usize << E::Fr::S, 
                "
                size is {}, that is larget than number of roots of unity 2^{}. Padded from {} gates and {} lookup table accesses", 
                new_size, 
                E::Fr::S,
                self.n(),
                total_number_of_table_entries,
            );
            assert!(
                new_size <= (1usize << E::Fr::S) / <Self as ConstraintSystem<E>>::Params::STATE_WIDTH, 
                "Circuit size is {}, that is larget than number of roots of unity 2^{} for copy-permutation over {} polys. Padded from {} gates and {} lookup table accesses", 
                new_size, 
                E::Fr::S, 
                <Self as ConstraintSystem<E>>::Params::STATE_WIDTH,
                self.n(),
                total_number_of_table_entries,
            );

            let pad_to = 1 << size_log_2;

            let new_size = if new_size <= pad_to {
                pad_to - 1
            } else {
                panic!("Requested padding to size 2^{}, but circuit already contains {} gates", size_log_2, new_size)
            };

            let dummy = Self::get_dummy_variable();

            let empty_vars = vec![dummy; <Self as ConstraintSystem<E>>::Params::STATE_WIDTH];
            let empty_witness = vec![E::Fr::zero(); <Self as ConstraintSystem<E>>::Params::WITNESS_WIDTH];

            for _ in self.n()..new_size {
                self.begin_gates_batch_for_step().unwrap();

                self.allocate_variables_without_gate(
                    &empty_vars,
                    &empty_witness
                ).expect("must add padding gate");

                self.end_gates_batch_for_step().unwrap();
            }
            assert_eq!(new_size, self.n());

            let new_size_for_aux = new_size - self.num_input_gates;

            if S::PRODUCE_SETUP {
                // pad gate selectors
                for (_, tracker) in self.aux_gate_density.0.iter_mut() {
                    tracker.grow(new_size_for_aux, false);
                }

                if self.num_table_lookups > 0{
                    // pad lookup selectors
                    for (_, selector) in self.table_selectors.iter_mut() {
                        selector.grow(new_size_for_aux, false);
                    }

                    // pad special purpose table selector poly
                    self.table_ids_poly.resize(new_size_for_aux, E::Fr::zero());
                }            
            }

            assert!((self.n()+1).is_power_of_two(), "padded circuit size is not power of two. self.n() = {}", self.n());
            
            self.is_finalized = true;
        }

        #[cfg(feature = "allocator")]
        fn get_storage_for_trace_step(&self, step: usize) -> &PolynomialStorage<E, A> {
            if step < self.num_input_gates {
                &self.inputs_storage
            } else {
                &self.aux_storage
            }
        }

        #[cfg(not(feature = "allocator"))]
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

        pub fn make_permutations(&self, worker: &Worker) -> Result<Vec<Polynomial::<E::Fr, Values>>, SynthesisError> {
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

        pub fn make_setup_polynomials(
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
            if self.sorted_gates.len() == 1 {
                return Ok(vec![]);
            }

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
                let entries = self.individual_table_canonical_sorted_entries.get(&table.functional_name()).unwrap();
                assert!(entries.len() == table.size(), "invalid number of elements in table {}", table.functional_name());
                // these are individual column vectors, so just copy
                for e in entries.iter() {
                    for idx in 0..3 {
                        column_contributions[idx].push(e[idx]);
                    }

                }

                let table_id = table.table_id();
                let pad_to_len = column_contributions[0].len();

                column_contributions.last_mut().unwrap().resize(pad_to_len, table_id);
            }


            Ok(column_contributions)
        }

        pub fn ensure_sorted_table(table: &LookupTableApplication<E>) -> Vec<Vec<E::Fr>> {
            let entries = table.get_table_values_for_polys();
            assert_eq!(entries.len(), 3);

            let mut uniqueness_checker = std::collections::HashSet::with_capacity(entries[0].len());

            // sort them in a standard lexicographic way, so our sorting is always simple
            let size = entries[0].len();
            let mut kv_set_entries = Vec::with_capacity(size);
            for i in 0..size {
                let entry = KeyValueSet::<E>::new([entries[0][i], entries[1][i], entries[2][i]]);
                let is_unique = uniqueness_checker.insert(entry);
                assert!(is_unique);
                kv_set_entries.push(entry);
            }

            kv_set_entries.sort();

            let mut result = vec![Vec::with_capacity(size); 3];

            for kv in kv_set_entries.iter() {
                for i in 0..3 {
                    result[i].push(kv.inner[i]);
                }
            }

            result
        }

        // pub fn calculate_interleaved_t_polys(&self) -> 
        //     Result<Vec<Vec<E::Fr>>, SynthesisError> {
        //     assert!(self.is_finalized);
            
        //     if self.tables.len() == 0 {
        //         return Ok(vec![])
        //     }
            
        //     // we should pass over every table and append it

        //     let mut width = 0;
        //     for table in self.tables.iter() {
        //         if width == 0 {
        //             width = table.width();
        //         } else {
        //             assert_eq!(width, table.width());
        //         }
        //     }

        //     assert_eq!(width, 3, "only support tables that span over 3 polynomials for now");

        //     assert!(self.is_finalized);
        //     if self.tables.len() == 0 {
        //         return Ok(vec![]);
        //     }

        //     // total number of gates, Input + Aux
        //     let size = self.n();

        //     let aux_gates_start = self.num_input_gates;

        //     let mut contributions = vec![vec![E::Fr::zero(); size]; 4];

        //     // make it shifted for ease of rotations later one
        //     let mut place_into_idx = aux_gates_start + 1;

        //     let lookup_selector = self.calculate_lookup_selector_values()?;

        //     for single_application in self.tables.iter() {
        //         // let entries = single_application.get_table_values_for_polys();
        //         println!("Sorting table {}", single_application.functional_name());
        //         let entries = Self::ensure_sorted_table(single_application);
        //         assert_eq!(entries.len(), 3);
        //         let table_id = single_application.table_id();
        //         let num_entries = single_application.size();

        //         assert_eq!(entries[0].len(), num_entries);

        //         for entry_idx in 0..num_entries {
        //             'inner: loop {
        //                 if lookup_selector[place_into_idx].is_zero() {
        //                     // we can place a table row into the poly
        //                     for column in 0..3 {
        //                         contributions[column][place_into_idx] = entries[column][entry_idx];
        //                     }
        //                     contributions[3][place_into_idx] = table_id;
        //                     place_into_idx += 1;

        //                     break 'inner;
        //                 } else {
        //                     // go for a next one 
        //                     place_into_idx += 1;
        //                 }
        //             }
        //         }
        //     }

        //     Ok(contributions)
        // }

        // pub fn calculate_s_poly_contributions_from_witness(&self) ->
        //     Result<Vec<Vec<E::Fr>>, SynthesisError> 
        // {
        //     if self.tables.len() == 0 {
        //         return Ok(vec![]);
        //     }
        //     // we first form a set of all occured witness values,
        //     // then include table entries to the set
        //     // and then sort this set

        //     let mut kv_set_entries = vec![];
        //     let mut contributions_per_column = vec![vec![]; 4];
        //     for (_table_idx, single_application) in self.tables.iter().enumerate() {
        //         // copy all queries from witness
        //         let table_name = single_application.functional_name();
        //         for kv_values in self.individual_table_entries.get(&table_name).unwrap().iter() {
        //             let entry = KeyValueSet::<E>::from_slice(&kv_values[..3]);
        //             kv_set_entries.push(entry);
        //         }

        //         // copy table elements themselves
        //         let entries = Self::ensure_sorted_table(single_application);
        //         // those are full values of polynomials, so we have to virtually transpose

        //         let size = entries[0].len();
        //         for i in 0..size {
        //             let entry = KeyValueSet::new([entries[0][i], entries[1][i], entries[2][i]]);
        //             kv_set_entries.push(entry)
        //         }

        //         kv_set_entries.sort();

        //         // now copy backward with addition of the table id

        //         for kv in kv_set_entries.iter() {
        //             for i in 0..3 {
        //                 contributions_per_column[i].push(kv.inner[i]);
        //             }
        //         }

        //         let table_id = single_application.table_id();
        //         let pad_to_len = contributions_per_column[0].len();
        //         contributions_per_column.last_mut().unwrap().resize(pad_to_len, table_id);

        //         kv_set_entries.truncate(0);
        //     }

        //     Ok(contributions_per_column)
        // }

        pub fn calculate_s_poly_contributions_from_witness(&self, delinearization_challenge: E::Fr) ->
            Result<Vec<E::Fr>, SynthesisError> 
        {
            if self.tables.len() == 0 {
                return Ok(vec![]);
            }

            // we can sort based on indexes only

            let mut kv_set_entries = vec![];
            let mut accumulated_contributions = vec![];
            for (_table_idx, single_application) in self.tables.iter().enumerate() {
                // copy all queries from witness
                let table_name = single_application.functional_name();
                for kv_index in self.individual_table_entries.get(&table_name).unwrap().iter() {
                    kv_set_entries.push(*kv_index);
                }

                let table_size = self.individual_table_canonical_sorted_entries.get(&table_name).unwrap().len();

                // copy table elements themselves
                for i in 0..table_size {
                    kv_set_entries.push(i as u32);
                }

                kv_set_entries.sort();

                // now copy backward with addition of the table id

                let canonical_rows = self.individual_table_canonical_sorted_entries.get(&table_name).unwrap();
                let table_id = single_application.table_id();
                let mut table_id_contribution = table_id;
                table_id_contribution.mul_assign(&delinearization_challenge);

                let mut collapsed_rows = Vec::with_capacity(table_size);
                for row in canonical_rows.iter() {
                    let mut value = table_id_contribution;
                    // horner rule for row[0] + alpha * row[1] + alpha^2 * row[2] + alpha^3 * table_id
                    let mut tmp = table_id_contribution;
                    tmp.add_assign(&row[2]);
                    tmp.mul_assign(&delinearization_challenge);
                    tmp.add_assign(&row[1]);
                    tmp.mul_assign(&delinearization_challenge);
                    tmp.add_assign(&row[0]);

                    collapsed_rows.push(tmp);
                }

                for el in kv_set_entries.iter() {
                    let value = collapsed_rows[*el as usize];
                    accumulated_contributions.push(value);
                }

                kv_set_entries.truncate(0);
            }

            Ok(accumulated_contributions)
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
                let table_id = single_application.table_id();
                let table_name = single_application.functional_name();
                let keys_and_values = single_application.applies_over();
                let selector_bitvec = self.table_selectors.get(&table_name).unwrap();

                // let num_non_empty = selector_bitvec.iter().filter(|x| *x).count();
                // println!("{} lookups for table {}", num_non_empty, table_name);

                assert!(selector_bitvec.len() >= num_aux_gates);
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

            for a in full_assignments.iter() {
                assert_eq!(a.len(), self.num_input_gates + self.num_aux_gates);
            }

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
    fn test_setup_and_prove_custom_gate_and_tables() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;
        use crate::plonk::better_better_cs::verifier::*;
        use crate::plonk::better_better_cs::setup::VerificationKey;

        let mut assembly = SetupAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4WithLookups::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        println!("Assembly contains {} gates", assembly.n());
        assert!(assembly.is_satisfied());

        assembly.finalize();

        println!("Finalized assembly contains {} gates", assembly.n());

        let worker = Worker::new();

        let setup = assembly.create_setup::<TestCircuit4WithLookups<Bn256>>(&worker).unwrap();

        let mut assembly = ProvingAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();
        circuit.synthesize(&mut assembly).expect("must work");
        assembly.finalize();

        let size = assembly.n().next_power_of_two();

        use crate::plonk::commitments::transcript::keccak_transcript::RollingKeccakTranscript;
        use crate::kate_commitment::*;

        let crs_mons = Crs::<Bn256, CrsForMonomialForm>::crs_42(size, &worker);

        let proof = assembly.create_proof::<TestCircuit4WithLookups<Bn256>, RollingKeccakTranscript<Fr>>(
            &worker, 
            &setup, 
            &crs_mons, 
            None
        ).unwrap();

        let vk = VerificationKey::from_setup(&setup, &worker, &crs_mons).unwrap();

        let valid = verify::<Bn256, TestCircuit4WithLookups<Bn256>, RollingKeccakTranscript<Fr>>(
            &vk,
            &proof,
            None,
        ).unwrap();

        assert!(valid);
        println!("Done!");
    }

    #[test]
    fn test_setup_and_prove_single_gate_and_tables() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;
        use crate::plonk::better_better_cs::verifier::*;
        use crate::plonk::better_better_cs::setup::VerificationKey;

        let mut assembly = SetupAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4WithLookupsManyGatesSmallTable::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        println!("Assembly contains {} gates", assembly.n());
        assert!(assembly.is_satisfied());

        assembly.finalize();

        println!("Finalized assembly contains {} gates", assembly.n());

        let worker = Worker::new();

        let setup = assembly.create_setup::<TestCircuit4WithLookupsManyGatesSmallTable<Bn256>>(&worker).unwrap();

        let mut assembly = ProvingAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();
        circuit.synthesize(&mut assembly).expect("must work");
        assembly.finalize();

        let size = assembly.n().next_power_of_two();

        use crate::plonk::commitments::transcript::keccak_transcript::RollingKeccakTranscript;
        use crate::kate_commitment::*;

        let crs_mons = Crs::<Bn256, CrsForMonomialForm>::crs_42(size, &worker);

        let proof = assembly.create_proof::<TestCircuit4WithLookupsManyGatesSmallTable<Bn256>, RollingKeccakTranscript<Fr>>(
            &worker, 
            &setup, 
            &crs_mons, 
            None
        ).unwrap();

        let vk = VerificationKey::from_setup(&setup, &worker, &crs_mons).unwrap();

        let valid = verify::<Bn256, TestCircuit4WithLookupsManyGatesSmallTable<Bn256>, RollingKeccakTranscript<Fr>>(
            &vk,
            &proof,
            None,
        ).unwrap();

        assert!(valid);
    }

    #[test]
    fn test_bench_long_synthesis() {
        use crate::pairing::bn256::{Bn256, Fr};
        use crate::worker::Worker;
        use crate::plonk::better_better_cs::verifier::*;
        use crate::plonk::better_better_cs::setup::VerificationKey;

        let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();

        let circuit = TestCircuit4::<Bn256> {
            _marker: PhantomData
        };

        circuit.synthesize(&mut assembly).expect("must work");

        dbg!(&assembly.n());

        assembly.finalize_to_size_log_2(26);
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

        fn all_queried_polynomials(&self) -> &'static [PolynomialInConstraint] {
            const A: [PolynomialInConstraint; 1] = [
                PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            ];

            &A
        }

        fn setup_polynomials(&self) -> &'static [PolyIdentifier] {
            &[]
        }

        fn variable_polynomials(&self) -> &'static [PolyIdentifier] {
            const A: [PolyIdentifier; 1] = [
                PolyIdentifier::VariablesPolynomial(0),
            ];

            &A
        }

        fn benefits_from_linearization(&self) -> bool {
            false
        }

        fn linearizes_over(&self) -> &'static [PolynomialInConstraint] {
            &[]
        }

        fn needs_opened_for_linearization(&self) -> &'static [PolynomialInConstraint] {
            &[]
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
           
            for &p in <Self as GateInternal<E>>::all_queried_polynomials(&self).into_iter() {
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
