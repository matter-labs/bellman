use super::cs::*;

use crate::pairing::ff::{Field, PrimeField};
use crate::pairing::{Engine, CurveAffine, EncodedPoint};

use crate::plonk::domains::*;
use crate::plonk::polynomials::*;
use crate::worker::Worker;
use crate::SynthesisError;

use crate::kate_commitment::*;

use std::marker::PhantomData;

use super::utils::*;
use super::LDE_FACTOR;

pub struct SetupPolynomials<E: Engine, P: PlonkConstraintSystemParams<E>> {
    pub n: usize,
    pub num_inputs: usize,
    pub selector_polynomials: Vec<Polynomial<E::Fr, Coefficients>>,
    pub next_step_selector_polynomials: Vec<Polynomial<E::Fr, Coefficients>>,
    pub permutation_polynomials: Vec<Polynomial<E::Fr, Coefficients>>,

    pub(crate) _marker: std::marker::PhantomData<P>,
}

use crate::byteorder::BigEndian;
use crate::byteorder::ReadBytesExt;
use crate::byteorder::WriteBytesExt;
use std::io::{Read, Write};


// pub trait EngineDataSerializationRead: Sized {
//     fn read<R: Read>(reader: R) -> std::io::Result<Self>; 
// }

// pub trait EngineDataSerializationWrite<E: Engine>: Sized {
//     fn write<W: Write>(&self, writer: W) -> std::io::Result<()>;
// }

// pub trait EngineDataSerialization<E: Engine>: EngineDataSerializationRead + EngineDataSerializationWrite<E> {}

// impl<E: Engine> EngineDataSerializationRead for E::Fr {
//     fn read<R: Read>(reader: R) -> std::io::Result<Self> {
//         read_fr(reader)
//     }
// }

// impl<E: Engine> EngineDataSerializationWrite<E> for E::Fr {
//     fn write<W: Write>(&self, writer: W) -> std::io::Result<()> {
//         write_fr(self, writer)
//     }
// }

// impl<E: Engine, T> EngineDataSerializationRead for T where E::G1Affine = T {
//     fn read<R: Read>(reader: R) -> std::io::Result<Self> {
//         let mut repr = <Self as CurveAffine>::Uncompressed::empty();
//         reader.read_exact(repr.as_mut())?;

//         let e = repr
//             .into_affine()
//             .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

//         Ok(e)
//     }
// }

// impl<E: Engine> EngineDataSerializationWrite<E> for E::G1Affine {
//     fn write<W: Write>(&self, writer: W) -> std::io::Result<()> {
//         writer.write_all(self.into_uncompressed().as_ref())?;
//     }
// }

pub fn read_curve_affine<G: CurveAffine, R: Read>(mut reader: R) -> std::io::Result<G> {
    let mut repr = G::Uncompressed::empty();
    reader.read_exact(repr.as_mut())?;

    let e = repr
        .into_affine()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

    Ok(e)
}

pub fn read_optional_curve_affine<G: CurveAffine, R: Read>(mut reader: R) -> std::io::Result<Option<G>> {
    use crate::ff::PrimeFieldRepr;

    let is_some = read_optional_flag(&mut reader)?;
    if is_some {
        let el = read_curve_affine(&mut reader)?;

        Ok(Some(el))
    } else {
        Ok(None)
    }
}

pub fn read_curve_affine_vector<G: CurveAffine, R: Read>(mut reader: R) -> std::io::Result<Vec<G>> {
    let num_elements = reader.read_u64::<BigEndian>()?;
    let mut elements = vec![];
    for _ in 0..num_elements {
        let el = read_curve_affine(&mut reader)?;
        elements.push(el);
    }

    Ok(elements)
}

pub fn write_curve_affine<G: CurveAffine, W: Write>(p: &G, mut writer: W) -> std::io::Result<()> {
    writer.write_all(p.into_uncompressed().as_ref())?;

    Ok(())
}

pub fn write_optional_curve_affine<G: CurveAffine, W: Write>(p: &Option<G>, mut writer: W) -> std::io::Result<()> {
    write_optional_flag(p.is_some(), &mut writer)?;
    if let Some(p) = p.as_ref() {
        write_curve_affine(p, &mut writer)?;
    }

    Ok(())
}

pub fn write_curve_affine_vec<G: CurveAffine, W: Write>(p: &[G], mut writer: W) -> std::io::Result<()> {
    writer.write_u64::<BigEndian>(p.len() as u64)?;
    for p in p.iter() {
        write_curve_affine(p, &mut writer)?;
    }
    Ok(())
}

pub fn write_fr<F: PrimeField, W: Write>(el: &F, mut writer: W) -> std::io::Result<()> {
    use crate::ff::PrimeFieldRepr;

    let repr = el.into_repr();
    repr.write_be(&mut writer)?;

    Ok(())
}

pub fn write_optional_fr<F: PrimeField, W: Write>(el: &Option<F>, mut writer: W) -> std::io::Result<()> {
    use crate::ff::PrimeFieldRepr;

    write_optional_flag(el.is_some(), &mut writer)?;
    if let Some(el) = el.as_ref() {
        write_fr(el, &mut writer)?;
    }

    Ok(())
}

pub fn write_fr_vec<F: PrimeField, W: Write>(p: &[F], mut writer: W) -> std::io::Result<()> {
    writer.write_u64::<BigEndian>(p.len() as u64)?;
    for p in p.iter() {
        write_fr(p, &mut writer)?;
    }
    Ok(())
}

pub fn write_fr_raw<F: PrimeField, W: Write>(el: &F, mut writer: W) -> std::io::Result<()> {
    use crate::ff::PrimeFieldRepr;

    let repr = el.into_raw_repr();
    repr.write_be(&mut writer)?;

    Ok(())
}

pub fn read_optional_fr<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Option<F>> {
    use crate::ff::PrimeFieldRepr;

    let is_some = read_optional_flag(&mut reader)?;
    if is_some {
        let el = read_fr(&mut reader)?;

        Ok(Some(el))
    } else {
        Ok(None)
    }
}

pub fn read_fr<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<F> {
    use crate::ff::PrimeFieldRepr;

    let mut repr = F::Repr::default();
    repr.read_be(&mut reader)?;

    F::from_repr(repr).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

pub fn read_fr_vec<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Vec<F>> {
    let num_elements = reader.read_u64::<BigEndian>()?;
    let mut elements = vec![];
    for _ in 0..num_elements {
        let el = read_fr(&mut reader)?;
        elements.push(el);
    }

    Ok(elements)
}

pub fn read_fr_raw<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<F> {
    use crate::ff::PrimeFieldRepr;
    let mut repr = F::Repr::default();
    repr.read_be(&mut reader)?;

    F::from_raw_repr(repr).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

pub fn read_optional_flag<R: Read>(mut reader: R) -> std::io::Result<bool> {
    let value = reader.read_u64::<BigEndian>()?;
    if value == 1 {
        return Ok(true);
    } else if value == 0 {
        return Ok(false);
    }
    
    panic!("invalid encoding of optional flag");
}

pub fn write_optional_flag<W: Write>(is_some: bool, mut writer: W) -> std::io::Result<()> {
    if is_some {
        writer.write_u64::<BigEndian>(1u64)?;
    } else {
        writer.write_u64::<BigEndian>(0u64)?;
    }
    
    Ok(())
}

pub fn read_optional_polynomial_coeffs<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Option<Polynomial<F, Coefficients>>> {
    let is_some = read_optional_flag(&mut reader)?;
    if is_some {
        let p = read_polynomial_coeffs(&mut reader)?;

        Ok(Some(p))
    } else {
        Ok(None)
    }
}

pub fn read_optional_polynomial_values_unpadded<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Option<Polynomial<F, Values>>> {
    let is_some = read_optional_flag(&mut reader)?;
    if is_some {
        let p = read_polynomial_values_unpadded(&mut reader)?;

        Ok(Some(p))
    } else {
        Ok(None)
    }
}

pub fn write_optional_polynomial<F: PrimeField, P: PolynomialForm, W: Write>(p: &Option<Polynomial<F, P>>, mut writer: W) -> std::io::Result<()> {
    write_optional_flag(p.is_some(), &mut writer)?;
    if let Some(p) = p.as_ref() {
        write_polynomial(p, &mut writer)?;
    }
    Ok(())
}

pub fn read_polynomials_coeffs_vec<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Vec<Polynomial<F, Coefficients>>> {
    let num_polys = reader.read_u64::<BigEndian>()?;
    let mut polys = vec![];
    for _ in 0..num_polys {
        let p = read_polynomial_coeffs(&mut reader)?;
        polys.push(p);
    }

    Ok(polys)
}

pub fn read_polynomials_values_unpadded_vec<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Vec<Polynomial<F, Values>>> {
    let num_polys = reader.read_u64::<BigEndian>()?;
    let mut polys = vec![];
    for _ in 0..num_polys {
        let p = read_polynomial_values_unpadded(&mut reader)?;
        polys.push(p);
    }

    Ok(polys)
}

pub fn write_polynomials_vec<F: PrimeField, P: PolynomialForm, W: Write>(p: &[Polynomial<F, P>], mut writer: W) -> std::io::Result<()> {
    writer.write_u64::<BigEndian>(p.len() as u64)?;
    for p in p.iter() {
        write_polynomial(p, &mut writer)?;
    }
    Ok(())
}

pub fn write_polynomial<F: PrimeField, P: PolynomialForm, W: Write>(p: &Polynomial<F, P>, mut writer: W) -> std::io::Result<()> {
    writer.write_u64::<BigEndian>(p.as_ref().len() as u64)?;
    for el in p.as_ref().iter() {
        write_fr(el, &mut writer)?;
    }
    Ok(())
}

pub fn read_polynomial_coeffs<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Polynomial<F, Coefficients>> {
    let num_values = reader.read_u64::<BigEndian>()?;
    let mut poly_coeffs = Vec::with_capacity(num_values as usize);
    for _ in 0..num_values {
        let el = read_fr(&mut reader)?;
        poly_coeffs.push(el);
    }

    Ok(Polynomial::from_coeffs(poly_coeffs).expect("must fit into some domain"))
}

pub fn read_polynomial_values_unpadded<F: PrimeField, R: Read>(mut reader: R) -> std::io::Result<Polynomial<F, Values>> {
    let num_values = reader.read_u64::<BigEndian>()?;
    let mut poly_values = Vec::with_capacity(num_values as usize);
    for _ in 0..num_values {
        let el = read_fr(&mut reader)?;
        poly_values.push(el);
    }

    Ok(Polynomial::from_values_unpadded(poly_values).expect("must fit into some domain"))
}

impl<E: Engine, P: PlonkConstraintSystemParams<E>> SetupPolynomials<E, P> {
    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        writer.write_u64::<BigEndian>(self.n as u64)?;
        writer.write_u64::<BigEndian>(self.num_inputs as u64)?;

        writer.write_u64::<BigEndian>(self.selector_polynomials.len() as u64)?;
        for p in self.selector_polynomials.iter() {
            write_polynomial(p, &mut writer)?;
        }

        writer.write_u64::<BigEndian>(self.next_step_selector_polynomials.len() as u64)?;
        for p in self.next_step_selector_polynomials.iter() {
            write_polynomial(p, &mut writer)?;
        }

        writer.write_u64::<BigEndian>(self.permutation_polynomials.len() as u64)?;
        for p in self.permutation_polynomials.iter() {
            write_polynomial(p, &mut writer)?;
        }

        Ok(())
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        let n = reader.read_u64::<BigEndian>()?;
        let num_inputs = reader.read_u64::<BigEndian>()?;

        let num_selectors = reader.read_u64::<BigEndian>()?;
        let mut selectors = Vec::with_capacity(num_selectors as usize);
        for _ in 0..num_selectors {
            let poly = read_polynomial_coeffs(&mut reader)?;
            selectors.push(poly);
        }

        let num_next_step_selectors = reader.read_u64::<BigEndian>()?;
        let mut next_step_selectors = Vec::with_capacity(num_next_step_selectors as usize);
        for _ in 0..num_next_step_selectors {
            let poly = read_polynomial_coeffs(&mut reader)?;
            next_step_selectors.push(poly);
        }

        let num_permutation_polys = reader.read_u64::<BigEndian>()?;
        let mut permutation_polys = Vec::with_capacity(num_permutation_polys as usize);
        for _ in 0..num_permutation_polys {
            let poly = read_polynomial_coeffs(&mut reader)?;
            permutation_polys.push(poly);
        }

        let new = Self {
            n: n as usize,
            num_inputs: num_inputs as usize,
            selector_polynomials: selectors,
            next_step_selector_polynomials: next_step_selectors,
            permutation_polynomials: permutation_polys,

            _marker: std::marker::PhantomData,
        };

        Ok(new)
    }
}

pub struct SetupPolynomialsPrecomputations<E: Engine, P: PlonkConstraintSystemParams<E>> {
    pub selector_polynomials_on_coset_of_size_4n_bitreversed: Vec<Polynomial<E::Fr, Values>>,
    pub next_step_selector_polynomials_on_coset_of_size_4n_bitreversed:
        Vec<Polynomial<E::Fr, Values>>,
    pub permutation_polynomials_on_coset_of_size_4n_bitreversed: Vec<Polynomial<E::Fr, Values>>,
    pub permutation_polynomials_values_of_size_n_minus_one: Vec<Polynomial<E::Fr, Values>>,
    pub inverse_divisor_on_coset_of_size_4n_bitreversed: Polynomial<E::Fr, Values>,
    pub x_on_coset_of_size_4n_bitreversed: Polynomial<E::Fr, Values>,

    pub(crate) _marker: std::marker::PhantomData<P>,
}

use crate::plonk::fft::cooley_tukey_ntt::{BitReversedOmegas, CTPrecomputations};

impl<E: Engine, P: PlonkConstraintSystemParams<E>> SetupPolynomialsPrecomputations<E, P> {
    pub fn from_setup_and_precomputations<CP: CTPrecomputations<E::Fr>>(
        setup: &SetupPolynomials<E, P>,
        worker: &Worker,
        omegas_bitreversed: &CP,
    ) -> Result<Self, SynthesisError> {
        let mut new = Self {
            selector_polynomials_on_coset_of_size_4n_bitreversed: vec![],
            next_step_selector_polynomials_on_coset_of_size_4n_bitreversed: vec![],
            permutation_polynomials_on_coset_of_size_4n_bitreversed: vec![],
            permutation_polynomials_values_of_size_n_minus_one: vec![],
            inverse_divisor_on_coset_of_size_4n_bitreversed: Polynomial::from_values(vec![
                E::Fr::one(),
            ])
            .unwrap(),
            x_on_coset_of_size_4n_bitreversed: Polynomial::from_values(vec![E::Fr::one()]).unwrap(),

            _marker: std::marker::PhantomData,
        };

        let required_domain_size = setup.selector_polynomials[0].size();

        assert!(required_domain_size.is_power_of_two());
        let coset_generator = E::Fr::multiplicative_generator();

        // let coset_generator = E::Fr::one();

        // we do not precompute q_const as we need to use it for public inputs;
        for p in setup.selector_polynomials[0..(setup.selector_polynomials.len() - 1)].iter() {
            let ext = p.clone().bitreversed_lde_using_bitreversed_ntt(
                &worker,
                LDE_FACTOR,
                omegas_bitreversed,
                &coset_generator,
            )?;

            new.selector_polynomials_on_coset_of_size_4n_bitreversed
                .push(ext);
        }

        for p in setup.next_step_selector_polynomials.iter() {
            let ext = p.clone().bitreversed_lde_using_bitreversed_ntt(
                &worker,
                LDE_FACTOR,
                omegas_bitreversed,
                &coset_generator,
            )?;

            new.next_step_selector_polynomials_on_coset_of_size_4n_bitreversed
                .push(ext);
        }

        for p in setup.permutation_polynomials.iter() {
            let lde = p.clone().bitreversed_lde_using_bitreversed_ntt(
                &worker,
                LDE_FACTOR,
                omegas_bitreversed,
                &coset_generator,
            )?;
            new.permutation_polynomials_on_coset_of_size_4n_bitreversed
                .push(lde);

            let as_values = p.clone().fft(&worker);
            let mut as_values = as_values.into_coeffs();
            as_values
                .pop()
                .expect("must shorted permutation polynomial values by one");

            let p = Polynomial::from_values_unpadded(as_values)?;

            new.permutation_polynomials_values_of_size_n_minus_one
                .push(p);
        }

        let mut vanishing_poly_inverse_bitreversed =
            evaluate_vanishing_polynomial_of_degree_on_domain_size::<E::Fr>(
                required_domain_size as u64,
                &E::Fr::multiplicative_generator(),
                (required_domain_size * LDE_FACTOR) as u64,
                &worker,
            )?;
        vanishing_poly_inverse_bitreversed.batch_inversion(&worker)?;
        vanishing_poly_inverse_bitreversed.bitreverse_enumeration(&worker);

        assert_eq!(
            vanishing_poly_inverse_bitreversed.size(),
            required_domain_size * LDE_FACTOR
        );

        // evaluate polynomial X on the coset
        let mut x_poly = Polynomial::from_values(vec![
            coset_generator;
            vanishing_poly_inverse_bitreversed.size()
        ])?;
        x_poly.distribute_powers(&worker, x_poly.omega);
        x_poly.bitreverse_enumeration(&worker);

        assert_eq!(x_poly.size(), required_domain_size * LDE_FACTOR);

        new.inverse_divisor_on_coset_of_size_4n_bitreversed = vanishing_poly_inverse_bitreversed;
        new.x_on_coset_of_size_4n_bitreversed = x_poly;

        Ok(new)
    }

    pub fn from_setup(
        setup: &SetupPolynomials<E, P>,
        worker: &Worker,
    ) -> Result<Self, SynthesisError> {
        let precomps =
            BitReversedOmegas::new_for_domain_size(setup.permutation_polynomials[0].size());

        Self::from_setup_and_precomputations(setup, worker, &precomps)
    }


    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        writer.write_u64::<BigEndian>(self.selector_polynomials_on_coset_of_size_4n_bitreversed.len() as u64)?;
        for p in &self.selector_polynomials_on_coset_of_size_4n_bitreversed {
            write_polynomial(p, &mut writer)?;
        }
        writer.write_u64::<BigEndian>(self.next_step_selector_polynomials_on_coset_of_size_4n_bitreversed.len() as u64)?;
        for p in &self.next_step_selector_polynomials_on_coset_of_size_4n_bitreversed {
            write_polynomial(p, &mut writer)?;
        }
        writer.write_u64::<BigEndian>(self.permutation_polynomials_on_coset_of_size_4n_bitreversed.len() as u64)?;
        for p in &self.permutation_polynomials_on_coset_of_size_4n_bitreversed {
            write_polynomial(p, &mut writer)?;
        }
        writer.write_u64::<BigEndian>(self.permutation_polynomials_values_of_size_n_minus_one.len() as u64)?;
        for p in &self.permutation_polynomials_values_of_size_n_minus_one {
            write_polynomial(p, &mut writer)?;
        }
        write_polynomial(&self.inverse_divisor_on_coset_of_size_4n_bitreversed, &mut writer)?;
        write_polynomial(&self.x_on_coset_of_size_4n_bitreversed, &mut writer)?;
        Ok(())
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        let num_selectors = reader.read_u64::<BigEndian>()?;
        let mut selector_polynomials_on_coset_of_size_4n_bitreversed =
            Vec::with_capacity(num_selectors as usize);
        for _ in 0..num_selectors {
            let poly = read_polynomial_values_unpadded(&mut reader)?;
            selector_polynomials_on_coset_of_size_4n_bitreversed.push(poly);
        }

        let num_next_step_selectors = reader.read_u64::<BigEndian>()?;
        let mut next_step_selector_polynomials_on_coset_of_size_4n_bitreversed =
            Vec::with_capacity(num_next_step_selectors as usize);
        for _ in 0..num_next_step_selectors {
            let poly = read_polynomial_values_unpadded(&mut reader)?;
            next_step_selector_polynomials_on_coset_of_size_4n_bitreversed.push(poly);
        }

        let num_permutation_polys = reader.read_u64::<BigEndian>()?;
        let mut permutation_polynomials_on_coset_of_size_4n_bitreversed =
            Vec::with_capacity(num_permutation_polys as usize);
        for _ in 0..num_permutation_polys {
            let poly = read_polynomial_values_unpadded(&mut reader)?;
            permutation_polynomials_on_coset_of_size_4n_bitreversed.push(poly);
        }

        let num_permutation_polys_size_minus_one = reader.read_u64::<BigEndian>()?;
        let mut permutation_polynomials_values_of_size_n_minus_one = Vec::with_capacity(num_permutation_polys as usize);
        for _ in 0..num_permutation_polys_size_minus_one {
            let poly = read_polynomial_values_unpadded(&mut reader)?;
            permutation_polynomials_values_of_size_n_minus_one.push(poly);
        }
        let inverse_divisor_on_coset_of_size_4n_bitreversed = read_polynomial_values_unpadded(&mut reader)?;
        let x_on_coset_of_size_4n_bitreversed = read_polynomial_values_unpadded(&mut reader)?;

        Ok(Self {
            selector_polynomials_on_coset_of_size_4n_bitreversed,
            next_step_selector_polynomials_on_coset_of_size_4n_bitreversed,
            permutation_polynomials_on_coset_of_size_4n_bitreversed,
            permutation_polynomials_values_of_size_n_minus_one,
            inverse_divisor_on_coset_of_size_4n_bitreversed,
            x_on_coset_of_size_4n_bitreversed,
            _marker: std::marker::PhantomData,
        })
    }
}

#[derive(Clone, Debug)]
pub struct Proof<E: Engine, P: PlonkConstraintSystemParams<E>> {
    pub num_inputs: usize,
    pub n: usize,
    pub input_values: Vec<E::Fr>,
    pub wire_commitments: Vec<E::G1Affine>,
    pub grand_product_commitment: E::G1Affine,
    pub quotient_poly_commitments: Vec<E::G1Affine>,

    pub wire_values_at_z: Vec<E::Fr>,
    pub wire_values_at_z_omega: Vec<E::Fr>,
    pub grand_product_at_z_omega: E::Fr,
    pub quotient_polynomial_at_z: E::Fr,
    pub linearization_polynomial_at_z: E::Fr,
    pub permutation_polynomials_at_z: Vec<E::Fr>,

    pub opening_at_z_proof: E::G1Affine,
    pub opening_at_z_omega_proof: E::G1Affine,

    pub _marker: std::marker::PhantomData<P>,
}

impl<E: Engine, P: PlonkConstraintSystemParams<E>> Proof<E, P> {
    pub fn empty() -> Self {
        use crate::pairing::CurveAffine;

        Self {
            num_inputs: 0,
            n: 0,
            input_values: vec![],
            wire_commitments: vec![],
            grand_product_commitment: E::G1Affine::zero(),
            quotient_poly_commitments: vec![],
            wire_values_at_z: vec![],
            wire_values_at_z_omega: vec![],
            grand_product_at_z_omega: E::Fr::zero(),
            quotient_polynomial_at_z: E::Fr::zero(),
            linearization_polynomial_at_z: E::Fr::zero(),
            permutation_polynomials_at_z: vec![],

            opening_at_z_proof: E::G1Affine::zero(),
            opening_at_z_omega_proof: E::G1Affine::zero(),

            _marker: std::marker::PhantomData,
        }
    }

    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        use crate::pairing::CurveAffine;

        assert_eq!(self.num_inputs, self.input_values.len());

        writer.write_u64::<BigEndian>(self.n as u64)?;
        writer.write_u64::<BigEndian>(self.num_inputs as u64)?;

        // writer.write_u64::<BigEndian>(self.input_values.len() as u64)?;
        for p in self.input_values.iter() {
            write_fr(p, &mut writer)?;
        }

        assert_eq!(self.wire_commitments.len(), P::STATE_WIDTH);
        writer.write_u64::<BigEndian>(self.wire_commitments.len() as u64)?;
        for p in self.wire_commitments.iter() {
            writer.write_all(p.into_uncompressed().as_ref())?;
        }

        writer.write_all(self.grand_product_commitment.into_uncompressed().as_ref())?;

        writer.write_u64::<BigEndian>(self.quotient_poly_commitments.len() as u64)?;
        for p in self.quotient_poly_commitments.iter() {
            writer.write_all(p.into_uncompressed().as_ref())?;
        }

        writer.write_u64::<BigEndian>(self.wire_values_at_z.len() as u64)?;
        for p in self.wire_values_at_z.iter() {
            write_fr(p, &mut writer)?;
        }

        writer.write_u64::<BigEndian>(self.wire_values_at_z_omega.len() as u64)?;
        for p in self.wire_values_at_z_omega.iter() {
            write_fr(p, &mut writer)?;
        }

        write_fr(&self.grand_product_at_z_omega, &mut writer)?;
        write_fr(&self.quotient_polynomial_at_z, &mut writer)?;
        write_fr(&self.linearization_polynomial_at_z, &mut writer)?;

        writer.write_u64::<BigEndian>(self.permutation_polynomials_at_z.len() as u64)?;
        for p in self.permutation_polynomials_at_z.iter() {
            write_fr(p, &mut writer)?;
        }

        writer.write_all(self.opening_at_z_proof.into_uncompressed().as_ref())?;
        writer.write_all(self.opening_at_z_omega_proof.into_uncompressed().as_ref())?;

        Ok(())
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        use crate::pairing::CurveAffine;
        use crate::pairing::EncodedPoint;

        let n = reader.read_u64::<BigEndian>()?;
        let num_inputs = reader.read_u64::<BigEndian>()?;

        let read_g1 = |reader: &mut R| -> std::io::Result<E::G1Affine> {
            let mut repr = <E::G1Affine as CurveAffine>::Uncompressed::empty();
            reader.read_exact(repr.as_mut())?;

            let e = repr
                .into_affine()
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

            Ok(e)
        };

        let mut inputs = Vec::with_capacity(num_inputs as usize);
        for _ in 0..num_inputs {
            let p = read_fr(&mut reader)?;
            inputs.push(p);
        }

        let num_wire_commitments = reader.read_u64::<BigEndian>()?;
        let mut wire_commitments = Vec::with_capacity(num_wire_commitments as usize);
        for _ in 0..num_wire_commitments {
            let p = read_g1(&mut reader)?;
            wire_commitments.push(p);
        }

        let grand_product_commitment = read_g1(&mut reader)?;

        let num_quotient_commitments = reader.read_u64::<BigEndian>()?;
        let mut quotient_poly_commitments = Vec::with_capacity(num_quotient_commitments as usize);
        for _ in 0..num_quotient_commitments {
            let p = read_g1(&mut reader)?;
            quotient_poly_commitments.push(p);
        }

        let num_wire_values_at_z = reader.read_u64::<BigEndian>()?;
        let mut wire_values_at_z = Vec::with_capacity(num_wire_values_at_z as usize);
        for _ in 0..num_wire_values_at_z {
            let p = read_fr(&mut reader)?;
            wire_values_at_z.push(p);
        }

        let num_wire_values_at_z_omega = reader.read_u64::<BigEndian>()?;
        let mut wire_values_at_z_omega = Vec::with_capacity(num_wire_values_at_z_omega as usize);
        for _ in 0..num_wire_values_at_z_omega {
            let p = read_fr(&mut reader)?;
            wire_values_at_z_omega.push(p);
        }

        let grand_product_at_z_omega = read_fr(&mut reader)?;
        let quotient_polynomial_at_z = read_fr(&mut reader)?;
        let linearization_polynomial_at_z = read_fr(&mut reader)?;

        let num_perm_at_z = reader.read_u64::<BigEndian>()?;
        let mut permutation_polynomials_at_z = Vec::with_capacity(num_perm_at_z as usize);
        for _ in 0..num_perm_at_z {
            let p = read_fr(&mut reader)?;
            permutation_polynomials_at_z.push(p);
        }

        let opening_at_z_proof = read_g1(&mut reader)?;
        let opening_at_z_omega_proof = read_g1(&mut reader)?;

        let new = Self {
            num_inputs: num_inputs as usize,
            n: n as usize,
            input_values: inputs,
            wire_commitments: wire_commitments,
            grand_product_commitment: grand_product_commitment,
            quotient_poly_commitments: quotient_poly_commitments,
            wire_values_at_z: wire_values_at_z,
            wire_values_at_z_omega: wire_values_at_z_omega,
            grand_product_at_z_omega,
            quotient_polynomial_at_z,
            linearization_polynomial_at_z,
            permutation_polynomials_at_z: permutation_polynomials_at_z,

            opening_at_z_proof: opening_at_z_proof,
            opening_at_z_omega_proof: opening_at_z_omega_proof,

            _marker: std::marker::PhantomData,
        };

        Ok(new)
    }
}

#[derive(Clone, Debug)]
pub struct VerificationKey<E: Engine, P: PlonkConstraintSystemParams<E>> {
    pub n: usize,
    pub num_inputs: usize,
    pub selector_commitments: Vec<E::G1Affine>,
    pub next_step_selector_commitments: Vec<E::G1Affine>,
    pub permutation_commitments: Vec<E::G1Affine>,
    pub non_residues: Vec<E::Fr>,

    pub g2_elements: [E::G2Affine; 2],

    pub _marker: std::marker::PhantomData<P>,
}

impl<E: Engine, P: PlonkConstraintSystemParams<E>> VerificationKey<E, P> {
    pub fn from_setup(
        setup: &SetupPolynomials<E, P>,
        worker: &Worker,
        crs: &Crs<E, CrsForMonomialForm>,
    ) -> Result<Self, SynthesisError> {
        assert_eq!(setup.selector_polynomials.len(), P::STATE_WIDTH + 2);
        if P::CAN_ACCESS_NEXT_TRACE_STEP == false {
            assert_eq!(setup.next_step_selector_polynomials.len(), 0);
        }
        assert_eq!(setup.permutation_polynomials.len(), P::STATE_WIDTH);

        let mut new = Self {
            n: setup.n,
            num_inputs: setup.num_inputs,
            selector_commitments: vec![],
            next_step_selector_commitments: vec![],
            permutation_commitments: vec![],
            non_residues: vec![],

            g2_elements: [crs.g2_monomial_bases[0], crs.g2_monomial_bases[1]],

            _marker: std::marker::PhantomData,
        };

        for p in setup.selector_polynomials.iter() {
            let commitment = commit_using_monomials(p, &crs, &worker)?;
            new.selector_commitments.push(commitment);
        }

        for p in setup.next_step_selector_polynomials.iter() {
            let commitment = commit_using_monomials(p, &crs, &worker)?;
            new.next_step_selector_commitments.push(commitment);
        }

        for p in setup.permutation_polynomials.iter() {
            let commitment = commit_using_monomials(p, &crs, &worker)?;
            new.permutation_commitments.push(commitment);
        }

        new.non_residues
            .extend(super::utils::make_non_residues::<E::Fr>(P::STATE_WIDTH - 1));

        Ok(new)
    }

    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        use crate::pairing::CurveAffine;

        writer.write_u64::<BigEndian>(self.n as u64)?;
        writer.write_u64::<BigEndian>(self.num_inputs as u64)?;

        writer.write_u64::<BigEndian>(self.selector_commitments.len() as u64)?;
        for p in self.selector_commitments.iter() {
            writer.write_all(p.into_uncompressed().as_ref())?;
        }

        writer.write_u64::<BigEndian>(self.next_step_selector_commitments.len() as u64)?;
        for p in self.next_step_selector_commitments.iter() {
            writer.write_all(p.into_uncompressed().as_ref())?;
        }

        writer.write_u64::<BigEndian>(self.permutation_commitments.len() as u64)?;
        for p in self.permutation_commitments.iter() {
            writer.write_all(p.into_uncompressed().as_ref())?;
        }

        writer.write_u64::<BigEndian>(self.non_residues.len() as u64)?;
        for p in self.non_residues.iter() {
            write_fr(p, &mut writer)?;
        }

        writer.write_all(self.g2_elements[0].into_uncompressed().as_ref())?;
        writer.write_all(self.g2_elements[1].into_uncompressed().as_ref())?;

        Ok(())
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        use crate::pairing::CurveAffine;
        use crate::pairing::EncodedPoint;

        let n = reader.read_u64::<BigEndian>()?;
        let num_inputs = reader.read_u64::<BigEndian>()?;

        let read_g1 = |reader: &mut R| -> std::io::Result<E::G1Affine> {
            let mut repr = <E::G1Affine as CurveAffine>::Uncompressed::empty();
            reader.read_exact(repr.as_mut())?;

            let e = repr
                .into_affine()
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

            Ok(e)
        };

        let read_g2_not_zero = |reader: &mut R| -> std::io::Result<E::G2Affine> {
            let mut repr = <E::G2Affine as CurveAffine>::Uncompressed::empty();
            reader.read_exact(repr.as_mut())?;

            let e = repr
                .into_affine()
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
                .and_then(|e| {
                    if e.is_zero() {
                        Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "point at infinity",
                        ))?
                    } else {
                        Ok(e)
                    }
                });

            e
        };

        let num_selectors = reader.read_u64::<BigEndian>()?;
        let mut selectors = Vec::with_capacity(num_selectors as usize);
        for _ in 0..num_selectors {
            let p = read_g1(&mut reader)?;
            selectors.push(p);
        }

        let num_next_step_selectors = reader.read_u64::<BigEndian>()?;
        let mut next_step_selectors = Vec::with_capacity(num_next_step_selectors as usize);
        for _ in 0..num_next_step_selectors {
            let p = read_g1(&mut reader)?;
            next_step_selectors.push(p);
        }

        let num_permutation_polys = reader.read_u64::<BigEndian>()?;
        let mut permutation_polys = Vec::with_capacity(num_permutation_polys as usize);
        for _ in 0..num_permutation_polys {
            let p = read_g1(&mut reader)?;
            permutation_polys.push(p);
        }

        let num_non_residues = reader.read_u64::<BigEndian>()?;
        let mut non_residues = Vec::with_capacity(num_non_residues as usize);
        for _ in 0..num_non_residues {
            let p = read_fr(&mut reader)?;
            non_residues.push(p);
        }

        let g2_points = [
            read_g2_not_zero(&mut reader)?,
            read_g2_not_zero(&mut reader)?,
        ];

        let new = Self {
            n: n as usize,
            num_inputs: num_inputs as usize,
            selector_commitments: selectors,
            next_step_selector_commitments: next_step_selectors,
            permutation_commitments: permutation_polys,
            non_residues: non_residues,

            g2_elements: g2_points,

            _marker: std::marker::PhantomData,
        };

        Ok(new)
    }
}