use super::cs::*;
use super::data_structures::*;
use crate::pairing::ff::*;
use crate::pairing::*;
use crate::plonk::polynomials::*;
use std::collections::HashMap;

use crate::plonk::domains::*;
use crate::worker::Worker;
use crate::SynthesisError;

use crate::kate_commitment::*;

use super::super::better_cs::utils::make_non_residues;

use crate::byteorder::BigEndian;
use crate::byteorder::ReadBytesExt;
use crate::byteorder::WriteBytesExt;
use std::io::{Read, Write};

use crate::plonk::better_cs::keys::*;

#[derive(Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Setup<E: Engine, C: Circuit<E>> {
    pub n: usize,
    pub num_inputs: usize,
    pub state_width: usize,
    pub num_witness_polys: usize,

    pub gate_setup_monomials: Vec<Polynomial<E::Fr, Coefficients>>,
    pub gate_selectors_monomials: Vec<Polynomial<E::Fr, Coefficients>>,
    pub permutation_monomials: Vec<Polynomial<E::Fr, Coefficients>>,

    pub total_lookup_entries_length: usize,
    pub lookup_selector_monomial: Option<Polynomial<E::Fr, Coefficients>>,
    pub lookup_tables_monomials: Vec<Polynomial<E::Fr, Coefficients>>,
    pub lookup_table_type_monomial: Option<Polynomial<E::Fr, Coefficients>>,

    pub non_residues: Vec<E::Fr>,

    #[serde(skip_serializing,skip_deserializing, default)]
    #[serde(bound(serialize = ""))]
    #[serde(bound(deserialize = ""))]
    _marker: std::marker::PhantomData<C>
}

impl<E: Engine, C: Circuit<E>> std::fmt::Debug for Setup<E, C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Setup")
            .field("n", &self.n)
            .field("num_inputs", &self.num_inputs)
            .field("gate_setup_monomials", &self.gate_setup_monomials)
            .field("gate_selectors_monomials", &self.gate_selectors_monomials)
            .field("permutation_monomials", &self.permutation_monomials)
            .field("total_lookup_entries_length", &self.total_lookup_entries_length)
            .field("lookup_selector_monomial", &self.lookup_selector_monomial)
            .field("lookup_tables_monomials", &self.lookup_tables_monomials)
            .field("lookup_table_type_monomial", &self.lookup_table_type_monomial)
            .finish()
    }
}

impl<E: Engine, C: Circuit<E>> Setup<E, C> {
    pub fn empty() -> Self {
        Self {
            n: 0,
            num_inputs: 0,
            state_width: 0,
            num_witness_polys: 0,
            gate_setup_monomials: vec![],
            gate_selectors_monomials: vec![],
            permutation_monomials: vec![],
        
            total_lookup_entries_length: 0,
            lookup_selector_monomial: None,
            lookup_tables_monomials: vec![],
            lookup_table_type_monomial: None,
            non_residues: vec![],
        
            _marker: std::marker::PhantomData
        }
    }

    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        writer.write_u64::<BigEndian>(self.n as u64)?;
        writer.write_u64::<BigEndian>(self.num_inputs as u64)?;
        writer.write_u64::<BigEndian>(self.state_width as u64)?;
        writer.write_u64::<BigEndian>(self.num_witness_polys as u64)?;

        write_polynomials_vec(&self.gate_setup_monomials, &mut writer)?;
        write_polynomials_vec(&self.gate_selectors_monomials, &mut writer)?;
        write_polynomials_vec(&self.permutation_monomials, &mut writer)?;

        writer.write_u64::<BigEndian>(self.total_lookup_entries_length as u64)?;
        write_optional_polynomial(&self.lookup_selector_monomial, &mut writer)?;
        write_polynomials_vec(&self.lookup_tables_monomials, &mut writer)?;
        write_optional_polynomial(&self.lookup_table_type_monomial, &mut writer)?;

        write_fr_vec(&self.non_residues, &mut writer)?;

        Ok(())
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        use crate::pairing::CurveAffine;
        use crate::pairing::EncodedPoint;

        let n = reader.read_u64::<BigEndian>()?;
        let num_inputs = reader.read_u64::<BigEndian>()?;
        let state_width = reader.read_u64::<BigEndian>()?;
        let num_witness_polys = reader.read_u64::<BigEndian>()?;

        let gate_setup_monomials = read_polynomials_coeffs_vec(&mut reader)?;
        let gate_selectors_monomials = read_polynomials_coeffs_vec(&mut reader)?;
        let permutation_monomials = read_polynomials_coeffs_vec(&mut reader)?;

        let total_lookup_entries_length = reader.read_u64::<BigEndian>()?;
        let lookup_selector_monomial = read_optional_polynomial_coeffs(&mut reader)?;
        let lookup_tables_monomials = read_polynomials_coeffs_vec(&mut reader)?;
        let lookup_table_type_monomial = read_optional_polynomial_coeffs(&mut reader)?;

        let non_residues = read_fr_vec(&mut reader)?;

        let new = Self {
            n: n as usize,
            num_inputs: num_inputs as usize,
            state_width: state_width as usize,
            num_witness_polys: num_witness_polys as usize,
            gate_setup_monomials,
            gate_selectors_monomials,
            permutation_monomials,
            total_lookup_entries_length: total_lookup_entries_length as usize,
            lookup_selector_monomial,
            lookup_tables_monomials,
            lookup_table_type_monomial,
            non_residues,

            _marker: std::marker::PhantomData,
        };

        Ok(new)
    }
}

#[derive(Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct VerificationKey<E: Engine, C: Circuit<E>> {
    pub n: usize,
    pub num_inputs: usize,
    pub state_width: usize,
    pub num_witness_polys: usize,

    pub gate_setup_commitments: Vec<E::G1Affine>,
    pub gate_selectors_commitments: Vec<E::G1Affine>,
    pub permutation_commitments: Vec<E::G1Affine>,

    pub total_lookup_entries_length: usize,
    pub lookup_selector_commitment: Option<E::G1Affine>,
    pub lookup_tables_commitments: Vec<E::G1Affine>,
    pub lookup_table_type_commitment: Option<E::G1Affine>,

    pub non_residues: Vec<E::Fr>,
    pub g2_elements: [E::G2Affine; 2],

    #[serde(skip_serializing,skip_deserializing, default)]
    #[serde(bound(serialize = ""))]
    #[serde(bound(deserialize = ""))]
    _marker: std::marker::PhantomData<C>
}

impl<E: Engine, C: Circuit<E>> std::fmt::Debug for VerificationKey<E, C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VerificationKey")
            .field("n", &self.n)
            .field("num_inputs", &self.num_inputs)
            .field("gate_setup_commitments", &self.gate_setup_commitments)
            .field("gate_selectors_commitments", &self.gate_selectors_commitments)
            .field("permutation_commitments", &self.permutation_commitments)
            .field("total_lookup_entries_length", &self.total_lookup_entries_length)
            .field("lookup_selector_commitment", &self.lookup_selector_commitment)
            .field("lookup_tables_commitments", &self.lookup_tables_commitments)
            .field("lookup_table_type_commitment", &self.lookup_table_type_commitment)
            .finish()
    }
}

impl<E: Engine, C: Circuit<E>> VerificationKey<E, C> {
    pub fn empty() -> Self {
        Self {
            n: 0,
            num_inputs: 0,
            state_width: 0,
            num_witness_polys: 0,
            gate_setup_commitments: vec![],
            gate_selectors_commitments: vec![],
            permutation_commitments: vec![],

            total_lookup_entries_length: 0,
            lookup_selector_commitment: None,
            lookup_tables_commitments: vec![],
            lookup_table_type_commitment: None,

            non_residues: vec![],
            g2_elements: [<E::G2Affine as pairing::CurveAffine>::zero(); 2],

            _marker: std::marker::PhantomData,
        }
    }

    pub fn from_setup(
        setup: &Setup<E, C>,
        worker: &Worker,
        crs: &Crs<E, CrsForMonomialForm>,
    ) -> Result<Self, SynthesisError> {
        let mut new = Self {
            n: setup.n,
            num_inputs: setup.num_inputs,
            state_width: setup.state_width,
            num_witness_polys: setup.num_witness_polys,
            gate_setup_commitments: vec![],
            gate_selectors_commitments: vec![],
            permutation_commitments: vec![],

            total_lookup_entries_length: setup.total_lookup_entries_length,
            lookup_selector_commitment: None,
            lookup_tables_commitments: vec![],
            lookup_table_type_commitment: None,
        
            non_residues: vec![],
            g2_elements: [crs.g2_monomial_bases[0], crs.g2_monomial_bases[1]],

            _marker: std::marker::PhantomData,
        };

        for (p, c) in vec![
            (&setup.gate_setup_monomials, &mut new.gate_setup_commitments),
            (&setup.gate_selectors_monomials, &mut new.gate_selectors_commitments),
            (&setup.permutation_monomials, &mut new.permutation_commitments),
            (&setup.lookup_tables_monomials, &mut new.lookup_tables_commitments),
        ].into_iter() {
            for p in p.iter() {
                let commitment = commit_using_monomials(p, &crs, &worker)?;
                c.push(commitment);
            }
        }

        if let Some(p) = setup.lookup_selector_monomial.as_ref() {
            let commitment = commit_using_monomials(p, &crs, &worker)?;
            new.lookup_selector_commitment = Some(commitment);
        }

        if let Some(p) = setup.lookup_table_type_monomial.as_ref() {
            let commitment = commit_using_monomials(p, &crs, &worker)?;
            new.lookup_table_type_commitment = Some(commitment);
        }

        new.non_residues = setup.non_residues.clone();

        // new.non_residues
        //     .extend(make_non_residues::<E::Fr>(state_width - 1));

        Ok(new)
    }

    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        writer.write_u64::<BigEndian>(self.n as u64)?;
        writer.write_u64::<BigEndian>(self.num_inputs as u64)?;
        writer.write_u64::<BigEndian>(self.state_width as u64)?;
        writer.write_u64::<BigEndian>(self.num_witness_polys as u64)?;

        write_curve_affine_vec(&self.gate_setup_commitments, &mut writer)?;
        write_curve_affine_vec(&self.gate_selectors_commitments, &mut writer)?;
        write_curve_affine_vec(&self.permutation_commitments, &mut writer)?;

        writer.write_u64::<BigEndian>(self.total_lookup_entries_length as u64)?;
        write_optional_curve_affine(&self.lookup_selector_commitment, &mut writer)?;
        write_curve_affine_vec(&self.lookup_tables_commitments, &mut writer)?;
        write_optional_curve_affine(&self.lookup_table_type_commitment, &mut writer)?;

        write_fr_vec(&self.non_residues, &mut writer)?;

        write_curve_affine(&self.g2_elements[0], &mut writer)?;
        write_curve_affine(&self.g2_elements[1], &mut writer)?;

        Ok(())
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        use crate::pairing::CurveAffine;
        use crate::pairing::EncodedPoint;

        let n = reader.read_u64::<BigEndian>()?;
        let num_inputs = reader.read_u64::<BigEndian>()?;
        let state_width = reader.read_u64::<BigEndian>()?;
        let num_witness_polys = reader.read_u64::<BigEndian>()?;

        let gate_setup_commitments = read_curve_affine_vector(&mut reader)?;
        let gate_selectors_commitments = read_curve_affine_vector(&mut reader)?;
        let permutation_commitments = read_curve_affine_vector(&mut reader)?;

        let total_lookup_entries_length = reader.read_u64::<BigEndian>()?;
        let lookup_selector_commitment = read_optional_curve_affine(&mut reader)?;
        let lookup_tables_commitments = read_curve_affine_vector(&mut reader)?;
        let lookup_table_type_commitment = read_optional_curve_affine(&mut reader)?;

        let non_residues = read_fr_vec(&mut reader)?;

        let h = read_curve_affine(&mut reader)?;
        let h_x = read_curve_affine(&mut reader)?;

        let new = Self {
            n: n as usize,
            num_inputs: num_inputs as usize,
            state_width: state_width as usize,
            num_witness_polys: num_witness_polys as usize,
            gate_setup_commitments,
            gate_selectors_commitments,
            permutation_commitments,
            total_lookup_entries_length: total_lookup_entries_length as usize,
            lookup_selector_commitment,
            lookup_tables_commitments,
            lookup_table_type_commitment,
            non_residues,

            g2_elements: [h, h_x],

            _marker: std::marker::PhantomData,
        };

        Ok(new)
    }
}

use super::data_structures::AssembledPolynomialStorageForMonomialForms;

impl<'a, E: Engine> AssembledPolynomialStorageForMonomialForms<'a, E> {
    pub fn extend_from_setup<C: Circuit<E>>(&mut self, setup: &'a Setup<E, C>) -> Result<(), SynthesisError> {
        // extend with gate setup polys, gate selectors, permutation polys
        // and lookup table setup polys if available
        let all_gates = C::declare_used_gates()?;

        let has_selectors = all_gates.len() > 1;

        let mut setup_gates_iter = setup.gate_setup_monomials.iter();
        for gate in all_gates.iter() {
            for &poly_id in gate.setup_polynomials().into_iter() {
                let poly_ref = setup_gates_iter.next().expect(&format!("must have gate setup poly {:?} for gate {:?} in setup", poly_id, gate));
                let proxy = PolynomialProxy::from_borrowed(poly_ref);
                self.setup_map.insert(poly_id, proxy);
            }
        }

        assert!(setup_gates_iter.next().is_none());

        if has_selectors {
            let mut selector_iter = setup.gate_selectors_monomials.iter();
            for gate in all_gates.into_iter() {
                let id = PolyIdentifier::GateSelector(gate.name());
                let poly_ref = selector_iter.next().expect(&format!("must have gate selector poly for gate {:?} in setup", gate));
                let proxy = PolynomialProxy::from_borrowed(poly_ref);
                self.gate_selectors.insert(id, proxy);
            }
            assert!(selector_iter.next().is_none());
        }

        for (idx, poly_ref) in setup.permutation_monomials.iter().enumerate() {
            let id = PolyIdentifier::PermutationPolynomial(idx);
            let proxy = PolynomialProxy::from_borrowed(poly_ref);
            self.setup_map.insert(id, proxy);
        }

        Ok(())
    }   
}