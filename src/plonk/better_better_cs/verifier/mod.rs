use crate::pairing::ff::*;
use crate::pairing::*;

use super::setup::VerificationKey;
use super::proof::{Proof, sort_queries_for_linearization};
use super::cs::*;
use super::data_structures::*;

use std::collections::HashMap;

use crate::plonk::domains::*;
use crate::SynthesisError;

use crate::plonk::commitments::transcript::*;
use crate::kate_commitment::*;
use crate::plonk::better_cs::utils::*;
use super::lookup_tables::LookupQuery;

pub const MAX_DILATION: usize = 1;

pub fn verify<E: Engine, C: Circuit<E>, T: Transcript<E::Fr>>(
    vk: &VerificationKey<E, C>, 
    proof: &Proof<E, C>,
    transcript_params: Option<T::InitializationParameters>,
) -> Result<bool, SynthesisError> {
    let ((pair_with_generator, pair_with_x), success) = aggregate::<_, _, T>(vk, proof, transcript_params)?;
    if !success {
        return Ok(false)
    }

    let valid = E::final_exponentiation(
        &E::miller_loop(&[
            (&pair_with_generator.prepare(), &vk.g2_elements[0].prepare()),
            (&pair_with_x.prepare(), &vk.g2_elements[1].prepare())
        ])
    ).ok_or(SynthesisError::Unsatisfiable)? == E::Fqk::one();

    Ok(valid)
}

fn safe_assert(must_be_true: bool) -> Result<(), SynthesisError> {
    if !must_be_true {
        return Err(SynthesisError::AssignmentMissing);
    }

    Ok(())
}

fn safe_assert_eq<T: Eq>(a: T, b: T) -> Result<(), SynthesisError> {
    safe_assert(a==b)
}

pub fn aggregate<E: Engine, C: Circuit<E>, T: Transcript<E::Fr>>(
    vk: &VerificationKey<E, C>, 
    proof: &Proof<E, C>,
    transcript_params: Option<T::InitializationParameters>,
) -> Result<((E::G1Affine, E::G1Affine), bool), SynthesisError> {
    let mut transcript = if let Some(params) = transcript_params {
        T::new_from_params(params)
    } else {
        T::new()
    };

    let sorted_gates = C::declare_used_gates()?;
    let num_different_gates = sorted_gates.len();

    safe_assert((vk.n+1).is_power_of_two())?;
    let required_domain_size = vk.n.next_power_of_two();

    let domain = Domain::<E::Fr>::new_for_size(required_domain_size as u64)?;

    for inp in proof.inputs.iter() {
        transcript.commit_field_element(inp);
    }

    for idx in 0..vk.state_width {
        let commitment = proof.state_polys_commitments.get(idx).ok_or(SynthesisError::AssignmentMissing)?;
        commit_point_as_xy::<E, T>(&mut transcript, commitment);
    }    

    for idx in 0..vk.num_witness_polys {
        let commitment = proof.witness_polys_commitments.get(idx).ok_or(SynthesisError::AssignmentMissing)?;
        commit_point_as_xy::<E, T>(&mut transcript, commitment);
    }   

    let mut eta = E::Fr::zero();
    if vk.total_lookup_entries_length > 0 {
        eta = transcript.get_challenge();

        let commitment = proof.lookup_s_poly_commitment.as_ref().ok_or(SynthesisError::AssignmentMissing)?;
        commit_point_as_xy::<E, T>(&mut transcript, commitment);
    }

    let beta_for_copy_permutation = transcript.get_challenge();
    let gamma_for_copy_permutation = transcript.get_challenge();

    let commitment = &proof.copy_permutation_grand_product_commitment;
    commit_point_as_xy::<E, T>(&mut transcript, commitment);

    let mut beta_for_lookup = None; 
    let mut gamma_for_lookup = None;

    if vk.total_lookup_entries_length > 0 {
        let beta_for_lookup_permutation = transcript.get_challenge();
        let gamma_for_lookup_permutation = transcript.get_challenge();

        beta_for_lookup = Some(beta_for_lookup_permutation);
        gamma_for_lookup = Some(gamma_for_lookup_permutation);

        let commitment = proof.lookup_grand_product_commitment.as_ref().ok_or(SynthesisError::AssignmentMissing)?;
        commit_point_as_xy::<E, T>(&mut transcript, commitment);
    }

    let alpha = transcript.get_challenge();

    let mut total_powers_of_alpha_for_gates = 0;
    for g in sorted_gates.iter() {
        total_powers_of_alpha_for_gates += g.num_quotient_terms();
    }

    // println!("Have {} terms from {} gates", total_powers_of_alpha_for_gates, sorted_gates.len());

    let mut current_alpha = E::Fr::one();
    let mut powers_of_alpha_for_gates = Vec::with_capacity(total_powers_of_alpha_for_gates);
    powers_of_alpha_for_gates.push(current_alpha);
    for _ in 1..total_powers_of_alpha_for_gates {
        current_alpha.mul_assign(&alpha);
        powers_of_alpha_for_gates.push(current_alpha);
    }

    safe_assert_eq(powers_of_alpha_for_gates.len(), total_powers_of_alpha_for_gates)?;

    let copy_grand_product_alphas;

    {
        current_alpha.mul_assign(&alpha);
        let alpha_0 = current_alpha;

        current_alpha.mul_assign(&alpha);

        let alpha_1 = current_alpha;

        copy_grand_product_alphas = Some([alpha_0, alpha_1]);
    }

    let mut lookup_grand_product_alphas = None;

    if vk.total_lookup_entries_length > 0 {
        current_alpha.mul_assign(&alpha);

        let alpha_0 = current_alpha;

        current_alpha.mul_assign(&alpha);

        let alpha_1 = current_alpha;

        current_alpha.mul_assign(&alpha);

        let alpha_2 = current_alpha;

        lookup_grand_product_alphas = Some([alpha_0, alpha_1, alpha_2]);
    }

    for commitment in proof.quotient_poly_parts_commitments.iter() {
        commit_point_as_xy::<E, T>(&mut transcript, commitment);
    }

    let z = transcript.get_challenge();

    let z_in_domain_size = z.pow(&[required_domain_size as u64]);

    let quotient_at_z = proof.quotient_poly_opening_at_z;
    transcript.commit_field_element(&quotient_at_z);

    // first reconstruct storage of all the commitments

    let mut setup_commitments_storage = HashMap::new();
    let mut gate_selectors_commitments_storage = HashMap::new();
    {
        let mut gate_setup_polys_commitments_iter = vk.gate_setup_commitments.iter();
        let mut gate_selectors_polys_commitments_iter = vk.gate_selectors_commitments.iter();
        println!("{:?}",vk.gate_selectors_commitments);
        for gate in sorted_gates.iter() {
            let key = PolyIdentifier::GateSelector(gate.name());
            let commitment = *gate_selectors_polys_commitments_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
            gate_selectors_commitments_storage.insert(key, commitment);

            let setup_polys = gate.setup_polynomials();
            for id in setup_polys.into_iter() {
                let commitment = *gate_setup_polys_commitments_iter.next().ok_or(SynthesisError::AssignmentMissing)?;

                setup_commitments_storage.insert(id, commitment);
            }
        }

        safe_assert(gate_setup_polys_commitments_iter.next().is_none())?;
        safe_assert(gate_selectors_polys_commitments_iter.next().is_none())?;
    }

    let queries_with_linearization = sort_queries_for_linearization(&sorted_gates, MAX_DILATION);

    let mut query_values_map = std::collections::HashMap::new();

    let mut state_polys_openings_at_z_iter = proof.state_polys_openings_at_z.iter();
    let mut state_polys_openings_at_dilations_iter = proof.state_polys_openings_at_dilations.iter();

    let mut all_values_queried_at_z = vec![];
    let mut all_values_queried_at_z_omega = vec![];

    let mut all_commitments_queried_at_z = vec![];
    let mut all_commitments_queried_at_z_omega = vec![];

    for (dilation_value, ids) in queries_with_linearization.state_polys.iter().enumerate() {
        safe_assert(dilation_value <= MAX_DILATION)?;
        for id in ids.into_iter() {
            let poly_idx = if let PolyIdentifier::VariablesPolynomial(idx) = id {
                idx
            } else {
                unreachable!();
            };

            let commitment = *proof.state_polys_commitments.get(*poly_idx).ok_or(SynthesisError::AssignmentMissing)?;

            let value = if dilation_value == 0 {
                let value = *state_polys_openings_at_z_iter.next().ok_or(SynthesisError::AssignmentMissing)?;

                all_values_queried_at_z.push(value);
                all_commitments_queried_at_z.push(commitment);

                value
            } else {
                let (dilation, state_poly_idx, value) = *state_polys_openings_at_dilations_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
                safe_assert_eq(dilation, dilation_value)?;
                safe_assert_eq(*poly_idx, state_poly_idx)?;
                safe_assert(state_poly_idx < vk.state_width)?;

                all_values_queried_at_z_omega.push(value);
                all_commitments_queried_at_z_omega.push(commitment);
                
                value
            };

            transcript.commit_field_element(&value);

            let key = PolynomialInConstraint::from_id_and_dilation(*id, dilation_value);

            query_values_map.insert(key, value);
        }
    }

    safe_assert(state_polys_openings_at_z_iter.next().is_none())?;
    safe_assert(state_polys_openings_at_dilations_iter.next().is_none())?;

    let mut witness_polys_openings_at_z_iter = proof.witness_polys_openings_at_z.iter();
    let mut witness_polys_openings_at_dilations_iter = proof.witness_polys_openings_at_dilations.iter();

    for (dilation_value, ids) in queries_with_linearization.witness_polys.iter().enumerate() {
        safe_assert(dilation_value <= MAX_DILATION)?;
        for id in ids.into_iter() {
            let poly_idx = if let PolyIdentifier::WitnessPolynomial(idx) = id {
                idx
            } else {
                unreachable!();
            };

            let commitment = *proof.witness_polys_commitments.get(*poly_idx).ok_or(SynthesisError::AssignmentMissing)?;

            let value = if dilation_value == 0 {
                let value = *witness_polys_openings_at_z_iter.next().ok_or(SynthesisError::AssignmentMissing)?;

                all_values_queried_at_z.push(value);
                all_commitments_queried_at_z.push(commitment);

                value
            } else {
                let (dilation, witness_poly_idx, value) = *witness_polys_openings_at_dilations_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
                safe_assert_eq(dilation, dilation_value)?;
                safe_assert_eq(*poly_idx, witness_poly_idx)?;
                safe_assert(witness_poly_idx < vk.num_witness_polys)?;
                
                all_values_queried_at_z_omega.push(value);
                all_commitments_queried_at_z_omega.push(commitment);
                
                value
            };

            transcript.commit_field_element(&value);

            let key = PolynomialInConstraint::from_id_and_dilation(*id, dilation_value);

            query_values_map.insert(key, value);
        }
    }

    safe_assert(witness_polys_openings_at_z_iter.next().is_none())?;
    safe_assert(witness_polys_openings_at_dilations_iter.next().is_none())?;

    let mut gate_setup_openings_at_z_iter = proof.gate_setup_openings_at_z.iter();

    for (gate_idx, queries) in queries_with_linearization.gate_setup_polys.iter().enumerate() {
        for (dilation_value, ids) in queries.iter().enumerate() {
            safe_assert(dilation_value <= MAX_DILATION)?;
            for id in ids.into_iter() {
                let poly_idx = if let PolyIdentifier::GateSetupPolynomial(_, idx) = id {
                    idx
                } else {
                    unreachable!();
                };

                let commitment = *setup_commitments_storage.get(&id).ok_or(SynthesisError::AssignmentMissing)?;

                let value = if dilation_value == 0 {
                    let (gate_index, setup_poly_index, value) = *gate_setup_openings_at_z_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
                    safe_assert_eq(gate_idx, gate_index)?;
                    safe_assert_eq(*poly_idx, setup_poly_index)?;

                    all_values_queried_at_z.push(value);
                    all_commitments_queried_at_z.push(commitment);
    
                    value
                } else {
                    unimplemented!("gate setup polynomials can not be time dilated");
                };

                transcript.commit_field_element(&value);

                let key = PolynomialInConstraint::from_id_and_dilation(*id, dilation_value);

                query_values_map.insert(key, value);
            }
        } 
    }

    safe_assert(gate_setup_openings_at_z_iter.next().is_none())?;

    // also open gate selectors
    let mut selector_values_iter = proof.gate_selectors_openings_at_z.iter();
    let mut selector_values = vec![];
    for s in queries_with_linearization.gate_selectors.iter() {
        let gate_index = sorted_gates.iter().position(|r| r == s).ok_or(SynthesisError::AssignmentMissing)?;

        let (gate_idx, value) = *selector_values_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
        safe_assert_eq(gate_index, gate_idx)?;
        transcript.commit_field_element(&value);

        let key = PolyIdentifier::GateSelector(s.name());

        // let commitment = *selector_commitments_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
        let commitment = *gate_selectors_commitments_storage.get(&key).ok_or(SynthesisError::AssignmentMissing)?;

        selector_values.push(value);
        all_values_queried_at_z.push(value);
        all_commitments_queried_at_z.push(commitment);
    }

    safe_assert(selector_values_iter.next().is_none())?;

    // copy-permutation polynomials queries

    let mut copy_permutation_polys_openings_at_z_iter = proof.copy_permutation_polys_openings_at_z.iter();
    let mut copy_permutation_polys_commitments_iter = vk.permutation_commitments.iter();

    let mut copy_permutation_queries = vec![];

    for _ in 0..(vk.state_width-1) {
        let value = *copy_permutation_polys_openings_at_z_iter.next().ok_or(SynthesisError::AssignmentMissing)?;

        transcript.commit_field_element(&value);

        copy_permutation_queries.push(value);
        all_values_queried_at_z.push(value);

        let commitment = *copy_permutation_polys_commitments_iter.next().ok_or(SynthesisError::AssignmentMissing)?;
        all_commitments_queried_at_z.push(commitment);
    }

    safe_assert(copy_permutation_polys_openings_at_z_iter.next().is_none())?;
    
    // copy-permutation grand product query

    let mut z_omega = z;
    z_omega.mul_assign(&domain.generator);

    // for polys below we will insert queried commitments manually into the corresponding lists
    let copy_permutation_z_at_z_omega = proof.copy_permutation_grand_product_opening_at_z_omega;
    transcript.commit_field_element(&copy_permutation_z_at_z_omega);

    if vk.total_lookup_entries_length > 0 {
        // first commit values at z, and then at z*omega
        transcript.commit_field_element(proof.lookup_t_poly_opening_at_z.as_ref().ok_or(SynthesisError::AssignmentMissing)?);
        transcript.commit_field_element(proof.lookup_selector_poly_opening_at_z.as_ref().ok_or(SynthesisError::AssignmentMissing)?);
        transcript.commit_field_element(proof.lookup_table_type_poly_opening_at_z.as_ref().ok_or(SynthesisError::AssignmentMissing)?);

        // now at z*omega
        transcript.commit_field_element(proof.lookup_s_poly_opening_at_z_omega.as_ref().ok_or(SynthesisError::AssignmentMissing)?);
        transcript.commit_field_element(proof.lookup_grand_product_opening_at_z_omega.as_ref().ok_or(SynthesisError::AssignmentMissing)?);
        transcript.commit_field_element(proof.lookup_t_poly_opening_at_z_omega.as_ref().ok_or(SynthesisError::AssignmentMissing)?);
    }

    let linearization_at_z = proof.linearization_poly_opening_at_z;
    transcript.commit_field_element(&linearization_at_z);

    // linearization is done, now perform sanity check
    // this is effectively a verification procedure

    let mut lookup_query = None;

    {
        let vanishing_at_z = evaluate_vanishing_for_size(&z, required_domain_size as u64);

        // first let's aggregate gates

        let mut t_num_on_full_domain = E::Fr::zero();

        let challenges_slice = &powers_of_alpha_for_gates[..];

        let mut all_gates = sorted_gates.clone();

        // we've suffered and linearization polynomial captures all the gates except the public input!

        {
            let mut tmp = linearization_at_z;
            // add input values

            let gate = all_gates.drain(0..1).into_iter().next().ok_or(SynthesisError::AssignmentMissing)?;
            safe_assert(gate.benefits_from_linearization())?;
            safe_assert(C::MainGate::default().into_internal() == gate)?;
            let gate = C::MainGate::default();
            let num_challenges = gate.num_quotient_terms();
            let (for_gate, _) = challenges_slice.split_at(num_challenges);

            let input_values = proof.inputs.clone();

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

            for idx in 0..(vk.state_width-1) {
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

            let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(vk.state_width-1));
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
            if vk.total_lookup_entries_length > 0 {

                let [alpha_0, alpha_1, alpha_2] = lookup_grand_product_alphas.expect("there must be powers of alpha for lookup permutation");

                let lookup_queries = LookupQuery::<E> {
                    s_at_z_omega: proof.lookup_s_poly_opening_at_z_omega.ok_or(SynthesisError::AssignmentMissing)?,
                    grand_product_at_z_omega: proof.lookup_grand_product_opening_at_z_omega.ok_or(SynthesisError::AssignmentMissing)?,
                    t_at_z: proof.lookup_t_poly_opening_at_z.ok_or(SynthesisError::AssignmentMissing)?,
                    t_at_z_omega: proof.lookup_t_poly_opening_at_z_omega.ok_or(SynthesisError::AssignmentMissing)?,
                    selector_at_z: proof.lookup_selector_poly_opening_at_z.ok_or(SynthesisError::AssignmentMissing)?,
                    table_type_at_z: proof.lookup_table_type_poly_opening_at_z.ok_or(SynthesisError::AssignmentMissing)?,
                };

                let beta_for_lookup_permutation = beta_for_lookup.ok_or(SynthesisError::AssignmentMissing)?;
                let gamma_for_lookup_permutation = gamma_for_lookup.ok_or(SynthesisError::AssignmentMissing)?;
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

                lookup_query = Some(lookup_queries);
            }
        }

        let mut lhs = quotient_at_z;
        lhs.mul_assign(&vanishing_at_z);

        let rhs = t_num_on_full_domain;

        if lhs != rhs {
            return Ok(((E::G1Affine::zero(), E::G1Affine::zero()), false));
        }
    }

    // now we need to reconstruct the effective linearization poly with homomorphic properties
    let linearization_commitment = {

        let mut challenges_slice = &powers_of_alpha_for_gates[..];

        let mut all_gates = sorted_gates.clone();

        let gate = all_gates.drain(0..1).into_iter().next().ok_or(SynthesisError::AssignmentMissing)?;
        safe_assert(gate.benefits_from_linearization())?;
        safe_assert(C::MainGate::default().into_internal() == gate)?;
        let gate = C::MainGate::default();
        let num_challenges = gate.num_quotient_terms();
        let (for_gate, rest) = challenges_slice.split_at(num_challenges);
        challenges_slice = rest;

        let input_values = proof.inputs.clone();

        let mut r = gate.contribute_into_linearization_commitment_for_public_inputs(
            required_domain_size,
            &input_values,
            z,
            &query_values_map,
            &setup_commitments_storage,
            for_gate,
        )?;

        let mut selectors_it = selector_values.clone().into_iter();

        if num_different_gates > 1 {
            // first multiply r by the selector value at z
            r.mul_assign(selectors_it.next().ok_or(SynthesisError::AssignmentMissing)?.into_repr());
        }

        // now proceed per gate
        for gate in all_gates.into_iter() {
            let num_challenges = gate.num_quotient_terms();
            let (for_gate, rest) = challenges_slice.split_at(num_challenges);
            challenges_slice = rest;

            if gate.benefits_from_linearization() {
                // gate benefits from linearization, so make temporary value
                let tmp = gate.contribute_into_linearization_commitment(
                    required_domain_size,
                    z,
                    &query_values_map,
                    &setup_commitments_storage,
                    for_gate,
                )?;

                let selector_value = selectors_it.next().ok_or(SynthesisError::AssignmentMissing)?;
                let mut scaled = tmp;
                scaled.mul_assign(selector_value.into_repr());

                r.add_assign(&scaled);
            } else {
                // we linearize over the selector, so take a selector and scale it
                let gate_value_at_z = gate.contribute_into_verification_equation(
                    required_domain_size,
                    z,
                    &query_values_map,
                    for_gate
                )?;

                let key = PolyIdentifier::GateSelector(gate.name());
                let gate_selector = gate_selectors_commitments_storage.get(&key).ok_or(SynthesisError::AssignmentMissing)?;
                let scaled = gate_selector.mul(gate_value_at_z.into_repr());

                r.add_assign(&scaled);
            }
        }

        safe_assert(selectors_it.next().is_none())?;
        safe_assert_eq(challenges_slice.len(), 0)?;

        // add contributions from copy-permutation and lookup-permutation

        // copy-permutation linearization comtribution
        {
            // + (a(z) + beta*z + gamma)*()*()*()*Z(x)

            let [alpha_0, alpha_1] = copy_grand_product_alphas.expect("there must be powers of alpha for copy permutation");

            let some_one = Some(E::Fr::one());
            let mut non_residues_iterator = some_one.iter().chain(&vk.non_residues);

            let mut factor = alpha_0;

            for idx in 0..vk.state_width {
                let key = PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(idx));
                let wire_value = query_values_map.get(&key)
                    .ok_or(SynthesisError::AssignmentMissing)?;
                let mut t = z;
                let non_res = non_residues_iterator.next().ok_or(SynthesisError::AssignmentMissing)?;
                t.mul_assign(&non_res);
                t.mul_assign(&beta_for_copy_permutation);
                t.add_assign(&wire_value);
                t.add_assign(&gamma_for_copy_permutation);

                factor.mul_assign(&t);
            }

            safe_assert(non_residues_iterator.next().is_none())?;

            let scaled = proof.copy_permutation_grand_product_commitment.mul(factor.into_repr());

            r.add_assign(&scaled);

            // - (a(z) + beta*perm_a + gamma)*()*()*z(z*omega) * beta * perm_d(X)

            let mut factor = alpha_0;
            factor.mul_assign(&beta_for_copy_permutation);
            factor.mul_assign(&copy_permutation_z_at_z_omega);

            for idx in 0..(vk.state_width-1) {
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

            let scaled = vk.permutation_commitments.get(vk.state_width - 1).ok_or(SynthesisError::AssignmentMissing)?.mul(factor.into_repr());
            
            r.sub_assign(&scaled);

            // + L_0(z) * Z(x)

            let mut factor = evaluate_l0_at_point(required_domain_size as u64, z)?;
            factor.mul_assign(&alpha_1);

            let scaled = proof.copy_permutation_grand_product_commitment.mul(factor.into_repr());

            r.add_assign(&scaled);
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

        if vk.total_lookup_entries_length > 0 {
            let [alpha_0, alpha_1, alpha_2] = lookup_grand_product_alphas.expect("there must be powers of alpha for lookup permutation");

            let lookup_queries = lookup_query.expect("lookup data must be constructed");

            // let s_at_z_omega = lookup_queries.s_at_z_omega;
            let grand_product_at_z_omega = lookup_queries.grand_product_at_z_omega;
            let t_at_z = lookup_queries.t_at_z;
            let t_at_z_omega = lookup_queries.t_at_z_omega;
            let selector_at_z = lookup_queries.selector_at_z;
            let table_type_at_z = lookup_queries.table_type_at_z;

            let l_0_at_z = evaluate_lagrange_poly_at_point(0, &domain, z)?;
            let l_n_minus_one_at_z = evaluate_lagrange_poly_at_point(required_domain_size - 1, &domain, z)?;

            let beta_for_lookup_permutation = beta_for_lookup.ok_or(SynthesisError::AssignmentMissing)?;
            let gamma_for_lookup_permutation = gamma_for_lookup.ok_or(SynthesisError::AssignmentMissing)?;

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

            let scaled = proof.lookup_s_poly_commitment.ok_or(SynthesisError::AssignmentMissing)?.mul(factor.into_repr());

            r.add_assign(&scaled);

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
            let eta = eta;
            // a,b,c
            safe_assert_eq(vk.state_width, 4)?;
            for idx in 0..(vk.state_width-1) {
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

            let scaled = proof.lookup_grand_product_commitment.ok_or(SynthesisError::AssignmentMissing)?.mul(factor.into_repr());

            r.add_assign(&scaled);
        }

        r.into_affine()
    };

    let v = transcript.get_challenge();

    // commit proofs

    commit_point_as_xy::<E, T>(&mut transcript, &proof.opening_proof_at_z);
    commit_point_as_xy::<E, T>(&mut transcript, &proof.opening_proof_at_z_omega);

    let u = transcript.get_challenge();

    // first perform naive verification at z
    // f(x) - f(z) = q(x)(x - z) =>
    // e(f(x) - f(z)*g + z*q(x), h)*e(-q(x), h^x) == 1
    // when we aggregate we need to aggregate f(x) part (commitments) and f(z) part (values)

    let mut values_queried_at_z = vec![quotient_at_z];
    values_queried_at_z.push(linearization_at_z);
    values_queried_at_z.extend(all_values_queried_at_z);

    let quotient_commitment_aggregated = {
        let mut quotient_commitments_iter = proof.quotient_poly_parts_commitments.iter();
        let mut result = quotient_commitments_iter.next().ok_or(SynthesisError::AssignmentMissing)?.into_projective();

        let mut current = z_in_domain_size;
        for part in quotient_commitments_iter {
            let tmp = *part;
            let tmp = tmp.mul(current.into_repr());

            result.add_assign(&tmp);
            current.mul_assign(&z_in_domain_size);
        }

        result.into_affine()
    };

    let mut commitments_queried_at_z = vec![];
    commitments_queried_at_z.push(quotient_commitment_aggregated);
    commitments_queried_at_z.push(linearization_commitment);
    commitments_queried_at_z.extend(all_commitments_queried_at_z);

    let mut reconstructed_lookup_t_poly_commitment = None;

    if vk.total_lookup_entries_length > 0 {
        // we need to add t(x), selector(x) and table type(x)
        values_queried_at_z.push(proof.lookup_t_poly_opening_at_z.ok_or(SynthesisError::AssignmentMissing)?);
        values_queried_at_z.push(proof.lookup_selector_poly_opening_at_z.ok_or(SynthesisError::AssignmentMissing)?);
        values_queried_at_z.push(proof.lookup_table_type_poly_opening_at_z.ok_or(SynthesisError::AssignmentMissing)?);

        // use eta to reconstruct t poly aggregated commitment
        let lookup_t_poly_commitment_aggregated = {
            let mut commitments_iter = vk.lookup_tables_commitments.iter();
            let mut result = commitments_iter.next().ok_or(SynthesisError::AssignmentMissing)?.into_projective();
    
            let mut current = eta;
            for part in commitments_iter {
                let tmp = *part;
                let tmp = tmp.mul(current.into_repr());
    
                result.add_assign(&tmp);
                current.mul_assign(&eta);
            }
    
            result.into_affine()
        };

        reconstructed_lookup_t_poly_commitment = Some(lookup_t_poly_commitment_aggregated);

        commitments_queried_at_z.push(lookup_t_poly_commitment_aggregated);
        commitments_queried_at_z.push(vk.lookup_selector_commitment.ok_or(SynthesisError::AssignmentMissing)?);
        commitments_queried_at_z.push(vk.lookup_table_type_commitment.ok_or(SynthesisError::AssignmentMissing)?);
    }

    let mut values_queried_at_z_omega = vec![copy_permutation_z_at_z_omega];
    values_queried_at_z_omega.extend(all_values_queried_at_z_omega);

    let mut commitments_queried_at_z_omega = vec![proof.copy_permutation_grand_product_commitment];
    commitments_queried_at_z_omega.extend(all_commitments_queried_at_z_omega);

    if vk.total_lookup_entries_length > 0 {
        // we need to add s(x), grand_product(x) and t(x)
        values_queried_at_z_omega.push(proof.lookup_s_poly_opening_at_z_omega.ok_or(SynthesisError::AssignmentMissing)?);
        values_queried_at_z_omega.push(proof.lookup_grand_product_opening_at_z_omega.ok_or(SynthesisError::AssignmentMissing)?);
        values_queried_at_z_omega.push(proof.lookup_t_poly_opening_at_z_omega.ok_or(SynthesisError::AssignmentMissing)?);

        commitments_queried_at_z_omega.push(proof.lookup_s_poly_commitment.ok_or(SynthesisError::AssignmentMissing)?);
        commitments_queried_at_z_omega.push(proof.lookup_grand_product_commitment.ok_or(SynthesisError::AssignmentMissing)?);
        commitments_queried_at_z_omega.push(reconstructed_lookup_t_poly_commitment.expect("t poly for lookup must be reconstructed"));
    }

    safe_assert_eq(commitments_queried_at_z.len(), values_queried_at_z.len())?;
    safe_assert_eq(commitments_queried_at_z_omega.len(), values_queried_at_z_omega.len())?;

    let mut aggregated_commitment_at_z = commitments_queried_at_z.drain(0..1).next().ok_or(SynthesisError::AssignmentMissing)?.into_projective();
    let mut aggregated_opening_at_z = values_queried_at_z.drain(0..1).next().ok_or(SynthesisError::AssignmentMissing)?;

    let mut aggregation_challenge = E::Fr::one();

    for (commitment, value) in commitments_queried_at_z.into_iter().zip(values_queried_at_z.into_iter()) {
        aggregation_challenge.mul_assign(&v);

        let scaled = commitment.mul(aggregation_challenge.into_repr());
        aggregated_commitment_at_z.add_assign(&scaled);
        // dbg!(aggregated_commitment_at_z.into_affine());

        let mut tmp = value;
        tmp.mul_assign(&aggregation_challenge);
        aggregated_opening_at_z.add_assign(&tmp);
    }

    aggregation_challenge.mul_assign(&v);

    let mut aggregated_commitment_at_z_omega = commitments_queried_at_z_omega.drain(0..1).next().ok_or(SynthesisError::AssignmentMissing)?.mul(aggregation_challenge.into_repr());
    let mut aggregated_opening_at_z_omega = values_queried_at_z_omega.drain(0..1).next().ok_or(SynthesisError::AssignmentMissing)?;
    aggregated_opening_at_z_omega.mul_assign(&aggregation_challenge);

    for (commitment, value) in commitments_queried_at_z_omega.into_iter().zip(values_queried_at_z_omega.into_iter()) {
        aggregation_challenge.mul_assign(&v);

        let scaled = commitment.mul(aggregation_challenge.into_repr());
        aggregated_commitment_at_z_omega.add_assign(&scaled);

        let mut tmp = value;
        tmp.mul_assign(&aggregation_challenge);
        aggregated_opening_at_z_omega.add_assign(&tmp);
    }

    // f(x)
    let mut pair_with_generator = aggregated_commitment_at_z;
    aggregated_commitment_at_z_omega.mul_assign(u.into_repr());
    pair_with_generator.add_assign(&aggregated_commitment_at_z_omega);

    // - f(z)*g
    let mut aggregated_value = aggregated_opening_at_z_omega;
    aggregated_value.mul_assign(&u);
    aggregated_value.add_assign(&aggregated_opening_at_z);

    let tmp = E::G1Affine::one().mul(aggregated_value.into_repr());
    pair_with_generator.sub_assign(&tmp);

    // +z * q(x)
    let mut tmp = proof.opening_proof_at_z.mul(z.into_repr());

    let mut t0 = z_omega;
    t0.mul_assign(&u);
    let t1 = proof.opening_proof_at_z_omega.mul(t0.into_repr());
    tmp.add_assign(&t1);
    pair_with_generator.add_assign(&tmp);

    // rhs
    let mut pair_with_x = proof.opening_proof_at_z_omega.mul(u.into_repr());
    pair_with_x.add_assign_mixed(&proof.opening_proof_at_z);

    let mut pair_with_x = pair_with_x.into_affine();
    pair_with_x.negate();

    let pair_with_generator = pair_with_generator.into_affine();
    
    Ok(((pair_with_generator, pair_with_x), true))
}
