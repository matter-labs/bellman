use crate::pairing::{
    Engine,
    CurveProjective,
    CurveAffine,
};

use crate::pairing::ff::{Field, PrimeField};

use super::{
    Proof,
    VerifyingKey,
    PreparedVerifyingKey,
};

use crate::{
    SynthesisError
};

pub fn prepare_verifying_key<E: Engine>(
    vk: &VerifyingKey<E>
) -> PreparedVerifyingKey<E>
{
    PreparedVerifyingKey {
        h_pc: vk.h_g2.prepare(),
        g_alpha: vk.alpha_g1.clone(),
        h_beta: vk.beta_g2.clone(),
        g_alpha_h_beta_ml: E::pairing(vk.alpha_g1, vk.beta_g2),
        g_gamma_pc: vk.gamma_g1.prepare(),
        h_gamma_pc: vk.gamma_g2.prepare(),
        query: vk.ic.clone(),
    }
}

pub fn verify_proof<'a, E: Engine>(
    pvk: &'a PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    public_inputs: &[E::Fr]
) -> Result<bool, SynthesisError>
{
    // checks input length against query length
    if (public_inputs.len() + 1) != pvk.query.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    // e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma}) *
    // e(C, H) where psi = \sum_{i=0}^l input_i pvk.query[i]

    let mut g_psi = pvk.query[0].into_projective();
    for (i, b) in public_inputs.iter().zip(pvk.query.iter().skip(1)) {
        g_psi.add_assign(&b.mul(i.into_repr()));
    }

    let mut test1_a_g_alpha = proof.a.into_projective();
    test1_a_g_alpha.add_assign_mixed(&pvk.g_alpha);
    let mut test1_a_g_alpha = test1_a_g_alpha.into_affine();
    test1_a_g_alpha.negate();

    let mut test1_b_h_beta = proof.b.into_projective();
    test1_b_h_beta.add_assign_mixed(&pvk.h_beta);
    let test1_b_h_beta = test1_b_h_beta.into_affine();


    let test1_r1 = pvk.g_alpha_h_beta_ml;
    let mut test1_r2 = E::miller_loop([
        (&test1_a_g_alpha.prepare(), &test1_b_h_beta.prepare()),
        (&g_psi.into_affine().prepare(), &pvk.h_gamma_pc),
        (&proof.c.prepare(), &pvk.h_pc),
    ].iter());
    test1_r2.mul_assign(&test1_r1);
    let test1 = E::final_exponentiation(&test1_r2).unwrap();

    // e(A, H^{gamma}) = e(G^{gamma}, B)
    let mut neg_b = proof.b;
    neg_b.negate();
    let test2 = E::final_exponentiation(&E::miller_loop([
        (&proof.a.prepare(), &pvk.h_gamma_pc),
        (&pvk.g_gamma_pc, &neg_b.prepare()),
    ].iter())).unwrap();

    let one = E::Fqk::one();
    Ok(test1 == one && test2 == one)
}
