use std::sync::{Arc, Mutex};

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField, UniformRand, Zero};
use ark_groth16::{ProvingKey, VerifyingKey};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use ark_std::{cfg_into_iter, cfg_iter};
use rand::Rng;

use crate::{
    error::Error,
    key::{FullKey, PartialKey},
    keypair::PublicKey,
    pot::Accumulator,
    ratio::RatioProof,
    utils::{
        batch_into_projective, batch_mul_fixed_scalar, merge_ratio_affine_vec, same_ratio_swap,
        seeded_rng,
    },
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, PartialEq)]
pub struct Transcript<E: PairingEngine> {
    pub key: FullKey<E>,
    pub initial_key: PartialKey<E>,
    pub contributions: Vec<PublicKey<E>>,
}

impl<E: PairingEngine> Transcript<E> {
    pub fn new_from_accumulator<C: ConstraintSynthesizer<E::Fr>>(
        accum: &Accumulator<E>,
        circuit: C,
    ) -> Result<Self, Error> {
        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone())?;
        cs.finalize();

        let num_constraints = cs.num_constraints();
        let num_instance_variables = cs.num_instance_variables();
        let total_constraints = num_constraints + num_instance_variables;
        let constraint_matrices = cs.to_matrices().ok_or(Error::MissingCSMatrices)?;

        let domain = Radix2EvaluationDomain::new(total_constraints)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let degree = domain.size();

        let (valid, g1_len, g2_len) = accum.check_pow_len();
        valid.then_some(()).ok_or(Error::InvalidPPOTPowers)?;
        (g2_len >= total_constraints && g1_len >= (degree << 1))
            .then_some(())
            .ok_or(Error::NotEnoughTauPowers(degree + 1))?;

        let h_query = cfg_into_iter!(0..degree)
            .map(|i| {
                let tmp = accum.tau_powers_g1[i + degree].into_projective();
                let tmp2 = accum.tau_powers_g1[i].into_projective();
                (tmp - tmp2).into_affine()
            })
            .collect::<Vec<_>>();

        let tau_lagrange_g1 = domain.ifft(&batch_into_projective(&accum.tau_powers_g1));
        let tau_lagrange_g2 = domain.ifft(&batch_into_projective(&accum.tau_powers_g2));
        let alpha_lagrange_g1 = domain.ifft(&batch_into_projective(&accum.alpha_tau_powers_g1));
        let beta_lagrange_g1 = domain.ifft(&batch_into_projective(&accum.beta_tau_powers_g1));
        let num_witnesses =
            constraint_matrices.num_witness_variables + constraint_matrices.num_instance_variables;
        let a_g1 = Arc::new(Mutex::new(vec![E::G1Projective::zero(); num_witnesses]));
        let b_g1 = Arc::new(Mutex::new(vec![E::G1Projective::zero(); num_witnesses]));
        let b_g2 = Arc::new(Mutex::new(vec![E::G2Projective::zero(); num_witnesses]));
        let ext = Arc::new(Mutex::new(vec![E::G1Projective::zero(); num_witnesses]));

        a_g1.lock()?[0..num_instance_variables]
            .clone_from_slice(&tau_lagrange_g1[num_constraints..total_constraints]);
        ext.lock()?[0..num_instance_variables]
            .clone_from_slice(&beta_lagrange_g1[num_constraints..total_constraints]);

        cfg_iter!(constraint_matrices.a)
            .zip(constraint_matrices.b)
            .zip(constraint_matrices.c)
            .zip(tau_lagrange_g1)
            .zip(tau_lagrange_g2)
            .zip(alpha_lagrange_g1)
            .zip(beta_lagrange_g1)
            .for_each(
                |((((((a_poly, b_poly), c_poly), tau_g1), tau_g2), alpha_tau), beta_tau)| {
                    cfg_iter!(a_poly).for_each(|(coeff, index)| {
                        a_g1.lock().unwrap()[*index] += tau_g1.mul(coeff.into_repr());
                        ext.lock().unwrap()[*index] += beta_tau.mul(coeff.into_repr());
                    });
                    cfg_iter!(b_poly).for_each(|(coeff, index)| {
                        b_g1.lock().unwrap()[*index] += tau_g1.mul(coeff.into_repr());
                        b_g2.lock().unwrap()[*index] += tau_g2.mul(coeff.into_repr());
                        ext.lock().unwrap()[*index] += alpha_tau.mul(coeff.into_repr());
                    });
                    cfg_iter!(c_poly).for_each(|(coeff, index)| {
                        ext.lock().unwrap()[*index] += tau_g1.mul(coeff.into_repr());
                    });
                },
            );

        let a_query = ProjectiveCurve::batch_normalization_into_affine(&a_g1.lock()?);
        let b_g1_query = ProjectiveCurve::batch_normalization_into_affine(&b_g1.lock()?);
        let b_g2_query = ProjectiveCurve::batch_normalization_into_affine(&b_g2.lock()?);
        let ext = ProjectiveCurve::batch_normalization_into_affine(&ext.lock()?);

        let public_cross_terms = Vec::from(&ext[..constraint_matrices.num_instance_variables]);
        let private_cross_terms = Vec::from(&ext[constraint_matrices.num_instance_variables..]);

        let key = ProvingKey::<E> {
            vk: VerifyingKey {
                alpha_g1: accum.alpha_tau_powers_g1[0],
                beta_g2: accum.beta_g2,
                gamma_g2: E::G2Affine::prime_subgroup_generator(),
                delta_g2: E::G2Affine::prime_subgroup_generator(),
                gamma_abc_g1: public_cross_terms,
            },
            beta_g1: accum.beta_tau_powers_g1[0],
            delta_g1: E::G1Affine::prime_subgroup_generator(),
            a_query,
            b_g1_query,
            b_g2_query,
            h_query: h_query.to_vec(),
            l_query: private_cross_terms,
        };

        Ok(Self {
            initial_key: (&key).into(),
            key: FullKey { key },
            contributions: vec![],
        })
    }

    pub fn contribute_seed(&mut self, seed: &[u8]) -> Result<(), Error> {
        self.contribute_rng(&mut seeded_rng(seed))
    }

    pub fn contribute_rng<R: Rng>(&mut self, rng: &mut R) -> Result<(), Error> {
        let delta = E::Fr::rand(rng);
        let delta_inverse = delta.inverse().expect("delta is not invertible");

        let proof = RatioProof::<E>::generate(delta, &self.key.challenge()?)?;

        batch_mul_fixed_scalar(&mut self.key.key.l_query, delta_inverse);
        batch_mul_fixed_scalar(&mut self.key.key.h_query, delta_inverse);
        self.key.key.delta_g1 = self.key.key.delta_g1.mul(delta).into_affine();
        self.key.key.vk.delta_g2 = self.key.key.vk.delta_g2.mul(delta).into_affine();
        self.contributions.push(PublicKey {
            delta_g2: self.key.key.vk.delta_g2,
            proof,
        });

        Ok(())
    }

    #[inline]
    pub fn verify(&self) -> Result<(), Error> {
        let mut challenge: (E::G2Affine, Vec<u8>) =
            (self.initial_key.delta_g2, self.initial_key.challenge()?);
        for contribution in self.contributions.iter() {
            contribution
                .proof
                .verify(&challenge.1)
                .map_err(|_| Error::InvalidRatioProof)?;

            same_ratio_swap::<E>(
                contribution.proof.get_g1(),
                (challenge.0, contribution.delta_g2),
            )
            .then_some(())
            .ok_or(Error::InvalidRatioProof)?;
            challenge = (contribution.delta_g2, contribution.challenge()?);
        }

        same_ratio_swap::<E>(
            (self.initial_key.delta_g1, self.key.key.delta_g1),
            (self.initial_key.delta_g2, challenge.0),
        )
        .then_some(())
        .ok_or(Error::InvalidRatioProof)?;

        same_ratio_swap::<E>(
            merge_ratio_affine_vec(&self.key.key.h_query, &self.initial_key.h_query),
            (self.initial_key.delta_g2, self.key.key.vk.delta_g2),
        )
        .then_some(())
        .ok_or(Error::InconsistentHChange)?;

        same_ratio_swap::<E>(
            merge_ratio_affine_vec(&self.key.key.l_query, &self.initial_key.l_query),
            (self.initial_key.delta_g2, self.key.key.vk.delta_g2),
        )
        .then_some(())
        .ok_or(Error::InconsistentLChange)?;
        Ok(())
    }

    #[inline]
    pub fn verify_key_transform(
        prev: &PartialKey<E>,
        next: &PartialKey<E>,
        proof: &RatioProof<E>,
    ) -> Result<(), Error> {
        proof
            .verify(&prev.challenge()?)
            .map_err(|_| Error::InvalidRatioProof)?;

        (same_ratio_swap::<E>(proof.get_g1(), (prev.delta_g2, next.delta_g2))
            && same_ratio_swap::<E>(
                (prev.delta_g1, next.delta_g1),
                (prev.delta_g2, next.delta_g2),
            ))
        .then_some(())
        .ok_or(Error::InvalidRatioProof)?;

        same_ratio_swap::<E>(
            merge_ratio_affine_vec(&next.h_query, &prev.h_query),
            (prev.delta_g2, next.delta_g2),
        )
        .then_some(())
        .ok_or(Error::InconsistentHChange)?;

        same_ratio_swap::<E>(
            merge_ratio_affine_vec(&next.l_query, &prev.l_query),
            (prev.delta_g2, next.delta_g2),
        )
        .then_some(())
        .ok_or(Error::InconsistentLChange)?;

        Ok(())
    }

    #[inline]
    pub fn verify_from_accumulator<C: ConstraintSynthesizer<E::Fr>>(
        &self,
        accum: &Accumulator<E>,
        circuit: C,
    ) -> Result<(), Error> {
        let initial_transcript = Transcript::new_from_accumulator(accum, circuit)?;
        (initial_transcript.initial_key == self.initial_key)
            .then_some(())
            .ok_or(Error::InvalidKey("initial_key"))?;
        (initial_transcript.key.key.beta_g1 == self.key.key.beta_g1)
            .then_some(())
            .ok_or(Error::InvalidKey("pk: beta_g1"))?;
        (initial_transcript.key.key.a_query == self.key.key.a_query)
            .then_some(())
            .ok_or(Error::InvalidKey("pk: a_query"))?;
        (initial_transcript.key.key.b_g1_query == self.key.key.b_g1_query)
            .then_some(())
            .ok_or(Error::InvalidKey("pk: b_g1"))?;
        (initial_transcript.key.key.b_g2_query == self.key.key.b_g2_query)
            .then_some(())
            .ok_or(Error::InvalidKey("pk: b_g2"))?;
        (initial_transcript.key.key.vk.alpha_g1 == self.key.key.vk.alpha_g1)
            .then_some(())
            .ok_or(Error::InvalidKey("vk: alpha_g1"))?;
        (initial_transcript.key.key.vk.beta_g2 == self.key.key.vk.beta_g2)
            .then_some(())
            .ok_or(Error::InvalidKey("vk: beta_g2"))?;
        (initial_transcript.key.key.vk.gamma_g2 == self.key.key.vk.gamma_g2)
            .then_some(())
            .ok_or(Error::InvalidKey("vk: gamma_g2"))?;
        (initial_transcript.key.key.vk.gamma_abc_g1 == self.key.key.vk.gamma_abc_g1)
            .then_some(())
            .ok_or(Error::InvalidKey("vk: gamma_abc_g1"))?;

        self.verify()?;

        Ok(())
    }
}
