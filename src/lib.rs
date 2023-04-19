pub mod error;
pub mod hasher;
pub mod key;
pub mod keypair;
pub mod pot;
pub mod ratio;
pub mod utils;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField, UniformRand, Zero};
use ark_groth16::{ProvingKey, VerifyingKey};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};
use ark_std::{cfg_into_iter, cfg_iter};
use key::{FullKey, PartialKey};
use keypair::PublicKey;
use pot::Accumulator;
use rand::Rng;
use ratio::RatioProof;
use utils::{merge_ratio_affine_vec, same_ratio_swap, seeded_rng};

use crate::{
    error::Error,
    utils::{batch_into_projective, batch_mul_fixed_scalar},
};

#[derive(Debug, Clone, PartialEq)]
pub struct Transcript<E: PairingEngine> {
    pub key: FullKey<E>,
    pub initial_key: PartialKey<E>,
    pub contributions: Vec<PublicKey<E>>,
}

impl<E: PairingEngine> Transcript<E> {
    pub fn new_from_accumulator<C: ConstraintSynthesizer<E::Fr>>(
        accum: Accumulator<E>,
        circuit: C,
    ) -> Result<Self, Error> {
        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone())?;
        cs.finalize();

        let num_constraints = cs.num_constraints();
        let num_instance_variables = cs.num_instance_variables();
        let constraint_matrices = cs.to_matrices().ok_or(Error::MissingCSMatrices)?;

        let domain = Radix2EvaluationDomain::new(num_constraints + num_instance_variables)
            .ok_or(Error::TooManyConstraints)?;
        let degree = domain.size as usize;

        let mut h_query = vec![E::G1Affine::zero(); degree];
        cfg_into_iter!(0..degree).for_each(|i| {
            let tmp = accum.tau_powers_g1[i + degree].into_projective();
            let tmp2 = accum.tau_powers_g1[i].into_projective();
            h_query[i] = (tmp - tmp2).into_affine();
        });

        let beta_g1 = accum.beta_tau_powers_g1[0];
        let tau_lagrange_g1 = domain.ifft(&batch_into_projective(&accum.tau_powers_g1));
        let tau_lagrange_g2 = domain.ifft(&batch_into_projective(&accum.tau_powers_g2));
        let alpha_lagrange_g1 = domain.ifft(&batch_into_projective(&accum.alpha_tau_powers_g1));
        let beta_lagrange_g1 = domain.ifft(&batch_into_projective(&accum.beta_tau_powers_g1));
        let num_witnesses =
            constraint_matrices.num_witness_variables + constraint_matrices.num_instance_variables;
        let mut a_g1 = vec![E::G1Projective::zero(); num_witnesses];
        let mut b_g1 = vec![E::G1Projective::zero(); num_witnesses];
        let mut b_g2 = vec![E::G2Projective::zero(); num_witnesses];
        let mut ext = vec![E::G1Projective::zero(); num_witnesses];

        let total = num_instance_variables + num_constraints;
        a_g1[0..num_instance_variables].clone_from_slice(&tau_lagrange_g1[num_constraints..total]);
        ext[0..num_instance_variables].clone_from_slice(&beta_lagrange_g1[num_constraints..total]);

        assert_eq!(a_g1.len(), b_g1.len());
        assert_eq!(a_g1.len(), b_g2.len());
        assert_eq!(a_g1.len(), ext.len());

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
                        a_g1[*index] += tau_g1.mul(coeff.into_repr());
                        ext[*index] += beta_tau.mul(coeff.into_repr());
                    });
                    cfg_iter!(b_poly).for_each(|(coeff, index)| {
                        b_g1[*index] += tau_g1.mul(coeff.into_repr());
                        b_g2[*index] += tau_g2.mul(coeff.into_repr());
                        ext[*index] += alpha_tau.mul(coeff.into_repr());
                    });
                    cfg_iter!(c_poly).for_each(|(coeff, index)| {
                        ext[*index] += tau_g1.mul(coeff.into_repr());
                    });
                },
            );

        let a_query = ProjectiveCurve::batch_normalization_into_affine(&a_g1);
        let b_g1_query = ProjectiveCurve::batch_normalization_into_affine(&b_g1);
        let b_g2_query = ProjectiveCurve::batch_normalization_into_affine(&b_g2);
        let ext = ProjectiveCurve::batch_normalization_into_affine(&ext);
        let public_cross_terms = Vec::from(&ext[..constraint_matrices.num_instance_variables]);
        let private_cross_terms = Vec::from(&ext[constraint_matrices.num_instance_variables..]);
        let key: FullKey<E> = ProvingKey::<E> {
            vk: VerifyingKey {
                alpha_g1: accum.alpha_tau_powers_g1[0],
                beta_g2: accum.beta_g2,
                gamma_g2: E::G2Affine::prime_subgroup_generator(),
                delta_g2: E::G2Affine::prime_subgroup_generator(),
                gamma_abc_g1: public_cross_terms,
            },
            beta_g1,
            delta_g1: E::G1Affine::prime_subgroup_generator(),
            a_query,
            b_g1_query,
            b_g2_query,
            h_query,
            l_query: private_cross_terms,
        }
        .into();

        Ok(Self {
            initial_key: (&key).into(),
            key,
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
            partial_key: (&self.key).into(),
            proof,
        });

        Ok(())
    }

    #[inline]
    pub fn verify(&self) -> Result<(), Error> {
        let mut key: &PartialKey<E> = &self.initial_key;
        for contribution in self.contributions.iter() {
            contribution
                .proof
                .verify(&key.challenge()?)
                .map_err(|_| Error::InvalidRatioProof)?;

            same_ratio_swap::<E>(
                contribution.proof.get_g1(),
                (key.delta_g2, contribution.partial_key.delta_g2),
            )
            .then_some(())
            .ok_or(Error::InvalidRatioProof)?;
            key = &contribution.partial_key;
        }

        same_ratio_swap::<E>(
            (self.initial_key.delta_g1, key.delta_g1),
            (self.initial_key.delta_g2, key.delta_g2),
        )
        .then_some(())
        .ok_or(Error::InvalidRatioProof)?;

        same_ratio_swap::<E>(
            merge_ratio_affine_vec(&key.h_query, &self.initial_key.h_query),
            (self.initial_key.delta_g2, key.delta_g2),
        )
        .then_some(())
        .ok_or(Error::InconsistentHChange)?;

        same_ratio_swap::<E>(
            merge_ratio_affine_vec(&key.l_query, &self.initial_key.l_query),
            (self.initial_key.delta_g2, key.delta_g2),
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
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use ark_bn254::{Bn254, Fr};
    use ark_ff::Zero;
    use ark_r1cs_std::{
        fields::fp::FpVar,
        prelude::{AllocVar, EqGadget},
    };
    use ark_relations::r1cs::ConstraintSynthesizer;

    use crate::{pot::Accumulator, utils::seeded_rng, Transcript};

    const NUM_CONSTRAINTS: usize = 10;

    struct DummyCircuit {
        pub a: Fr,
        pub b: Fr,
    }

    impl ConstraintSynthesizer<Fr> for DummyCircuit {
        fn generate_constraints(
            self,
            cs: ark_relations::r1cs::ConstraintSystemRef<Fr>,
        ) -> ark_relations::r1cs::Result<()> {
            let a = FpVar::new_witness(cs.clone(), || Ok(self.a))?;
            let b = FpVar::new_input(cs, || Ok(self.b))?;

            for _ in 0..NUM_CONSTRAINTS {
                a.enforce_equal(&b)?;
            }

            Ok(())
        }
    }

    #[test]
    fn groth16_domain() -> Result<(), Box<dyn Error>> {
        let accum = Accumulator::<Bn254>::generate(
            NUM_CONSTRAINTS.pow(2) + 1,
            1,
            &mut seeded_rng(b"hehehe"),
        );
        let mut transcript = Transcript::new_from_accumulator(
            accum,
            DummyCircuit {
                a: Fr::zero(),
                b: Fr::zero(),
            },
        )?;

        transcript.contribute_seed(b"ofekwoo")?;
        transcript.verify()?;

        transcript.contribute_seed(b"aewrog")?;
        transcript.verify()?;

        transcript.contribute_seed(b"aewrog")?;
        transcript.contribute_seed(b"agwi")?;
        transcript.contribute_seed(b"23op2j43gk3")?;
        transcript.contribute_seed(b"zdsf9vaobjk")?;
        transcript.contribute_seed(b"k12999")?;
        transcript.contribute_seed(b"193u9hxcv")?;
        transcript.contribute_seed(b"1200xcvl")?;
        transcript.contribute_seed(b"z0cvixkkk")?;
        transcript.verify()?;

        Ok(())
    }
}
