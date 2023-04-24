use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, Read, Seek, SeekFrom},
    iter,
    path::Path,
};

use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::{BigInteger, FpParameters, One, PrimeField, UniformRand};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_relations::r1cs::SynthesisError;
use ark_std::{add_to_trace, cfg_iter_mut, end_timer, start_timer};
use rand::Rng;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{error::Error, reader::PairingReader};

#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub struct Accumulator<E: PairingEngine> {
    pub tau_powers_g1: Vec<E::G1Affine>,
    pub tau_powers_g2: Vec<E::G2Affine>,
    pub alpha_tau_powers_g1: Vec<E::G1Affine>,
    pub beta_tau_powers_g1: Vec<E::G1Affine>,
    pub beta_g2: E::G2Affine,
}

impl<E: PairingEngine> Accumulator<E> {
    pub fn check_pow_len(&self) -> (bool, usize, usize) {
        let g1_len = self.tau_powers_g1.len();
        let g2_len = self.tau_powers_g2.len();

        (
            (g1_len << 1) >= g2_len
                && g2_len == self.alpha_tau_powers_g1.len()
                && g2_len == self.beta_tau_powers_g1.len(),
            g1_len,
            g2_len,
        )
    }

    pub fn contribute<R: Rng>(&mut self, rng: &mut R) {
        let tau = E::Fr::rand(rng);
        let alpha = E::Fr::rand(rng);
        let beta = E::Fr::rand(rng);

        let g1_powers = self.tau_powers_g1.len();
        let g2_powers = self.tau_powers_g2.len();

        let mut tau_powers = iter::successors(Some(E::Fr::one()), |x| Some(*x * tau))
            .take(g1_powers)
            .collect::<Vec<_>>();
        let remaining_tau_powers = tau_powers.split_off(g2_powers);

        cfg_iter_mut!(self.tau_powers_g1)
            .zip(cfg_iter_mut!(self.tau_powers_g2))
            .zip(cfg_iter_mut!(self.alpha_tau_powers_g1))
            .zip(cfg_iter_mut!(self.beta_tau_powers_g1))
            .zip(tau_powers)
            .for_each(
                |((((tau_g1, tau_g2), alpha_tau_g1), beta_tau_g1), tau_power)| {
                    *tau_g1 = tau_g1.mul(tau_power).into();
                    *tau_g2 = tau_g2.mul(tau_power).into();
                    *alpha_tau_g1 = alpha_tau_g1.mul(tau_power * alpha).into();
                    *beta_tau_g1 = beta_tau_g1.mul(tau_power * beta).into();
                },
            );
        cfg_iter_mut!(self.tau_powers_g1)
            .skip(g2_powers)
            .zip(remaining_tau_powers)
            .for_each(|(tau_g1, tau_power)| *tau_g1 = tau_g1.mul(tau_power).into());
        self.beta_g2 = self.beta_g2.mul(beta).into();
    }

    pub fn empty(g1_powers: usize, g2_powers: usize) -> Self {
        Self {
            tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g1_powers],
            tau_powers_g2: vec![E::G2Affine::prime_subgroup_generator(); g2_powers],
            alpha_tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g2_powers],
            beta_tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g2_powers],
            beta_g2: E::G2Affine::prime_subgroup_generator(),
        }
    }

    pub fn empty_from_max_constraints(max_constraints: usize) -> Result<Self, Error> {
        let domain = Radix2EvaluationDomain::<E::Fr>::new(max_constraints)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let degree = domain.size();
        let g2_powers = match ((max_constraints << 1) - 1) >= (degree << 1) {
            true => max_constraints,
            false => degree + 1,
        };
        Ok(Self::empty((g2_powers << 1) - 1, g2_powers))
    }
}

impl<E: PairingEngine + PairingReader> Accumulator<E> {
    pub fn from_ptau_file<P: AsRef<Path>>(path: P) -> Result<Accumulator<E>, Error> {
        let timer = start_timer!(|| "Reading from ptau file");

        let mut reader = BufReader::new(File::open(path)?);
        let section = {
            let timer = start_timer!(|| "Reading file header");

            let buf: &mut BufReader<File> = &mut reader;
            let mut b = [0u8; 4];
            buf.read_exact(&mut b)?;
            let ty = b.into_iter().map(char::from).collect::<String>();
            add_to_trace!(|| "Type: ", || ty);

            let mut b = [0u8; 4];
            buf.read_exact(&mut b)?;
            let version = u32::from_le_bytes(b);
            add_to_trace!(|| "Version: ", || version.to_string());

            let mut b = [0u8; 4];
            buf.read_exact(&mut b)?;
            let n_section = u32::from_le_bytes(b);
            add_to_trace!(|| "Sections : ", || n_section.to_string());

            let mut section = BTreeMap::new();
            for _ in 0..n_section {
                let mut b = [0u8; 4];
                buf.read_exact(&mut b)?;
                let ht = u32::from_le_bytes(b);

                let mut b = [0u8; 8];
                buf.read_exact(&mut b)?;
                let hl = u64::from_le_bytes(b);

                let current_pos = buf.stream_position()?;
                section
                    .entry(ht)
                    .or_insert(Vec::new())
                    .push((current_pos, hl));

                buf.seek_relative(hl.try_into()?)?;
            }

            end_timer!(timer);

            Result::<_, Error>::Ok(section)
        }?;

        let header = section
            .get(&1)
            .and_then(|e| e.first())
            .expect("Header not found");
        reader.seek(SeekFrom::Start(header.0))?;
        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let n = u32::from_le_bytes(b) as usize;

        let mut buffer = vec![0u8; n];
        reader.read_exact(&mut buffer)?;
        (buffer == <<E::Fq as PrimeField>::Params as FpParameters>::MODULUS.to_bytes_le())
            .then_some(())
            .ok_or(Error::PtauFieldNotMatching)?;

        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let power = u32::from_le_bytes(b);
        add_to_trace!(|| "Power: ", || power.to_string());

        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let ceremony_power = u32::from_le_bytes(b);
        add_to_trace!(|| "Ceremony Power: ", || ceremony_power.to_string());

        let ceremony = section
            .get(&7)
            .and_then(|e| e.first())
            .expect("Ceremony not found");
        reader.seek(SeekFrom::Start(ceremony.0))?;
        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let n_contribution = u32::from_le_bytes(b);
        add_to_trace!(|| "Ceremony Length: ", || n_contribution.to_string());

        let n_tau_g1s = 2u64.pow(power) * 2 - 1;
        let n_tau_g2s = 2u64.pow(power);

        let element_timer = start_timer!(|| "Reading pot elements");

        add_to_trace!(|| "Tau G1 Length: ", || n_tau_g1s.to_string());
        add_to_trace!(|| "Tau G2 Length: ", || n_tau_g2s.to_string());

        let tau_g1s = section
            .get(&2)
            .and_then(|e| e.first())
            .expect("tau_g1 not found");
        reader.seek(SeekFrom::Start(tau_g1s.0))?;

        let tau_g1_timer = start_timer!(|| "Reading tau g1");
        let tau_g1 = (0..n_tau_g1s)
            .map(|_| E::read_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(tau_g1_timer);

        let tau_g2_timer = start_timer!(|| "Reading tau g2");
        let tau_g2 = (0..n_tau_g2s)
            .map(|_| E::read_g2(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(tau_g2_timer);

        let alpha_tau_g1_timer = start_timer!(|| "Reading alpha tau g1");
        let alpha_tau_g1 = (0..n_tau_g2s)
            .map(|_| E::read_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(alpha_tau_g1_timer);

        let beta_tau_g1_timer = start_timer!(|| "Reading beta tau g1");
        let beta_tau_g1 = (0..n_tau_g2s)
            .map(|_| E::read_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(beta_tau_g1_timer);

        let beta_g2_timer = start_timer!(|| "Reading beta g2");
        let beta_g2 = E::read_g2(&mut reader)?;
        end_timer!(beta_g2_timer);

        end_timer!(element_timer);
        end_timer!(timer);

        Ok(Accumulator {
            tau_powers_g1: tau_g1,
            tau_powers_g2: tau_g2,
            alpha_tau_powers_g1: alpha_tau_g1,
            beta_tau_powers_g1: beta_tau_g1,
            beta_g2,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use ark_bn254::{Bn254, Fr};
    use ark_ff::Zero;
    use ark_groth16::Groth16;
    use ark_r1cs_std::{
        fields::fp::FpVar,
        prelude::{AllocVar, EqGadget},
    };
    use ark_relations::r1cs::ConstraintSynthesizer;
    use ark_snark::SNARK;
    use rand::rngs::OsRng;

    use crate::transcript::Transcript;

    use super::Accumulator;

    struct DummyCircuit {
        pub a: Fr,
        pub b: Fr,
        pub c: Fr,
    }

    impl ConstraintSynthesizer<Fr> for DummyCircuit {
        fn generate_constraints(
            self,
            cs: ark_relations::r1cs::ConstraintSystemRef<Fr>,
        ) -> ark_relations::r1cs::Result<()> {
            let a = FpVar::new_witness(cs.clone(), || Ok(self.a))?;
            let b = FpVar::new_witness(cs.clone(), || Ok(self.b))?;
            let c = FpVar::new_input(cs, || Ok(self.c))?;
            let d = &a + &b;

            c.enforce_equal(&d)?;

            Ok(())
        }
    }

    #[test]
    fn from_ptau_file_works() -> Result<(), Box<dyn Error>> {
        let rng = &mut OsRng;
        let ptau_path = "pot8.ptau";

        let accum = Accumulator::<Bn254>::from_ptau_file(ptau_path)?;
        let transcript = Transcript::new_from_accumulator(
            &accum,
            DummyCircuit {
                a: Fr::zero(),
                b: Fr::zero(),
                c: Fr::zero(),
            },
        )?;

        let proof = Groth16::prove(
            &transcript.key.key,
            DummyCircuit {
                a: Fr::from(1),
                b: Fr::from(1),
                c: Fr::from(2),
            },
            rng,
        )?;

        let valid = Groth16::verify(&transcript.key.key.vk, &[Fr::from(2)], &proof)?;
        assert!(valid, "Proof must be valid");

        let valid = Groth16::verify(&transcript.key.key.vk, &[Fr::from(4)], &proof)?;
        assert!(!valid, "Proof must be not valid");

        Ok(())
    }
}
