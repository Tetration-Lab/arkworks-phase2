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
use ark_std::{add_to_trace, cfg_into_iter, cfg_iter_mut, end_timer, start_timer};
use rand::Rng;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{
    error::Error,
    reader::PairingReader,
    utils::{batch_into_affine, batch_into_projective},
};

#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub struct Accumulator<E: PairingEngine> {
    pub tau_powers_g1: Vec<E::G1Affine>,
    pub tau_powers_g2: Vec<E::G2Affine>,
    pub alpha_tau_powers_g1: Vec<E::G1Affine>,
    pub beta_tau_powers_g1: Vec<E::G1Affine>,
    pub beta_g2: E::G2Affine,
}

impl<E: PairingEngine> Accumulator<E> {
    pub fn prepare_with_size(&self, size: usize) -> Result<PreparedAccumulator<E>, Error> {
        let timer = start_timer!(|| "Preparing accumulator");

        let (_, g1_len, g2_len) = self.check_pow_len();

        let domain =
            Radix2EvaluationDomain::new(size).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let size = domain.size();

        (g2_len >= size && g1_len >= (size << 1) - 1)
            .then_some(())
            .ok_or(Error::NotEnoughPOTDegree(ark_std::log2(size)))?;

        let h_query_timer = start_timer!(|| "Computing h_query");
        let h_query = cfg_into_iter!(0..size - 1)
            .map(|i| (self.tau_powers_g1[i + size] + -self.tau_powers_g1[i]))
            .collect::<Vec<_>>();
        end_timer!(h_query_timer);

        let tau_lagrange_g1_timer = start_timer!(|| "Computing inverse fft to tau_lagrange_g1");
        let tau_lagrange_g1 = domain.ifft(&batch_into_projective(&self.tau_powers_g1[..size]));
        end_timer!(tau_lagrange_g1_timer);

        let tau_lagrange_g2_timer = start_timer!(|| "Computing inverse fft to tau_lagrange_g2");
        let tau_lagrange_g2 = domain.ifft(&batch_into_projective(&self.tau_powers_g2[..size]));
        end_timer!(tau_lagrange_g2_timer);

        let alpha_lagrange_g1_timer = start_timer!(|| "Computing inverse fft to alpha_lagrange_g1");
        let alpha_lagrange_g1 =
            domain.ifft(&batch_into_projective(&self.alpha_tau_powers_g1[..size]));
        end_timer!(alpha_lagrange_g1_timer);

        let beta_lagrange_g1_timer = start_timer!(|| "Computing inverse fft to beta_lagrange_g1");
        let beta_lagrange_g1 =
            domain.ifft(&batch_into_projective(&self.beta_tau_powers_g1[..size]));
        end_timer!(beta_lagrange_g1_timer);

        end_timer!(timer);

        Ok(PreparedAccumulator {
            alpha: self.alpha_tau_powers_g1[0],
            beta: self.beta_tau_powers_g1[0],
            tau_lagrange_g1,
            tau_lagrange_g2,
            alpha_lagrange_g1,
            beta_lagrange_g1,
            beta_g2: self.beta_g2,
            h_query,
        })
    }

    pub fn prepare(&self) -> Result<PreparedAccumulator<E>, Error> {
        let (valid, _, g2_len) = self.check_pow_len();

        (valid && g2_len.is_power_of_two())
            .then_some(())
            .ok_or(Error::InvalidPOTSize)?;

        self.prepare_with_size(g2_len)
    }

    pub fn check_pow_len(&self) -> (bool, usize, usize) {
        let g1_len = self.tau_powers_g1.len();
        let g2_len = self.tau_powers_g2.len();

        (
            ((g1_len << 1) - 1) >= g2_len
                && g2_len == self.alpha_tau_powers_g1.len()
                && g2_len == self.beta_tau_powers_g1.len(),
            g1_len,
            g2_len,
        )
    }

    pub fn contribute<R: Rng>(&mut self, rng: &mut R) {
        let timer = start_timer!(|| "Contributing to accumulator");

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

        end_timer!(timer);
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

    pub fn empty_from_degree(degree: usize) -> Result<Self, Error> {
        (degree <= 28)
            .then_some(())
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

        Ok(Self::empty(((1 << degree) << 1) - 1, 1 << degree))
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

        let n_tau_g2s = 1 << power;
        let n_tau_g1s = (n_tau_g2s << 1) - 1;

        let element_timer = start_timer!(|| "Reading pot elements");

        add_to_trace!(|| "Tau G1 Length: ", || n_tau_g1s.to_string());
        add_to_trace!(|| "Tau G2 Length: ", || n_tau_g2s.to_string());

        let tau_g1_sec = section
            .get(&2)
            .and_then(|e| e.first())
            .expect("tau_g1 not found");
        reader.seek(SeekFrom::Start(tau_g1_sec.0))?;

        let tau_g1_timer = start_timer!(|| "Reading tau g1");
        let tau_g1 = (0..n_tau_g1s)
            .map(|_| E::read_ptau_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(tau_g1_timer);

        let tau_g2_sec = section
            .get(&3)
            .and_then(|e| e.first())
            .expect("tau_g2 not found");
        reader.seek(SeekFrom::Start(tau_g2_sec.0))?;

        let tau_g2_timer = start_timer!(|| "Reading tau g2");
        let tau_g2 = (0..n_tau_g2s)
            .map(|_| E::read_ptau_g2(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(tau_g2_timer);

        let alpha_tau_g1_sec = section
            .get(&4)
            .and_then(|e| e.first())
            .expect("alpha_tau_g1 not found");
        reader.seek(SeekFrom::Start(alpha_tau_g1_sec.0))?;

        let alpha_tau_g1_timer = start_timer!(|| "Reading alpha tau g1");
        let alpha_tau_g1 = (0..n_tau_g2s)
            .map(|_| E::read_ptau_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(alpha_tau_g1_timer);

        let beta_tau_g1_sec = section
            .get(&5)
            .and_then(|e| e.first())
            .expect("beta_tau_g1 not found");
        reader.seek(SeekFrom::Start(beta_tau_g1_sec.0))?;

        let beta_tau_g1_timer = start_timer!(|| "Reading beta tau g1");
        let beta_tau_g1 = (0..n_tau_g2s)
            .map(|_| E::read_ptau_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(beta_tau_g1_timer);

        let beta_g2_sec = section
            .get(&6)
            .and_then(|e| e.first())
            .expect("beta_g2 not found");
        reader.seek(SeekFrom::Start(beta_g2_sec.0))?;

        let beta_g2_timer = start_timer!(|| "Reading beta g2");
        let beta_g2 = E::read_ptau_g2(&mut reader)?;
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

    #[deprecated = "Only use this for lesser degree constraint generation"]
    pub fn from_radix_file<P: AsRef<Path>>(path: P) -> Result<Accumulator<E>, Error> {
        let timer = start_timer!(|| "Reading from phase1radix file");
        let path: &Path = path.as_ref();

        let n_taus = 1usize
            << path
                .iter()
                .last()
                .and_then(|p| {
                    p.to_str()
                        .and_then(|e| e.split("phase1radix2m").last())
                        .and_then(|e| e.parse::<usize>().ok())
                })
                .ok_or(crate::error::Error::Read(String::from(
                    "Unable to read evaluation domain from radix file name",
                )))?;
        let domain =
            Radix2EvaluationDomain::new(n_taus).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

        let mut reader = BufReader::new(File::open(path)?);

        let aux_timer = start_timer!(|| "Reading auxilary");
        let _a = E::read_radix_g1(&mut reader)?;
        let _b = E::read_radix_g1(&mut reader)?;
        end_timer!(aux_timer);

        let beta_g2_timer = start_timer!(|| "Reading beta g2");
        let beta_g2 = E::read_radix_g2(&mut reader)?;
        end_timer!(beta_g2_timer);

        let tau_g1_timer = start_timer!(|| "Reading tau g1");
        let mut tau_powers_g1 = batch_into_affine(
            &domain.fft(&batch_into_projective(
                &(0..n_taus)
                    .map(|_| E::read_radix_g1(&mut reader))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
        );
        end_timer!(tau_g1_timer);

        let tau_g2_timer = start_timer!(|| "Reading tau g2");
        let tau_powers_g2 = batch_into_affine(
            &domain.fft(&batch_into_projective(
                &(0..n_taus)
                    .map(|_| E::read_radix_g2(&mut reader))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
        );
        end_timer!(tau_g2_timer);

        let alpha_tau_g1_timer = start_timer!(|| "Reading alpha tau g1");
        let alpha_tau_powers_g1 = batch_into_affine(
            &domain.fft(&batch_into_projective(
                &(0..n_taus)
                    .map(|_| E::read_radix_g1(&mut reader))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
        );
        end_timer!(alpha_tau_g1_timer);

        let beta_tau_g1_timer = start_timer!(|| "Reading beta tau g1");
        let beta_tau_powers_g1 = batch_into_affine(
            &domain.fft(&batch_into_projective(
                &(0..n_taus)
                    .map(|_| E::read_radix_g1(&mut reader))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
        );
        end_timer!(beta_tau_g1_timer);

        let tau_g1_aux_timer = start_timer!(|| "Reading tau g1 auxilary");
        tau_powers_g1.append(
            &mut (0..n_taus - 1)
                .map(|i| -> Result<_, Error> {
                    Ok(E::read_radix_g1(&mut reader)? + tau_powers_g1[i])
                })
                .collect::<Result<Vec<_>, _>>()?,
        );
        end_timer!(tau_g1_aux_timer);

        end_timer!(timer);

        Ok(Accumulator {
            tau_powers_g1,
            tau_powers_g2,
            alpha_tau_powers_g1,
            beta_tau_powers_g1,
            beta_g2,
        })
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
#[must_use = "Prepared accumulator must be used in constraint generation"]
pub struct PreparedAccumulator<E: PairingEngine> {
    pub alpha: E::G1Affine,
    pub beta: E::G1Affine,
    pub tau_lagrange_g1: Vec<E::G1Projective>,
    pub tau_lagrange_g2: Vec<E::G2Projective>,
    pub alpha_lagrange_g1: Vec<E::G1Projective>,
    pub beta_lagrange_g1: Vec<E::G1Projective>,
    pub beta_g2: E::G2Affine,
    pub h_query: Vec<E::G1Affine>,
}

impl<E: PairingEngine> PreparedAccumulator<E> {
    pub fn check_pow_len(&self) -> (bool, usize) {
        let len = self.tau_lagrange_g1.len();

        (
            len == self.tau_lagrange_g2.len()
                && len == self.alpha_lagrange_g1.len()
                && len == self.beta_lagrange_g1.len()
                && (len - 1) == self.h_query.len(),
            len,
        )
    }
}

impl<E: PairingEngine + PairingReader> PreparedAccumulator<E> {
    pub fn from_radix_file<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        let timer = start_timer!(|| "Reading from phase1radix file");
        let path: &Path = path.as_ref();

        let n_taus = 1usize
            << path
                .iter()
                .last()
                .and_then(|p| {
                    p.to_str()
                        .and_then(|e| e.split("phase1radix2m").last())
                        .and_then(|e| e.parse::<usize>().ok())
                })
                .ok_or(crate::error::Error::Read(String::from(
                    "Unable to read evaluation domain from radix file name",
                )))?;

        let mut reader = BufReader::new(File::open(path)?);

        let aux_timer = start_timer!(|| "Reading auxilary");
        let alpha = E::read_radix_g1(&mut reader)?;
        let beta = E::read_radix_g1(&mut reader)?;
        end_timer!(aux_timer);

        let beta_g2_timer = start_timer!(|| "Reading beta g2");
        let beta_g2 = E::read_radix_g2(&mut reader)?;
        end_timer!(beta_g2_timer);

        let tau_g1_timer = start_timer!(|| "Reading tau g1");
        let tau_lagrange_g1 = (0..n_taus)
            .map(|_| E::read_radix_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(tau_g1_timer);

        let tau_g2_timer = start_timer!(|| "Reading tau g2");
        let tau_lagrange_g2 = (0..n_taus)
            .map(|_| E::read_radix_g2(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(tau_g2_timer);

        let alpha_tau_g1_timer = start_timer!(|| "Reading alpha tau g1");
        let alpha_lagrange_g1 = (0..n_taus)
            .map(|_| E::read_radix_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(alpha_tau_g1_timer);

        let beta_tau_g1_timer = start_timer!(|| "Reading beta tau g1");
        let beta_lagrange_g1 = (0..n_taus)
            .map(|_| E::read_radix_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(beta_tau_g1_timer);

        let h_query_timer = start_timer!(|| "Reading h query");
        let h_query = (0..n_taus - 1)
            .map(|_| E::read_radix_g1(&mut reader))
            .collect::<Result<Vec<_>, _>>()?;
        end_timer!(h_query_timer);

        end_timer!(timer);

        Ok(PreparedAccumulator {
            alpha,
            beta,
            tau_lagrange_g1: batch_into_projective(&tau_lagrange_g1),
            tau_lagrange_g2: batch_into_projective(&tau_lagrange_g2),
            alpha_lagrange_g1: batch_into_projective(&alpha_lagrange_g1),
            beta_lagrange_g1: batch_into_projective(&beta_lagrange_g1),
            beta_g2,
            h_query,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use ark_bls12_381::Bls12_381;
    use ark_bn254::{Bn254, Fr};
    use ark_ff::{PrimeField, Zero};
    use ark_groth16::Groth16;
    use ark_r1cs_std::{
        fields::fp::FpVar,
        prelude::{AllocVar, EqGadget},
    };
    use ark_relations::r1cs::ConstraintSynthesizer;
    use ark_snark::SNARK;
    use rand::rngs::OsRng;

    use crate::{accumulator::PreparedAccumulator, transcript::Transcript};

    use super::Accumulator;

    struct DummyCircuit<F: PrimeField> {
        pub a: F,
        pub b: F,
        pub c: F,
    }

    impl<F: PrimeField> ConstraintSynthesizer<F> for DummyCircuit<F> {
        fn generate_constraints(
            self,
            cs: ark_relations::r1cs::ConstraintSystemRef<F>,
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
                a: ark_bn254::Fr::zero(),
                b: ark_bn254::Fr::zero(),
                c: ark_bn254::Fr::zero(),
            },
        )?;

        let proof = Groth16::prove(
            &transcript.key.key,
            DummyCircuit {
                a: ark_bn254::Fr::from(1),
                b: ark_bn254::Fr::from(1),
                c: ark_bn254::Fr::from(2),
            },
            rng,
        )?;

        let valid = Groth16::verify(&transcript.key.key.vk, &[Fr::from(2)], &proof)?;
        assert!(valid, "Proof must be valid");

        let valid = Groth16::verify(&transcript.key.key.vk, &[Fr::from(4)], &proof)?;
        assert!(!valid, "Proof must be not valid");

        Ok(())
    }

    #[test]
    fn from_radix_file_works() -> Result<(), Box<dyn Error>> {
        let rng = &mut OsRng;
        let radix_path = "phase1radix2m2";

        let accum = PreparedAccumulator::<Bls12_381>::from_radix_file(radix_path)?;
        let transcript = Transcript::new_from_prepared_accumulator(
            &accum,
            DummyCircuit {
                a: ark_bls12_381::Fr::zero(),
                b: ark_bls12_381::Fr::zero(),
                c: ark_bls12_381::Fr::zero(),
            },
        )?;

        let proof = Groth16::prove(
            &transcript.key.key,
            DummyCircuit {
                a: ark_bls12_381::Fr::from(1),
                b: ark_bls12_381::Fr::from(1),
                c: ark_bls12_381::Fr::from(2),
            },
            rng,
        )?;

        let valid = Groth16::verify(
            &transcript.key.key.vk,
            &[ark_bls12_381::Fr::from(2)],
            &proof,
        )?;
        assert!(valid, "Proof must be valid");

        let valid = Groth16::verify(
            &transcript.key.key.vk,
            &[ark_bls12_381::Fr::from(4)],
            &proof,
        )?;
        assert!(!valid, "Proof must be not valid");

        Ok(())
    }

    #[test]
    fn from_radix_bn254_file_works() -> Result<(), Box<dyn Error>> {
        let rng = &mut OsRng;
        let radix_path = "bn254phase1radix2m2";

        let accum = PreparedAccumulator::<Bn254>::from_radix_file(radix_path)?;
        let transcript = Transcript::new_from_prepared_accumulator(
            &accum,
            DummyCircuit {
                a: ark_bn254::Fr::zero(),
                b: ark_bn254::Fr::zero(),
                c: ark_bn254::Fr::zero(),
            },
        )?;

        let proof = Groth16::prove(
            &transcript.key.key,
            DummyCircuit {
                a: ark_bn254::Fr::from(1),
                b: ark_bn254::Fr::from(1),
                c: ark_bn254::Fr::from(2),
            },
            rng,
        )?;

        let valid = Groth16::verify(&transcript.key.key.vk, &[ark_bn254::Fr::from(2)], &proof)?;
        assert!(valid, "Proof must be valid");

        let valid = Groth16::verify(&transcript.key.key.vk, &[ark_bn254::Fr::from(4)], &proof)?;
        assert!(!valid, "Proof must be not valid");

        Ok(())
    }
}
