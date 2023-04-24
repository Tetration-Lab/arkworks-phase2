pub mod error;
pub mod hasher;
pub mod key;
pub mod keypair;
pub mod pot;
pub mod ratio;
pub mod reader;
pub mod transcript;
pub mod utils;

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

    use crate::{pot::Accumulator, transcript::Transcript};

    const NUM_CONSTRAINTS: usize = 50;

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

            for _ in 0..(NUM_CONSTRAINTS - 5) {
                c.enforce_equal(&d)?;
            }

            Ok(())
        }
    }

    #[test]
    fn groth16_domain() -> Result<(), Box<dyn Error>> {
        let rng = &mut OsRng;
        let mut accum = Accumulator::<Bn254>::empty_from_max_constraints(NUM_CONSTRAINTS)?;
        accum.contribute(rng);

        let mut transcript = Transcript::new_from_accumulator(
            &accum,
            DummyCircuit {
                a: Fr::zero(),
                b: Fr::zero(),
                c: Fr::zero(),
            },
        )?;

        transcript.contribute_seed(b"ofekwoo")?;
        transcript.verify()?;

        transcript.contribute_seed(b"aewrog")?;
        transcript.verify()?;

        let pk = transcript.key.key;
        let proof = Groth16::prove(
            &pk,
            DummyCircuit {
                a: Fr::from(1),
                b: Fr::from(2),
                c: Fr::from(3),
            },
            rng,
        )?;
        let valid = Groth16::verify(&pk.vk, &[Fr::from(3)], &proof)?;
        assert!(valid, "Proof must be valid");

        Ok(())
    }
}
