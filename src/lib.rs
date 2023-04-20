pub mod error;
pub mod hasher;
pub mod key;
pub mod keypair;
pub mod pot;
pub mod ratio;
pub mod transcript;
pub mod utils;

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

    use crate::{pot::Accumulator, transcript::Transcript, utils::seeded_rng};

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
