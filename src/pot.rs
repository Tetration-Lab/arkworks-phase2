use ark_ec::PairingEngine;

#[cfg(test)]
use rand::Rng;

#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub struct Accumulator<E: PairingEngine> {
    pub tau_powers_g1: Vec<E::G1Affine>,
    pub tau_powers_g2: Vec<E::G2Affine>,
    pub alpha_tau_powers_g1: Vec<E::G1Affine>,
    pub beta_tau_powers_g1: Vec<E::G1Affine>,
    pub beta_g2: E::G2Affine,
}

impl<E: PairingEngine> Accumulator<E> {
    #[cfg(test)]
    pub(crate) fn empty(g1_powers: usize, g2_powers: usize) -> Self {
        use ark_ec::AffineCurve;

        Self {
            tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g1_powers],
            tau_powers_g2: vec![E::G2Affine::prime_subgroup_generator(); g2_powers],
            alpha_tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g2_powers],
            beta_tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g2_powers],
            beta_g2: E::G2Affine::prime_subgroup_generator(),
        }
    }

    #[cfg(test)]
    pub(crate) fn contribute<R: Rng>(&mut self, rng: &mut R) {
        use std::iter;

        use ark_ec::AffineCurve;
        use ark_ff::{One, UniformRand};
        use ark_std::cfg_iter_mut;

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
}

#[cfg(test)]
mod tests {
    use std::{
        collections::BTreeMap,
        error::Error,
        fs::File,
        io::{BufReader, Seek, SeekFrom},
    };

    use ark_bn254::Bn254;
    use ark_ec::{AffineCurve, PairingEngine};
    use ark_ff::{BigInteger, FpParameters, PrimeField, Zero};
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read};

    type Section = BTreeMap<u32, Vec<(u64, u64)>>;
    type Curve = Bn254;

    fn read_bin(buf: &mut BufReader<File>) -> Result<Section, Box<dyn Error>> {
        let mut b = [0u8; 4];
        buf.read_exact(&mut b)?;
        let ty = b.into_iter().map(|x| char::from(x)).collect::<String>();
        println!("Type: {ty}");

        let mut b = [0u8; 4];
        buf.read_exact(&mut b)?;
        let version = u32::from_le_bytes(b);
        println!("Version: {version}");

        let mut b = [0u8; 4];
        buf.read_exact(&mut b)?;
        let n_section = u32::from_le_bytes(b);
        println!("Sections: {n_section}");

        let mut section = Section::new();
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

        println!("Section {:?}", section);

        Ok(section)
    }

    #[test]
    fn from_ptau_file() -> Result<(), Box<dyn Error>> {
        let f = File::open("pot8.ptau")?;
        let mut reader = BufReader::new(f);
        let section = read_bin(&mut reader)?;
        //let mut buffer = Vec::new();

        let header = section
            .get(&1)
            .and_then(|e| e.get(0))
            .expect("Header not found");
        reader.seek(SeekFrom::Start(header.0))?;
        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let n = u32::from_le_bytes(b) as usize;
        println!("Fr length {}", n);

        let mut buffer = vec![0u8; n];
        reader.read_exact(&mut buffer)?;
        println!(
            "Fq bytes matched: {}",
            buffer
                == <<<Curve as PairingEngine>::Fq as PrimeField>::Params as FpParameters>::MODULUS
                    .to_bytes_le()
        );

        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let power = u32::from_le_bytes(b);
        println!("Power {}", power);

        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let ceremony_power = u32::from_le_bytes(b);
        println!("Ceremony power {}", ceremony_power);

        let ceremony = section
            .get(&7)
            .and_then(|e| e.get(0))
            .expect("Ceremony not found");
        reader.seek(SeekFrom::Start(ceremony.0))?;
        let mut b = [0u8; 4];
        reader.read_exact(&mut b)?;
        let n_contribution = u32::from_le_bytes(b);
        println!("Contributions {}", n_contribution);

        let g1_f_size =
            <<<Curve as PairingEngine>::G1Affine as AffineCurve>::BaseField as Zero>::zero()
                .serialized_size();
        let _g2_f_size =
            <<<Curve as PairingEngine>::G2Affine as AffineCurve>::BaseField as Zero>::zero()
                .serialized_size();

        let mut b = vec![0u8; g1_f_size * 2];
        reader.read_exact(&mut b)?;
        let _tau_g1 =
            <<Curve as PairingEngine>::G1Affine as CanonicalDeserialize>::deserialize(&*b)?;

        Ok(())
    }
}
