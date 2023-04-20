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
    pub(crate) fn generate<R: Rng>(g1_powers: usize, g2_powers: usize, rng: &mut R) -> Self {
        use std::iter;

        use ark_ec::AffineCurve;
        use ark_ff::{One, UniformRand};
        use ark_std::cfg_iter_mut;

        let tau = E::Fr::rand(rng);
        let alpha = E::Fr::rand(rng);
        let beta = E::Fr::rand(rng);

        let mut instance = Self {
            tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g1_powers],
            tau_powers_g2: vec![E::G2Affine::prime_subgroup_generator(); g2_powers],
            alpha_tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g2_powers],
            beta_tau_powers_g1: vec![E::G1Affine::prime_subgroup_generator(); g2_powers],
            beta_g2: E::G2Affine::prime_subgroup_generator(),
        };

        let mut tau_powers = iter::successors(Some(E::Fr::one()), |x| Some(*x * tau))
            .take(g1_powers)
            .collect::<Vec<_>>();
        let remaining_tau_powers = tau_powers.split_off(g2_powers);
        cfg_iter_mut!(instance.tau_powers_g1)
            .zip(cfg_iter_mut!(instance.tau_powers_g2))
            .zip(cfg_iter_mut!(instance.alpha_tau_powers_g1))
            .zip(cfg_iter_mut!(instance.beta_tau_powers_g1))
            .zip(tau_powers)
            .for_each(
                |((((tau_g1, tau_g2), alpha_tau_g1), beta_tau_g1), tau_power)| {
                    *tau_g1 = tau_g1.mul(tau_power).into();
                    *tau_g2 = tau_g2.mul(tau_power).into();
                    *alpha_tau_g1 = alpha_tau_g1.mul(tau_power * alpha).into();
                    *beta_tau_g1 = beta_tau_g1.mul(tau_power * beta).into();
                },
            );
        cfg_iter_mut!(instance.tau_powers_g1)
            .skip(g2_powers)
            .zip(remaining_tau_powers)
            .for_each(|(tau_g1, tau_power)| *tau_g1 = tau_g1.mul(tau_power).into());
        instance.beta_g2 = instance.beta_g2.mul(beta).into();

        instance
    }
}
