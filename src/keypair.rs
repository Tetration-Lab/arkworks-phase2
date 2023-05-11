use ark_ec::PairingEngine;

use crate::{error::Error, ratio::RatioProof, utils::serialize};

pub type PrivateKey<E> = <E as PairingEngine>::Fr;

#[derive(Debug, Clone, PartialEq)]
pub struct PublicKey<E: PairingEngine> {
    pub delta_g2: E::G2Affine,
    pub proof: RatioProof<E>,
}

impl<E: PairingEngine> PublicKey<E> {
    pub fn challenge(&self) -> Result<Vec<u8>, Error> {
        serialize(&self.delta_g2)
    }
}
