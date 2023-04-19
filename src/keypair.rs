use ark_ec::PairingEngine;

use crate::{key::PartialKey, ratio::RatioProof};

pub type PrivateKey<E> = <E as PairingEngine>::Fr;

#[derive(Debug, Clone, PartialEq)]
pub struct PublicKey<E: PairingEngine> {
    pub partial_key: PartialKey<E>,
    pub proof: RatioProof<E>,
}
