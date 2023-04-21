use std::borrow::Borrow;

use ark_ec::PairingEngine;
use ark_groth16::ProvingKey;

use crate::{error::Error, utils::serialize};

#[derive(Debug, Clone, PartialEq)]
pub struct FullKey<E: PairingEngine> {
    pub key: ProvingKey<E>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PartialKey<E: PairingEngine> {
    pub delta_g1: E::G1Affine,
    pub delta_g2: E::G2Affine,
    pub h_query: Vec<E::G1Affine>,
    pub l_query: Vec<E::G1Affine>,
}

impl<E: PairingEngine> Borrow<PartialKey<E>> for FullKey<E> {
    fn borrow(&self) -> &PartialKey<E> {
        todo!()
    }
}

impl<E: PairingEngine> PartialEq<PartialKey<E>> for FullKey<E> {
    fn eq(&self, other: &PartialKey<E>) -> bool {
        self.key.delta_g1 == other.delta_g1
            && self.key.vk.delta_g2 == other.delta_g2
            && self.key.h_query == other.h_query
            && self.key.l_query == other.l_query
    }
}

impl<E: PairingEngine> From<ProvingKey<E>> for FullKey<E> {
    fn from(val: ProvingKey<E>) -> Self {
        FullKey { key: val }
    }
}

impl<E: PairingEngine> From<FullKey<E>> for ProvingKey<E> {
    fn from(val: FullKey<E>) -> Self {
        val.key
    }
}

impl<E: PairingEngine> From<ProvingKey<E>> for PartialKey<E> {
    fn from(val: ProvingKey<E>) -> Self {
        PartialKey {
            delta_g1: val.delta_g1,
            delta_g2: val.vk.delta_g2,
            h_query: val.h_query,
            l_query: val.l_query,
        }
    }
}

impl<E: PairingEngine> From<&'_ FullKey<E>> for PartialKey<E> {
    fn from(val: &'_ FullKey<E>) -> Self {
        PartialKey {
            delta_g1: val.key.delta_g1,
            delta_g2: val.key.vk.delta_g2,
            h_query: val.key.h_query.clone(),
            l_query: val.key.l_query.clone(),
        }
    }
}

impl<E: PairingEngine> From<&'_ ProvingKey<E>> for PartialKey<E> {
    fn from(val: &'_ ProvingKey<E>) -> Self {
        PartialKey {
            delta_g1: val.delta_g1,
            delta_g2: val.vk.delta_g2,
            h_query: val.h_query.clone(),
            l_query: val.l_query.clone(),
        }
    }
}

impl<E: PairingEngine> PartialKey<E> {
    pub fn challenge(&self) -> Result<Vec<u8>, Error> {
        Ok(serialize(&self.delta_g1)?
            .into_iter()
            .chain(serialize(&self.delta_g2)?)
            .collect())
    }
}

impl<E: PairingEngine> FullKey<E> {
    pub fn challenge(&self) -> Result<Vec<u8>, Error> {
        Ok(serialize(&self.key.delta_g1)?
            .into_iter()
            .chain(serialize(&self.key.vk.delta_g2)?)
            .collect())
    }
}
