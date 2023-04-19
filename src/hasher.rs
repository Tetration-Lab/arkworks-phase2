use std::marker::PhantomData;

use ark_ec::PairingEngine;
use ark_ff::UniformRand;
use rand::Rng;

pub struct HashToCurve<E: PairingEngine>(PhantomData<E>);

impl<E: PairingEngine> HashToCurve<E> {
    pub fn hash_g1<R: Rng>(rng: &mut R) -> E::G1Affine {
        <E::G1Projective as UniformRand>::rand(rng).into()
    }
    pub fn hash_g2<R: Rng>(rng: &mut R) -> E::G2Affine {
        <E::G2Projective as UniformRand>::rand(rng).into()
    }
}
