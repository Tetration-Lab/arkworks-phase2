#![allow(clippy::double_must_use)]

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::{cfg_into_iter, cfg_iter, cfg_iter_mut};
use rand::thread_rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::error::Error;

#[inline]
#[must_use]
pub fn batch_into_projective<A: AffineCurve>(points: &[A]) -> Vec<A::Projective> {
    cfg_iter!(points).map(A::into_projective).collect()
}

#[inline]
pub fn batch_mul_fixed_scalar<A: AffineCurve>(points: &mut [A], scalar: A::ScalarField) {
    cfg_iter_mut!(points).for_each(|point| *point = point.mul(scalar).into_affine())
}

#[inline]
#[must_use]
pub fn same_ratio<E: PairingEngine>(
    lhs: (E::G1Affine, E::G2Affine),
    rhs: (E::G1Affine, E::G2Affine),
) -> bool {
    E::pairing(lhs.0, lhs.1) == E::pairing(rhs.0, rhs.1)
}

#[inline]
#[must_use]
pub fn same_ratio_swap<E: PairingEngine>(
    lhs: (E::G1Affine, E::G1Affine),
    rhs: (E::G2Affine, E::G2Affine),
) -> bool {
    E::pairing(lhs.0, rhs.1) == E::pairing(lhs.1, rhs.0)
}

#[must_use]
pub fn seeded_rng(bytes: &[u8]) -> ChaCha20Rng {
    let seed = blake3::hash(bytes);
    ChaCha20Rng::from_seed(*seed.as_bytes())
}

#[must_use]
pub fn serialize<T: CanonicalSerialize>(value: &T) -> Result<Vec<u8>, Error> {
    let mut output = vec![];
    value
        .serialize(&mut output)
        .map_err(|e| Error::Custom(e.to_string()))?;
    Ok(output)
}

#[must_use]
pub fn serialize_uncompressed<T: CanonicalSerialize>(value: &T) -> Result<Vec<u8>, Error> {
    let mut output = vec![];
    value
        .serialize_uncompressed(&mut output)
        .map_err(|e| Error::Custom(e.to_string()))?;
    Ok(output)
}

#[must_use]
#[inline]
pub fn merge_ratio_affine_vec<A: AffineCurve>(lhs: &[A], rhs: &[A]) -> (A, A) {
    assert_eq!(lhs.len(), rhs.len(), "lhs.len() != rhs.len()");
    let (mut l, mut r) = (A::zero(), A::zero());
    cfg_into_iter!(0..lhs.len())
        .map(|_| A::ScalarField::rand(&mut thread_rng()))
        .zip(lhs)
        .zip(rhs)
        .for_each(|((s, l1), r1)| {
            l = l + l1.mul(s).into_affine();
            r = r + r1.mul(s).into_affine();
        });
    (l, r)
}
