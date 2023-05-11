#![allow(clippy::double_must_use)]

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::{cfg_iter, cfg_iter_mut};
use rand::rngs::OsRng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::error::Error;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[inline]
#[must_use]
pub fn batch_into_projective<A: AffineCurve>(points: &[A]) -> Vec<A::Projective> {
    cfg_iter!(points).map(A::into_projective).collect()
}

#[inline]
#[must_use]
pub fn batch_into_affine<P: ProjectiveCurve>(points: &[P]) -> Vec<P::Affine> {
    let mut points = points.to_vec();
    P::batch_normalization(&mut points);
    cfg_iter!(points).map(P::into_affine).collect()
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
    #[cfg(not(feature = "parallel"))]
    let result = {
        let (mut l, mut r) = (A::zero(), A::zero());
        (0..lhs.len())
            .map(|_| A::ScalarField::rand(&mut OsRng))
            .zip(lhs)
            .zip(rhs)
            .for_each(|((s, l1), r1)| {
                l = l + l1.mul(s).into_affine();
                r = r + r1.mul(s).into_affine();
            });
        (l, r)
    };

    #[cfg(feature = "parallel")]
    let result = {
        lhs.into_par_iter()
            .zip(rhs)
            .map(|(lhs, rhs)| {
                let s = A::ScalarField::rand(&mut OsRng);
                (lhs.mul(s).into_affine(), rhs.mul(s).into_affine())
            })
            .reduce(
                || (A::zero(), A::zero()),
                |(l, r), (l1, r1)| (l + l1, r + r1),
            )
    };

    result
}
