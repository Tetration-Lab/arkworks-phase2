use std::io::Read;

use ark_ec::{
    bls12::{Bls12, Bls12Parameters},
    bn::{Bn, BnParameters},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::{PrimeField, QuadExtField, Zero};
use ark_serialize::CanonicalSerialize;

use crate::error::Error;

pub trait PairingReader: PairingEngine {
    fn read_g1<R: Read>(reader: &mut R) -> Result<Self::G1Affine, Error>;
    fn read_g2<R: Read>(reader: &mut R) -> Result<Self::G2Affine, Error>;
}

impl<T: BnParameters> PairingReader for Bn<T> {
    fn read_g1<R: Read>(reader: &mut R) -> Result<Self::G1Affine, Error> {
        let zero = Self::Fq::zero();
        let len = zero.uncompressed_size();
        let mut b = vec![0; len * 2];
        reader.read_exact(&mut b)?;
        let x = Self::Fq::from_le_bytes_mod_order(&b[..len]);
        let y = Self::Fq::from_le_bytes_mod_order(&b[len..]);
        let g1 = Self::G1Projective::new(x, y, zero).into_affine();

        assert!(g1.is_on_curve(), "G1 must be on curve");

        Ok(g1)
    }

    fn read_g2<R: Read>(reader: &mut R) -> Result<Self::G2Affine, Error> {
        let zero = Self::Fq::zero();
        let fp2_zero = <Self::G2Affine as AffineCurve>::BaseField::zero();
        let len = zero.uncompressed_size();
        let mut b = vec![0; len * 4];
        reader.read_exact(&mut b)?;
        let x0 = Self::Fq::from_le_bytes_mod_order(&b[..len]);
        let x1 = Self::Fq::from_le_bytes_mod_order(&b[len..len * 2]);
        let x = QuadExtField::new(x0, x1);
        let y0 = Self::Fq::from_le_bytes_mod_order(&b[len * 2..len * 3]);
        let y1 = Self::Fq::from_le_bytes_mod_order(&b[len * 3..]);
        let y = QuadExtField::new(y0, y1);
        let g2 = Self::G2Projective::new(x, y, fp2_zero).into_affine();

        assert!(g2.is_on_curve(), "G2 must be on curve");

        Ok(g2)
    }
}

impl<T: Bls12Parameters> PairingReader for Bls12<T> {
    fn read_g1<R: Read>(reader: &mut R) -> Result<Self::G1Affine, Error> {
        let zero = Self::Fq::zero();
        let len = zero.uncompressed_size();
        let mut b = vec![0; len * 2];
        reader.read_exact(&mut b)?;
        let x = Self::Fq::from_le_bytes_mod_order(&b[..len]);
        let y = Self::Fq::from_le_bytes_mod_order(&b[len..]);
        let g1 = Self::G1Projective::new(x, y, zero).into_affine();

        assert!(g1.is_on_curve(), "G1 must be on curve");

        Ok(g1)
    }

    fn read_g2<R: Read>(reader: &mut R) -> Result<Self::G2Affine, Error> {
        let zero = Self::Fq::zero();
        let fp2_zero = <Self::G2Affine as AffineCurve>::BaseField::zero();
        let len = zero.uncompressed_size();
        let mut b = vec![0; len * 4];
        reader.read_exact(&mut b)?;
        let x0 = Self::Fq::from_le_bytes_mod_order(&b[..len]);
        let x1 = Self::Fq::from_le_bytes_mod_order(&b[len..len * 2]);
        let x = QuadExtField::new(x0, x1);
        let y0 = Self::Fq::from_le_bytes_mod_order(&b[len * 2..len * 3]);
        let y1 = Self::Fq::from_le_bytes_mod_order(&b[len * 3..]);
        let y = QuadExtField::new(y0, y1);
        let g2 = Self::G2Projective::new(x, y, fp2_zero).into_affine();

        assert!(g2.is_on_curve(), "G2 must be on curve");

        Ok(g2)
    }
}
