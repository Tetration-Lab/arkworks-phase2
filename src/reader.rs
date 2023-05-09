use std::{io::Read, ops::Div};

use ark_ec::{
    bls12::{Bls12, Bls12Parameters},
    bn::{Bn, BnParameters},
    PairingEngine,
};
use ark_ff::{BigInteger, FpParameters, FromBytes, PrimeField, QuadExtField};

use crate::error::Error;

pub trait PairingReader: PairingEngine {
    fn read_ptau_fq<R: Read>(reader: &mut R) -> Result<Self::Fq, Error>;
    fn read_ptau_g1<R: Read>(reader: &mut R) -> Result<Self::G1Affine, Error>;
    fn read_ptau_g2<R: Read>(reader: &mut R) -> Result<Self::G2Affine, Error>;
}

impl<T: BnParameters> PairingReader for Bn<T> {
    fn read_ptau_fq<R: Read>(reader: &mut R) -> Result<Self::Fq, Error> {
        let mut b = vec![0; <Self::Fq as PrimeField>::BigInt::NUM_LIMBS * 8];
        reader.read_exact(&mut b)?;
        Ok(
            Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[..])?)
                .ok_or(Error::Read(String::from("Unable to read Fq")))?
                .div(Self::Fq::from_repr(<Self::Fq as PrimeField>::Params::R).unwrap()),
        )
    }

    fn read_ptau_g1<R: Read>(reader: &mut R) -> Result<Self::G1Affine, Error> {
        let r = Self::Fq::from_repr(<Self::Fq as PrimeField>::Params::R).unwrap();
        let len = <Self::Fq as PrimeField>::BigInt::NUM_LIMBS * 8;
        let mut b = vec![0; len * 2];
        reader.read_exact(&mut b)?;

        let x = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[..len])?)
            .ok_or(Error::Read(String::from("Unable to read Fq x")))?
            .div(r);
        let y = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[len..])?)
            .ok_or(Error::Read(String::from("Unable to read Fq y")))?
            .div(r);

        let g1 = Self::G1Affine::new(x, y, false);

        assert!(g1.is_on_curve(), "G1 must be on curve");

        Ok(g1)
    }

    fn read_ptau_g2<R: Read>(reader: &mut R) -> Result<Self::G2Affine, Error> {
        let r = Self::Fq::from_repr(<Self::Fq as PrimeField>::Params::R).unwrap();
        let len = <Self::Fq as PrimeField>::BigInt::NUM_LIMBS * 8;
        let mut b = vec![0; len * 4];
        reader.read_exact(&mut b)?;

        let x0 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[..len])?)
            .ok_or(Error::Read(String::from("Unable to read Fq x0")))?
            .div(r);
        let x1 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[len..(len * 2)])?)
            .ok_or(Error::Read(String::from("Unable to read Fq x1")))?
            .div(r);
        let x = QuadExtField::new(x0, x1);
        let y0 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(
            &b[(len * 2)..(len * 3)],
        )?)
        .ok_or(Error::Read(String::from("Unable to read Fq y0")))?
        .div(r);
        let y1 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[(len * 3)..])?)
            .ok_or(Error::Read(String::from("Unable to read Fq y1")))?
            .div(r);
        let y = QuadExtField::new(y0, y1);
        let g2 = Self::G2Affine::new(x, y, false);

        assert!(g2.is_on_curve(), "G2 must be on curve");

        Ok(g2)
    }
}

impl<T: Bls12Parameters> PairingReader for Bls12<T> {
    fn read_ptau_fq<R: Read>(reader: &mut R) -> Result<Self::Fq, Error> {
        let mut b = vec![0; <Self::Fq as PrimeField>::BigInt::NUM_LIMBS * 8];
        reader.read_exact(&mut b)?;
        Ok(
            Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[..])?)
                .ok_or(Error::Read(String::from("Unable to read Fq")))?
                .div(Self::Fq::from_repr(<Self::Fq as PrimeField>::Params::R).unwrap()),
        )
    }

    fn read_ptau_g1<R: Read>(reader: &mut R) -> Result<Self::G1Affine, Error> {
        let r = Self::Fq::from_repr(<Self::Fq as PrimeField>::Params::R).unwrap();
        let len = <Self::Fq as PrimeField>::BigInt::NUM_LIMBS * 8;
        let mut b = vec![0; len * 2];
        reader.read_exact(&mut b)?;

        let x = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[..len])?)
            .ok_or(Error::Read(String::from("Unable to read Fq x")))?
            .div(r);
        let y = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[len..])?)
            .ok_or(Error::Read(String::from("Unable to read Fq y")))?
            .div(r);

        let g1 = Self::G1Affine::new(x, y, false);

        assert!(g1.is_on_curve(), "G1 must be on curve");

        Ok(g1)
    }

    fn read_ptau_g2<R: Read>(reader: &mut R) -> Result<Self::G2Affine, Error> {
        let r = Self::Fq::from_repr(<Self::Fq as PrimeField>::Params::R).unwrap();
        let len = <Self::Fq as PrimeField>::BigInt::NUM_LIMBS * 8;
        let mut b = vec![0; len * 4];
        reader.read_exact(&mut b)?;

        let x0 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[..len])?)
            .ok_or(Error::Read(String::from("Unable to read Fq x0")))?
            .div(r);
        let x1 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[len..(len * 2)])?)
            .ok_or(Error::Read(String::from("Unable to read Fq x1")))?
            .div(r);
        let x = QuadExtField::new(x0, x1);
        let y0 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(
            &b[(len * 2)..(len * 3)],
        )?)
        .ok_or(Error::Read(String::from("Unable to read Fq y0")))?
        .div(r);
        let y1 = Self::Fq::from_repr(<Self::Fq as PrimeField>::BigInt::read(&b[(len * 3)..])?)
            .ok_or(Error::Read(String::from("Unable to read Fq y1")))?
            .div(r);
        let y = QuadExtField::new(y0, y1);
        let g2 = Self::G2Affine::new(x, y, false);

        assert!(g2.is_on_curve(), "G2 must be on curve");

        Ok(g2)
    }
}
