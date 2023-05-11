use std::sync::PoisonError;

use ark_relations::r1cs::SynthesisError;
use thiserror::Error;

#[derive(Error, Clone, Debug, Eq, PartialEq)]
pub enum Error {
    #[error("Constraint System: {0}")]
    ConstraintSystem(#[from] SynthesisError),

    #[error("Try From Int: {0}")]
    TryFromInt(#[from] std::num::TryFromIntError),

    #[error("IO Error: {0}")]
    IO(String),

    #[error("Serialization Error: {0}")]
    Serialization(String),

    #[error("Lock Already Held")]
    Locked,

    #[error("Missing Constraint System Matrices")]
    MissingCSMatrices,

    #[error("Constraint System Hashes Mismatch")]
    ConstraintSystemHashesDiffer,

    #[error("Invalid Ratio Proof")]
    InvalidRatioProof,

    #[error("Inconsistent Delta Change Error")]
    InconsistentDeltaChange,

    #[error("Inconsistent L Change")]
    InconsistentLChange,

    #[error("Inconsistent H Change")]
    InconsistentHChange,

    #[error("Invariant Violation {0}")]
    InvariantViolated(&'static str),

    #[error("Invalid key element {0}")]
    InvalidKey(&'static str),

    #[error("Invalid initial partial key")]
    InvalidPartialKey,

    #[error("Ptau Field Check Not Matching")]
    PtauFieldNotMatching,

    #[error("Too Enough POT Degree, Needed Atleast {0}")]
    NotEnoughPOTDegree(u32),

    #[error("Invalid Power of Tau Degree, Need {0}")]
    InvalidPOTDegree(u32),

    #[error("Invalid Power of Tau Size")]
    InvalidPOTSize,

    #[error("Read Error: {0}")]
    Read(String),

    #[error("{0}")]
    Custom(String),
}

impl<E> From<PoisonError<E>> for Error {
    fn from(_err: PoisonError<E>) -> Self {
        Error::Locked
    }
}

impl From<std::io::Error> for Error {
    fn from(value: std::io::Error) -> Self {
        Self::IO(value.to_string())
    }
}

impl From<ark_serialize::SerializationError> for Error {
    fn from(value: ark_serialize::SerializationError) -> Self {
        Self::Serialization(value.to_string())
    }
}
