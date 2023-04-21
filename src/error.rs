use ark_relations::r1cs::SynthesisError;

use thiserror::Error;

#[derive(Error, Clone, Debug, Eq, PartialEq)]
pub enum Error {
    #[error("Constraint System: {0}")]
    ConstraintSystem(#[from] SynthesisError),

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

    #[error("Too Enough Tau Powers G2, Needed Atleast {0}")]
    NotEnoughTauPowers(usize),

    #[error("Invalid Perpetual Power of Tau Powers")]
    InvalidPPOTPowers,

    #[error("{0}")]
    Custom(String),
}
