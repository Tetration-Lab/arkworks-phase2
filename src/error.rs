use ark_relations::r1cs::SynthesisError;

use thiserror::Error;

#[derive(Error, Clone, Debug, Eq, PartialEq)]
pub enum Error {
    #[error("Too Many Constraints")]
    TooManyConstraints,

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

    #[error("Invariant Violation")]
    InvariantViolated(&'static str),

    #[error("{0}")]
    Custom(String),
}
