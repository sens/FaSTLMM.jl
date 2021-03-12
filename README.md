# FaST-LMM: *Fa*ctored *S*pectrally *T*ransformed *L*inear *M*ixed *M*odels

[![Build Status](https://travis-ci.org/sens/FaSTLMM.jl.svg?branch=master)](https://travis-ci.org/sens/FaSTLMM.jl)

Genetic analysis in structured populations used mixed linear models
where the variance matrix of the error term is a linear combination of
an identity matrix and a positive definite matrix.

The linear model is of the familiar form: 𝑦 = 𝑋 β + ϵ.

- 𝑦: phenotype
-𝑋: covariates
- β: fixed effects
- ϵ: error term

Further, V(ϵ) = τ²𝐾+ σ²𝐼, where τ² is
the genetic variance, σ² is the environmental variance, 𝐾
is the kinship matrix, and 𝐼 is the identity matrix.

The key idea in speeding up computations here is that by rotating the
phenotypes by the eigenvectors of $K$ we can transform estimation to a
weighted least squares problem.

This implementation is my attempt to learn Julia and numerical linear
algebra.  The code is being tested.

Guide to the directories:

- `src`: Julia source code
- `data`: Example data for development and testing
- `test`: Code for testing
- `docs`: Notes on comparisons with other implementations
