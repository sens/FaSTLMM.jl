# FaST-LMM: *Fa*ctored *S*pectrally *T*ransformed *L*inear *M*ixed *M*odels

Genetic analysis in structured populations used mixed linear models
where the variance matrix of the error term is a linear combination of
an identity matrix and a positive definite matrix.

The linear model is of the familiar form: $$y = X \beta + e$$

- $y$: phenotype
- $X$: covariates
- $\beta$: fixed effects
- $e$: error term

 Further $V(e) = \sigma_G^2 K + \sigma_E^2 I$, where $\sigma_G^2$ is
 the genetic variance, $\sigma_E^2$ is the environmental variance, $K$
 is the kinship matrix, and $I$ is the identity matrix.

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
