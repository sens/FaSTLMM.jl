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

This implementation is my attempt to learn Julia and numerical linear
algebra.
