# FaST-LMM: *Fa*ctored *S*pectrally *T*ransformed *L*inear *M*ixed *M*odels

Genetic analysis in structured populations used mixed linear models
where the variance matrix of the error term is a linear combination of
an identity matrix and a positive definite matrix.

The linear model is $$y = X\beta + e,$$ where $V(e) = \sigma_G^2 K +
\sigma_E^2 I$.

This implementation is my attempt to learn Julia and numerical linear
algebra.
