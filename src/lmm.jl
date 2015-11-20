##################################################################
# Fast linear mixed models
##################################################################
#
# We implement linear mixed models for data which has covariance of
# the form tau2*K + sigma2*I, where sigma2 and tau2 are positive
# scalars, K is a symmetric positive definite "kinship" matrix and I
# is the identity matrix.
#
# Models of this type are used for GWAS and QTL mapping in structured
# populations.
#
# ################################################################

"""
rotateData: Rotates data with respect to the kinship matrix

    y = phenotype matrix
    X = predictor matrix
    K = kinship matrix, expected to be symmetric and positive definite
"""

function rotateData(y,X,K)
    # check dimensions
    n = size(y,1)
    if( size(X,1) != n )
        error("Dimension mismatch.")
    elseif( size(K,1) != n )
        error("Dimension mismatch.")
    end

    # check symmetry and positive definiteness of K
    if( !issym(K) )
        error("K is not symmetric.")
    elseif( !isposdef(K) )
        error("K is not positive definite.")
    end

    # spectral decomposition of a symmetric matrix
    EF = eigfact(K)

    # return rotated phenotype, covariates, and eigenvalues
    return EF[:vectors]'y, EF[:vectors]'X, EF[:values]

end



"""
wls: Weighted least squares estimation

    y = outcome
    X = predictors
    w = weights (should be positive)

The variance estimate is maximum likelihood
"""

function wls(y,X,w)

    # number ofindividuals
    n = size(y,1)

    # square root of the weights
    sqrtw = sqrt(w)
    # scale by weights
    yy = diagm(sqrtw)*y
    XX = diagm(sqrtw)*X

    # QR decomposition of the transformed data
    (q,r) = qr(XX)
    yy = At_mul_B(q,yy)
    b = r\yy

    # estimate y and calculate rss
    yyhat = XX*b
    rss = norm(yy-yyhat)^2

    # return coefficient and variance estimate
    return b, rss/n

end

# ################################################################
# function to calculate log likelihood of data given fixed effects
# ################################################################



# ################################################################
# function to estimate heritability given fixed effects
# ################################################################



