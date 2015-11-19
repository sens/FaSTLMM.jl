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

# y = phenotype matrix
# X = predictor matrix
# K = kinship matrix, expected to be symmetric and positive definite

# ################################################################
# function to rotate the data with respect to the kinship matrix
# ################################################################

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


# ################################################################
# function to perform weighted least squares estimation
# ################################################################

# y = outcome
# X = predictors
# w = weights (should be positive)

function wls(y,X,w)

    sqrtw = sqrt(w)
    y = diagm(sqrtw)*y
    X = diagm(sqrtw)*X

    (q,r) = qr(X)
    y = At_mul_B(q,y)
    b = r\y

    yhat = X*b
    rss = norm(y-yhat)^2
    
    return b, rss
    
end

# ################################################################
# function to estimate error variance
# ################################################################


# ################################################################
# function to estimate heritability
# ################################################################
