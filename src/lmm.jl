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

using Base.LinAlg.LAPACK

# y = phenotype matrix
# X = predictor matrix
# K = kinship matrix, expected to be symmetric and positive definite

function rotateData(y,X,K)
    # check dimensions
    n = size(y,1)
    if( size(X,1) != n )
        error("Dimension mismatch.")
    elseif( size(K,1) != n )
        error("Dimension mismatch.")
    end

    # check symmetry and positive definiteness of K
    if( !issymm(K) )
        error("K is not symmetric.")
    elseif( !isposdef(K) )
        error("K is not positive definite.")
    end

    # spectral decomposition of a symmetric matrix
    # this rewrites the K matrix, so need to check what to do
    # allocate future eigenvector matrix to kinship matrix
    U = K
    S = syev!(jobz='N',uplo='U',U);

    # return rotated phenotype, covariates, and eigenvalues
    return U'y, U'X, S
    
end


