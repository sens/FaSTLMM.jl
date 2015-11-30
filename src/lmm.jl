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

function rotateData(y::Array{Float64,2},X::Array{Float64,2},
                    K::Array{Float64,2})

    # check dimensions
    n = size(y,1)
    if( ( size(X,1) != n ) | ( size(K,1) != n ))
        error("Dimension mismatch.")
    end

    # check symmetry and positive definiteness of K
    if( !(issym(K)) )
        error("K is not symmetric.")
    end

    if( !(isposdef(K)) )
        error("K is not positive definite.")
    end

    # spectral decomposition of a symmetric matrix
    EF = eigfact(K)

    # return rotated phenotype, covariates, and eigenvalues
    return EF[:vectors]'y, EF[:vectors]'X, EF[:values]

end



"""
wls: Weighted least squares estimation

    y = outcome, matrix
    X = predictors, matrix
    w = weights (should be positive), one-dim vector

The variance estimate is maximum likelihood
"""

function wls(y::Array{Float64,2},X::Array{Float64,2},w::Array{Float64,1})

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
    yhat = X*b
    rss = sum((diagm(1./sqrtw)*(y-yhat)).^2,1)

    # return coefficient and variance estimate
    return b, rss

end

function invlogit(x::Float64)
    exp(x)/(1+exp(x))
end
    
# ################################################################
# function to calculate log likelihood of data given fixed effects
# ################################################################
"""
logLik: log likelihood of data
"""
function logLik(logsigma2::Float64,logith2::Float64,
                y::Array{Float64,2},
                X::Array{Float64,2},
                d::Array{Float64,1})
    # weights
    w =exp(logsigma2) * ( invlogit(logith2)*d + (1-invlogit(logith2)) )

    (b,rss) = wls(y,X,w)
    yhat = X*b
    sqrtw = sqrt(w)

    n = size(w,1)
    # get normal pdfs
    lp = -rss/2 - 0.5 * sum(log(w)) 
    # sum to get log likelihood
    return lp[1,1]
end

##################################################################
# function to estimate variance components and heritability
##################################################################

"""
estVarComp: estimate variance components
"""


function estVarComp(y::Array{Float64,2},
                    X::Array{Float64,2},
                    d::Array{Float64,1},logsigma2::Float64,logith2::Float64)

    function logLik0(z::Array{Float64,1})
        -logLik(z[1],z[2],y,X,d)
    end

    est = optimize(logLik0,[logsigma2,logith2])
    return exp(est.minimum[1]) , invlogit(est.minimum[2])
end

##################################################################
# function to fit mixed model
##################################################################

function lmm( y::Array{Float64,2}, X::Array{Float64,2},
                    d::Array{Float64,1},logsigma2::Float64,logith2::Float64)

