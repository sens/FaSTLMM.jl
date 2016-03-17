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

##################################################################
# rotateData: rotate by orthogonal transformation
##################################################################

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

#    # check symmetry and positive definiteness of K
#    if( !(issym(K)) )
#        error("K is not symmetric.")
#    end

#    if( !(isposdef(K)) )
#        error("K is not positive definite.")
#    end

    # spectral decomposition of a symmetric matrix
    EF = eigfact(K)

    # return rotated phenotype, covariates, and eigenvalues
    return EF[:vectors]'y, EF[:vectors]'X, EF[:values]

end

##################################################################
# wls: weighted least squares        
##################################################################        

type Wls
    b::Array{Float64,2}
    sigma2::Float64
    ell::Float64
end
            
            
"""
wls: Weighted least squares estimation

y = outcome, matrix
X = predictors, matrix
w = weights (positive, inversely proportional to variance), one-dim vector

The variance estimate is maximum likelihood
"""

function wls(y::Array{Float64,2},X::Array{Float64,2},w::Array{Float64,1},
             reml::Bool=false,loglik=false)

    # number of individuals
    n = size(y,1)
    # number of covariates
    p = size(X,2)
    
    # check if weights are positive
    if(any(w.<=.0))
        error("Some weights are not positive.")
    end
        
    # square root of the weights
    sqrtw = sqrt(w)
    # scale by weights
    yy = y.*sqrtw
    XX = diagm(sqrtw)*X

        
    # QR decomposition of the transformed data
    (q,r) = qr(XX)
    b = r\At_mul_B(q,yy)
    # estimate yy and calculate rss
    yyhat = XX*b
    rss = sum((yy-yyhat).^2)

    if( reml )        
        sigma2 = rss/(n-p)
    else
        sigma2 = rss/n
    end
        
    # return coefficient and variance estimate
    if(loglik)
        logdetSigma = n*log(sigma2) - sum(log(w))
        ell = -0.5 * ( logdetSigma + rss/sigma2 )
        if ( reml )
            ell += - log(det(r))
        end
        return Wls(b,sigma2,ell)
    else
        return b, sigma2
    end
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
        -logLik(z[1],z[2],y,X,d,false)
    end

    est = optimize(logLik0,[logsigma2,logith2])
    return exp(est.minimum[1]) , invlogit(est.minimum[2])
end

##################################################################
# function to fit mixed model
##################################################################

#=
function lmm( y::Array{Float64,2}, X::Array{Float64,2},
                    d::Array{Float64,1},logsigma2::Float64,logith2::Float64)
=#
