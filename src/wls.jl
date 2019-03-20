##################################################################
# wls: weighted least squares        
##################################################################        

mutable struct Wls
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
             reml::Bool=false,loglik::Bool=false)

    # check if weights are positive
    if(any(w.<=.0))
        error("Some weights are not positive.")
    end
         
    # square root of the weights
    sqrtw = sqrt.(w)
    # scale by weights
    # yy = y.*sqrtw
    yy = Diagonal(sqrtw)*y
    # XX = diagm(sqrtw)*X
    XX = Diagonal(sqrtw)*X

    out = ls(yy,XX,reml,loglik)

    if(loglik)
        out.ell = out.ell + sum(log.(w))/2 
    end

    return out
        
end

function ls(y::Array{Float64,2},X::Array{Float64,2},
             reml::Bool=false,loglik=false)

    # number of individuals
    n = size(y,1)
    # number of covariates
    p = size(X,2)
    
    # least squares solution
    fct = qr(X)
    b = fct\y
    
    # estimate yy and calculate rss
    yhat = X*b
    # yhat = q*At_mul_B(q,yy)
    rss = norm((y-yhat))^2

    if( reml )        
        sigma2 = rss/(n-p)
    else
        sigma2 = rss/n
    end

    # return coefficient and variance estimate
    # logdetSigma = n*log(sigma2) - sum(log.(w))
    logdetSigma = n*log(sigma2)
    ell = -0.5 * ( logdetSigma + rss/sigma2 ) 
    if ( reml )
        ell -=  log(abs(det(fct.R))) - (p/2)*(log(sigma2))
    end
        
    return Wls(b,sigma2,ell)
        
end

