###########################################################
# univariate genome scan functions for many traits
# no covariates, two genotype groups, no missing data
###########################################################



################################################################
# plain genome scan
################################################################

function bulkscan(y::Array{Float64,2},g::Array{Float64,2})

    y0 = deepcopy(y)
    g0 = deepcopy(g)

    colStandardize!(y0)
    colStandardize!(g0)

    lod = calcLod(y0,g0)
    return lod
end

# bulk weighted linear regression
# weights inversely proportional to variance

function bulkscan(y::Array{Float64,2},g::Array{Float64,2},w::Vector{Float64})

    n = size(y,1)
    p = size(y,2)
    m = size(g0,2)

    X0 = ones(n,1)
    y0=deepcopy(y)
    g0 = deepcopy(g)

    # convert problem to least squares
    sqrtw = sqrt.(w)
    rowScale!(X0,1.0/.sqrtw)
    rowScale!(y0,1.0/.sqrtw)
    rowScale!(g0,1.0/.sqrtw)

    y0 = bulkresid(y0,X0)
    g0 = bulkresid(g0,X0)

    # baseline sum of squares
    rss0 = sum(y0^.2,dims=1)

    # normalize genotypes
    g0ss = sum(g0^.2,dims=1)
    # rowScale!(g0,g0norm)

    b = (y0'*g0)

    lod = zeros(p,m)

    for i=1:p
        for j=1:m
            lod[i,j] = (n/2)*log10(b[i,j]^2/(rss0[i]*g0ss[j]))
        end
    end

    return lod
end


function bulkls(y::Matrix{Float64},X::Matrix{Float64},loglik::Bool=false)

    n = size(y,1)

    fct = qr(X)
    b = fct\y
    yhat = X*b
    resid = y-yhat
    rss = sum(resid.^2,dims=1)

    sigma2 = rss./n
    logdetSigma = n*.log.(sigma2)
    ell = -0.5 *. (logdetSigma +. rss./sigma2)

    return b, sigma2, ell

end

# residual from least squares

function bulkresid(y::Matrix{Float64},X::Matrix{Float64})

    n = size(y,1)

    fct = qr(X)
    b = fct\y
    yhat = X*b
    resid = y-yhat

    return resid

end

function bulkWls(y::Matrix{Float64},X::Matrix{Float64},
    w::Vector{Float64},loglik::Bool=false)

    # check if weights are positive
    if(any(w.<=.0))
        error("Some weights are not positive.")
    end

    sqrtw = sqrt.(w)
    y0 = rowScale!(deepcopy(y),1.0/.sqrtw))
    X0 = rowScale!(deepcopy(X),1.0/.sqrtw))

    (b,sigma2,ell) = bulkls(y0,X0,loglik)
     ell = ell + sum(log.(w))/2

    return b,sigma2,ell
end

# estimate heritability

function esth2(y::Matrix{Float64},X::Matrix{Float64},
               lambda::Vector{Float64})
    est = flmm(y,X,lambda,false)
    return est.h2
end

function esth2_grid(y::Matrix{Float64},X::Matrix{Float64},
                    lambda::Vector{Float64},ngrid::Int64)

    h2grid = (0:(ngrid-1)+0.5)/(ngrid)
    loglik = zeros(ngrid)


    for i=1:ngrid
        wts = makeweights(h2grid[i],lambda)
        # rowScale!(y0,sqrt.(wts))
        # rowScale!(X0,sqrt.(wts))
        out = wls(y,X,1.0./wts,false,true)
        loglik[i] = out.ell
    end
    return h2grid[argmax(loglik)]
end


# estimate heritability for multiple phenotypes

function bulkesth2(y::Matrix{Float64},X::Matrix{Float64},
                   lambda::Vector{Float64})
    p = size(y,2)
    h2 = zeros(p)
    Threads.@threads for i=1:p
    # for i=1:p
        h2[i] = esth2(reshape(y[:,i],:,1),X,lambda)
    end
    return h2
end

function bulkesth2_grid(y::Matrix{Float64},X::Matrix{Float64},
                   lambda::Vector{Float64},ngrid::Int64)
    p = size(y,2)
    h2 = zeros(p)
    Threads.@threads for i=1:p
    # for i=1:p
        h2[i] = esth2_grid(reshape(y[:,i],:,1),X,lambda,ngrid)
    end
    return h2
end

# calculate LOD scores

function calcLod(y::Matrix{Float64},g::Matrix{Float64})

    n = size(y,1)
    r = (y'*g) ./ n
    # r = LinearAlgebra.BLAS.gemm('T','N',y,g) ./ n
    r2lod!(r,n)
    return r

end

# convert correlations to LOD scores

function r2lod!(r::Matrix{Float64},n::Int64)

    c = -n/(2*log(10))
    (nr,nc) = size(r)
    # lod = zeros(nr,nc)
    Threads.@threads for j=1:nc
        for i=1:nr
            r[i,j] = c*log1p(-r[i,j]^2)
        end
    end
    # return lod
end

function r2lod_dist!(r::Matrix{Float64},n::Int64)

    c = -n/(2*log(10))
    (nr,nc) = size(r)
    r0 = SharedArray(r)
    # lod = zeros(nr,nc)
    @distributed for j=1:nc
        for i=1:nr
            r0[i,j] = c*log1p(-r0[i,j]^2)
        end
    end
    r[:] = Matrix(r0)
    # return lod
end
