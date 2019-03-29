###########################################################
# genome scan function; no covariates, two genotype groups
###########################################################

using Distributed
using Random
using LinearAlgebra
using SharedArrays

include("lmm.jl")


function scan(y::Array{Float64,2},g::Array{Float64,2},
              K::Array{Float64,2},reml::Bool)

    # number of markers
    (n,m) = size(g)
    # make intercept
    intcpt = ones(n,1)
    # rotate data
    (y0,X0,lambda0) = rotateData(y,[intcpt g],K)
    # fit null lmm
    out00 = flmm(y0,reshape(X0[:,1], :, 1),lambda0,reml)
    # weights proportional to the variances
    wts = makeweights( out00.sigma2,out00.h2,lambda0 )
    # rescale by weights
    # scale!(sqrt.(1 ./wts),y0)
    # scale!(sqrt.(1 ./wts),X0)
    weight = sqrt.(1 ./wts)
    scale = diagm(0 => weight)
    y0 = scale * y0
    X0 = scale * X0

    # perform genome scan
    # rss0 = sum(y0.^2)
    out0 = rss(y0,reshape(X0[:,1],n,1))
    lod = zeros(m)
    X = zeros(n,2)
    X[:,1] = X0[:,1]
    for i = 1:m 
        X[:,2] = X0[:,i+1]
        out1 = rss(y0,X)
        lod[i] = (n/2)*(log10(out0[1]) - log10(out1[1]))
    end

    return lod
        
end


## genome scan with permutations

function scan(y::Array{Float64,2},g::Array{Float64,2},
              K::Array{Float64,2},nperm::Int64=1024,rndseed::Int64=0,reml::Bool=true)

    # number of markers
    (n,m) = size(g)
    # make intercept
    intcpt = ones(n,1)
    # rotate data
    (y0,X0,lambda0) = rotateData(y,[intcpt g],K)
    # fit null lmm
    out00 = flmm(y0,reshape(X0[:,1], :, 1),lambda0,reml)
    # weights proportional to the variances
    wts = makeweights( out00.sigma2,out00.h2,lambda0 )
    # rescale by weights
    # scale!(sqrt.(1 ./wts),y0)
    # scale!(sqrt.(1 ./wts),X0)
    weight = sqrt.(1 ./wts)
    scale = diagm(0 => weight)
    y0 = scale * y0
    X0 = scale * X0

    ## random permutations
    rng = MersenneTwister(rndseed);
    y0perm = zeros(n,nperm+1)
    y0perm[:,1] = y0
    for i=1:nperm
        y0perm[:,i+1] = shuffle(rng,y0)
    end
    
    # perform genome scan
    out0 = rss(y0perm,reshape(X0[:,1],n,1))
    lod = SharedArray(zeros(nperm+1,m))
    X = zeros(n,2)
    X[:,1] = X0[:,1]
    @sync @distributed for i = 1:m 
        X[:,2] = X0[:,i+1]
        out1 = rss(y0perm,X)
        lod[:,i] = (n/2)*(log10.(out0) .- log10.(out1))
    end

    return lod
        
end

function makeweights( sigma2::Float64, h2::Float64,
                      lambda::Array{Float64,1} )
    return h2*lambda .+ (1.0-h2)
end
