###########################################################
# genome scan function; no covariates, two genotype groups
###########################################################

using Distributed
using Random
using LinearAlgebra
using SharedArrays

include("lmm.jl")
include("util.jl")

###
# this function may not be needed
###

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
    wts = makeweights( out00.h2,lambda0 )
    # rescale by weights
    rowScale!(y0,sqrt.(wts))
    rowScale!(X0,sqrt.(wts))

    # perform genome scan
    out0 = rss(y0,reshape(X0[:,1],n,1))
    lod = zeros(m)
    X = zeros(n,2)
    X[:,1] = X0[:,1]
    for i = 1:m
        X[:,2] = X0[:,i+1]
        out1 = rss(y0,X)
        lod[i] = (n/2)*(log10(out0[1]) - log10(out1[1]))
    end

    return ( out00.sigma2, out00.h2, lod )

end


## genome scan with permutations

function scan(y::Array{Float64,2},g::Array{Float64,2},
              K::Array{Float64,2},nperm::Int64=1024,
              rndseed::Int64=0,reml::Bool=true)

    # number of markers
    (n,m) = size(g)
    # make intercept
    intcpt = ones(n,1)
    # rotate data
    (y0,X0,lambda0) = rotateData(y,[intcpt g],K)
    # fit null lmm
    vc = flmm(y0,reshape(X0[:,1], :, 1),lambda0,reml)
    # weights proportional to the variances
    wts = makeweights( vc.h2,lambda0 )
    # rescale by weights
    rowScale!(y0,sqrt.(wts))
    rowScale!(X0,sqrt.(wts))

    ## random permutations
    rng = MersenneTwister(rndseed);
    y0perm = shuffleVector(rng,y0[:,1],nperm)

    ## null rss
    out0 = rss(y0perm,reshape(X0[:,1],n,1))
    ## make shared array to hold LOD scores
    lod = SharedArray(zeros(nperm+1,m))
    ## initialize covariate matrix
    X = zeros(n,2)
    X[:,1] = X0[:,1]
    ## loop over markers
    for i = 1:m
        ## change the second column of covariate matrix X
        X[:,2] = X0[:,i+1]
        ## alternative rss
        out1 = rss(y0perm,X)
        ## calculate LOD score and assign
        lod[:,i] = (n/2)*(log10.(out0) .- log10.(out1))
    end

    return lod

end

function makeweights( h2::Float64,
                      lambda::Array{Float64,1} )
    return h2*lambda .+ (1.0-h2)
end


