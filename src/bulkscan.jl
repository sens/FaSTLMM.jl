###########################################################
# univariate genome scan functions for many traits
# no covariates, two genotype groups, no missing data
###########################################################

using Distributed
using Random
using LinearAlgebra
using SharedArrays

include("scan.jl")
include("util.jl")

################################################################
# plain genome scan
################################################################

function bulkscan(y::Array{Float64,2},g::Array{Float64,2})

    # number of markers
    (n,m) = size(g)
    
    y0 = deepcopy(y)
    g0 = deepcopy(g)

    colStandardize!(y0)
    colStandardize!(g0)

    lod = calcLod(y0,g0)
    return lod
end

# calculate lod scores

function calcLod(y::Matrix{Float64},g::Matrix{Float64})

    n = size(y,1)
    r = (y'*g) ./ n
    # r = LinearAlgebra.BLAS.gemm('T','N',y,g) ./ n    
    r2lod!(r,n)
    return r
    
end

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
