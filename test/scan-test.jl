include("../src/scan.jl")
include("../src/kinship.jl")
include("../src/readData.jl")
include("../src/wls.jl")

using DelimitedFiles
using LinearAlgebra
using Optim
# using FaSTLMM

using BenchmarkTools

pheno = readBXDpheno("../data/bxdData/traits.csv")
geno = readGenoProb("../data/bxdData/geno_prob.csv")
k = calcKinship(geno)

## genome scan
@btime scan(reshape(pheno[:,1], :, 1), geno, k, true)
## genome scan permutation
@btime scan(reshape(pheno[:,1], :, 1), geno, k, 1024,1,true);
