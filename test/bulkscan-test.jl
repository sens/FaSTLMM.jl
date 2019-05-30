using BenchmarkTools

include("../src/bulkscan.jl")

# number of individuals
N = 100
# number of markers
M = 7000
# number of phenotypes/traits
P = 20000

# phenotypes and genotypes
y = randn(N,P);
g = randn(N,M);

@btime lod = bulkscan(y,g);
