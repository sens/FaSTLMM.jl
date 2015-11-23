# include the function
include("../src/lmm.jl")

# make a random kinship matrix

N = 100
M = 50

# random genotypes
g = rand( Binomial(2,0.5), (N,M) );
# calculate pairwise distances
d = pairwise( Cityblock(),g );
# calculate kinship
K = 1-d;

