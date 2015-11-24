using Distributions
using Distances

# include the function
include("../src/lmm.jl")

# make a random kinship matrix

N = 100
M = 500
P = 2

# random genotypes
g = rand( Binomial(1,0.5), (N,M) );
# calculate pairwise distances
d = pairwise( Cityblock(),g' );
# calculate kinship
K = 1-d/(M);

e = rand( Normal(), (N,1) );
X = rand( Normal(), (N,2) );
beta =[2 3]';

(yy,XX,dd) =rotateData(y,X,K);
X = [ones(div(N,2),1);zeros(div(N,2),1)];
w = 1 + 99*X
e = e.*sqrt(w[:,1])
X = [ones(N,1) X];
y = X*beta + e;
