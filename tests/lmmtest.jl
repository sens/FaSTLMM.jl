using Distributions
using Distances
using Optim

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

X = [ones(div(N,2),1);zeros(div(N,2),1)];
w = 1 + 99*X
e = e.*sqrt(w[:,1])
X = [ones(N,1) X];
y = X*beta + e;

(yy,XX,dd) =rotateData(y,X,K);
estVarComp(yy,XX,dd,log(var(y)),0.0)

##############

using DataFrames
using DataArrays

K = readtable("../data/kinship.csv");
K = K[2:size(K,2)];
K = DataArray(K);

pheno = readtable("../data/pheno.csv");
pheno = pheno[2:size(pheno,2)];
pheno = DataArray(pheno);

covar = readtable("../data/covar.csv");
covar = covar[2:size(covar,2)];
covar = DataArray(covar);

y = Array{Float64}(size(pheno,1),1);
z = convert(Array{Float64,1},pheno[:,1]);
y[:,1] = z;
X = convert(Array{Float64,2},covar);
K = convert(Array{Float64,2},K);

(yy,XX,dd) = rotateData(y,X,K)
