using Distributions
using Distances
using Optim
using DataArrays
using DataFrames

# include the function
include("../src/lmm.jl")


K = readtable("../data/kinship.csv");
K = K[2:size(K,2)];
K = DataArray(K);

pheno = readtable("../data/pheno.csv");
pheno = pheno[2:size(pheno,2)];
pheno = DataArray(pheno);

covar = readtable("../data/covar.csv");
covar = covar[2:size(covar,2)];
covar = DataArray(covar);

X = convert(Array{Float64,2},covar);
K = convert(Array{Float64,2},K);

y = Array{Float64}(size(pheno,1),1);
y[:,1] = convert(Array{Float64,1},pheno[:,1]);
(yy,XX,dd) = rotateData(y,X,K)

y[:,1] = convert(Array{Float64,1},pheno[:,12]);
(yy,XX,dd) = rotateData(y,X,K)


# NPOINTS = 1000;
# loglik = Array{Float64}(NPOINTS);
# p = ((1:NPOINTS)-0.5)/NPOINTS;
# for i = 1:NPOINTS
#     loglik[i] = wls(yy,XX,1./(p[i]*dd+(1-p[i])),false,true).ell
# end


# plot(x=p,y=loglik,Geom.line)

@time for j=1:100
    for i = 1:size(pheno,2)
        whichKeep = !isna(pheno[:,i])
        y = Array{Float64}(sum(whichKeep),1)
        y[:,1] = convert(Array{Float64,1},pheno[whichKeep,i]);
        (yy,XX,lambda) = rotateData(y,X[whichKeep,:],
                                    K[whichKeep,whichKeep])
        out = lmm(yy,XX,lambda,true)
    end
end
