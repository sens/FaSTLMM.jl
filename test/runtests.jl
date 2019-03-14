using Test
using Optim
# using DataArrays
using DataFrames
# using FaSTLMM
using DelimitedFiles
using Missings 


# include the function
# include("../src/lmm.jl")


# K = readtable("../data/kinship.csv");
# K = K[2:size(K,2)];
# K = DataArray(K);

K = readdlm("../data/kinship.csv", ','; skipstart=1)[:, 2:end]

# pheno = readtable("../data/pheno.csv");
# pheno = pheno[2:size(pheno,2)];
# pheno = DataArray(pheno);

pheno = readdlm("../data/pheno.csv", ','; skipstart=1)[:, 2:end]

# covar = readtable("../data/covar.csv");
# covar = covar[2:size(covar,2)];
# covar = DataArray(covar);

X = readdlm("../data/pheno.csv", ','; skipstart=1)[:, 2:end]

# X = convert(Array{Float64,2},covar);
# K = convert(Array{Float64,2},K);

# initialize matrix to store results
res = Array{Float64}(size(pheno,2)*2,size(covar,2)+4)

# loop through the phenotypes
for i = 1:size(pheno,2)
    # keep only those individuals without missing phenotypes
    whichKeep = .!(ismissing.(pheno[:,i]))
    y = Array{Float64}(undef, sum(whichKeep),1)
    y[:,1] = convert(Array{Float64,1},pheno[whichKeep,i]);
    # perform rotation
    (yy,XX,lambda) = rotateData(y,X[whichKeep,:],
                                K[whichKeep,whichKeep])
    out0 = flmm(yy,XX,lambda,false)
    out1 = flmm(yy,XX,lambda,true)    
    res[2*i-1,:] = [out0.b; out0.sigma2; out0.h2; out0.ell; 0]
    res[2*i,:]   = [out1.b; out1.sigma2; out1.h2; out1.ell; 1]    
end

cnames =["b0";"b1";"sigma2";"h2";"loglik";"reml"];
resDF = DataFrame(res);
names!(resDF,convert(Array{Symbol},cnames));
writetable("julia_results.csv",resDF);

pylmm = readtable("../data/pylmm_results.csv");

# test for ml
@test (pylmm[:hsq])[2*(1:26)] ≈ (resDF[:h2])[2*(1:26)-1] atol=0.1
@test maximum(abs(log(resDF[:sigma2][2*(1:26)-1] ./
                      pylmm[:sigmasq][2*(1:26)]))) <= 0.012

# test for ml
@test  pylmm[:hsq][2*(1:26)-1] ≈ resDF[:h2][2*(1:26)] atol=0.1
@test maximum(abs(log(resDF[:sigma2][2*(1:26)] ./
                      pylmm[:sigmasq][2*(1:26)-1]))) <= 0.01
