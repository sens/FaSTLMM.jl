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




# NPOINTS = 1000;
# loglik = Array{Float64}(NPOINTS);
# p = ((1:NPOINTS)-0.5)/NPOINTS;
# for i = 1:NPOINTS
#     loglik[i] = wls(yy,XX,1./(p[i]*dd+(1-p[i])),false,true).ell
# end


# plot(x=p,y=loglik,Geom.line)

# initialize matrix to store results
res = Array{Float64}(size(pheno,2)*2,size(covar,2)+4)

# loop through the phenotypes
for i = 1:size(pheno,2)
    # keep only those individuals without missing phenotypes
    whichKeep = !isna(pheno[:,i])
    y = Array{Float64}(sum(whichKeep),1)
    y[:,1] = convert(Array{Float64,1},pheno[whichKeep,i]);
    # perform rotation
    (yy,XX,lambda) = rotateData(y,X[whichKeep,:],
                                K[whichKeep,whichKeep])
    out0 = lmm(yy,XX,lambda,false)
    out1 = lmm(yy,XX,lambda,true)    
    res[2*i-1,:] = [out0.b; out0.sigma2; out0.h2; out0.ell; 0]
    res[2*i,:]   = [out1.b; out1.sigma2; out1.h2; out1.ell; 1]    
end

cnames =["b0";"b1";"sigma2";"h2";"loglik";"reml"];
resDF = DataFrame(res);
names!(resDF,convert(Array{Symbol},cnames));
writetable("results.csv",resDF);

###################################################################

function microbenchmark(nrep::Int64,f::Function,x...)
    res = Array{Float64}(nrep)

    for i=1:nrep
        tic()
        f(x...)
        res[i] = toq()
    end
    return res
end


function analyzeAllPheno(pheno::DataArray{Real,2},X::Array{Float64,2},
                         K::Array{Float64,2})
    for i = 1:size(pheno,2)
        # keep only those individuals without missing phenotypes
        whichKeep = !isna(pheno[:,i])
        y = Array{Float64}(sum(whichKeep),1)
        y[:,1] = convert(Array{Float64,1},pheno[whichKeep,i]);
        # perform rotation
        (yy,XX,lambda) = rotateData(y,X[whichKeep,:],
                                    K[whichKeep,whichKeep])
        out = lmm(yy,XX,lambda,true)
    end
end
