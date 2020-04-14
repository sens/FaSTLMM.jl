using Optim
using DelimitedFiles
using Missings
using LinearAlgebra

# include the function
include(joinpath(@__DIR__,"..","src","lmm.jl"))
include(joinpath(@__DIR__,"..","src","wls.jl"))

K = readdlm(joinpath(@__DIR__,"..","data","kinship.csv"), ','; skipstart=1)[:, 2:end]

pheno = readdlm(joinpath(@__DIR__, "..","data", "pheno.csv"), ','; skipstart=1)[:, 2:end]

covar = readdlm(joinpath(@__DIR__, "..","data", "covar.csv"), ','; skipstart=1)[:, 2:end]

X = convert(Array{Float64,2},covar);
K = convert(Array{Float64,2},K);

# initialize matrix to store results
res = Array{Float64}(undef,size(pheno,2)*2,size(covar,2)+4);

# loop through the phenotypes
for i = 1:size(pheno,2)
    # println("Column $i in $(size(pheno,2))")
    # replace all "NA" with missing type
    for j = 1:size(pheno,1)
        if pheno[j,i] == "NA"
            pheno[j,i] = missing
        end
    end

    # keep only those individuals without missing phenotypes
    whichKeep = .!ismissing.(pheno[:,i])
    y = Array{Float64}(undef,sum(whichKeep),1)
    y[:,1] = convert(Array{Float64,1},pheno[whichKeep,i]);


    # perform rotation
    (yy,XX,lambda) = rotateData(y,X[whichKeep,:],
                                K[whichKeep,whichKeep])
    out0 = flmm(yy,XX,lambda,false)
    out1 = flmm(yy,XX,lambda,true)
    res[2*i-1,:] = [out0.b; out0.sigma2; out0.h2; out0.ell; 0]
    res[2*i,:]   = [out1.b; out1.sigma2; out1.h2; out1.ell; 1]
end

cnames =["b0" "b1" "sigma2" "h2" "loglik" "reml"];
writedlm("julia_results.csv", [cnames; res], ",")

###################################################################

function benchmark(nrep::Int64,f::Function,x...;results::Bool=false)

    res = Array{Float64}(undef, nrep)

    for i=1:nrep
        start = time_ns()
        f(x...)
        res[i] = time_ns() - start
    end

    if(results)
        return res
    else
        return  [minimum(res) quantile(res,[0.25  0.5 0.75]) maximum(res)]
    end
end


function analyzeAllPheno(pheno::Array{Any,2},X::Array{Float64,2},
                         K::Array{Float64,2})
    for i = 1:size(pheno,2)
    #     # replace all "NA" with missing typ. UPDATE: Does not require replacing anymore because it is done already.
        # for j = 1:size(pheno,1)
        #     if pheno[j,i] == "NA"
        #         pheno[j,i] = missing
        #     end
        # end

        # keep only those individuals without missing phenotypes
        whichKeep = .!ismissing.(pheno[:,i])
        y = Array{Float64}(undef, sum(whichKeep),1)
        y[:,1] = convert(Array{Float64,1},pheno[whichKeep,i]);
        # perform rotation
        (yy,XX,lambda) = rotateData(y,X[whichKeep,:],
                                    K[whichKeep,whichKeep])
        out = flmm(yy,XX,lambda,true)
    end
end

###################################################################

res = benchmark(100,analyzeAllPheno,pheno,X,K,results=true)
writedlm("julia_time.csv",res, ",")
