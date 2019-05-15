include("../src/scan.jl")
include("../src/kinship.jl")
include("../src/readData.jl")
include("../src/wls.jl")
include("../src/lmm.jl")

using DelimitedFiles
using LinearAlgebra
using Optim
using Distributions
using Test
# using FaSTLMM

using BenchmarkTools


function run_julia(pheno::Array{Float64,1}, geno::Array{Float64,2}, k::Array{Float64,2}, reml::Bool, method::String )
    ## genome scan
    (sigma2, h2, lod) = scan(reshape(pheno, :, 1), geno, k, reml, method)
    ## genome scan permutation
    # scan(reshape(pheno[:,1], :, 1), geno, k, 1024,1,true);

    ## transform LOD to -log10(p) (univariate)
    result = -log.(10,(ccdf.(Chisq(1),2*log(10)*lod)));
    result = result[1:2:end]
    return (result, sigma2, h2)
end

function run_gemma(pheno_output_file::AbstractString, geno_output_file::AbstractString, iter::Int64)
    ## Converting data sets to format usable by gemma: 
    transform_bxd_pheno_to_gemma(pheno_file,pheno_output_file, iter);

    # Run this command in terminal to get gemma result, scan_result is the output file.  
    run(`gemma -g $geno_output_file -p $pheno_output_file -k ./output/result.cXX.txt -lmm 2 -o scan_result_$iter -no-check`)

    ##for gemma ouput (LRT) :-log10(p) transformation
    gemma_scan = readdlm("./output/scan_result_$iter.assoc.txt";header=true)
    lrtp=gemma_scan[1][:,end]
    return -log.(10,lrtp)
end

pheno_file = "../data/bxdData/traits.csv"
pheno = readBXDpheno(pheno_file)
geno_file = "../data/bxdData/geno_prob.csv"
geno = readGenoProb(geno_file)
# k = calcKinship(geno)

geno_output_file = "../data/bxdData/bxd_geno_for_gemma.txt"
transform_bxd_geno_to_gemma(geno_file, geno_output_file);

## run gemma:
# Run this command in terminal to get kinship matrix from gemma. 
run(`gemma -g $geno_output_file -p ../data/bxdData/pheno_for_gemma.txt -gk -no-check`)

# getting kinship matrix from gemma 
k = convert(Array{Float64,2},readdlm("./output/result.cXX.txt", '\t'))

# testing result captures the comparison result between lmm and gemma. 
# One row is for one phenotype, it contains the # of agreement, # of exeed threshold, # agreed and exeed threshold, sigma2, h2
num_run = 100
run_count = 0
testing_result = Array{Float64}(undef, num_run,5)#size(pheno)[2], 5)

julia_gemma = Array{Float64}(undef, Int64(size(geno)[2]/2), 2)

rescan_pass_rate = Array{Float64}(undef, num_run,2)

#looping over all phenotype. 
for i in 81:81#1:num_run
    global run_count ## added this line because of the stupidity of julia scope 

    #################################################################
    #                              gemma                            #
    #################################################################
    pheno_output_file = "../data/bxdData/pheno_for_gemma_$i.txt"
    gemma_result = run_gemma(pheno_output_file, geno_output_file, i)

    #################################################################
    #                              julia                            #
    #################################################################
    julia_time = @elapsed (julia_result, sigma2, h2) = run_julia(pheno[:,i], geno, k, false, "alt")
    println("Julia scan ran $julia_time seconds. ")

    
    #################################################################
    #                              compare                          #
    #################################################################

    (total_pass, total_th, total_pass_th) = compareValues(julia_result, gemma_result, 1e-2, 2.0)
    #columb name of testingresult is: # of agreement, # of exeed threshold, # agreed and exeed threshold, sigma2, h2
    testing_result[i,:] = [total_pass, total_th, total_pass_th, sigma2, h2]

    if (total_pass_th/total_th < 1.0)
        #do rescan 
        run_count += 1
        rescan_time = @elapsed (rescan_result, sigma2, h2) = run_julia(pheno[:,i], geno, k, false, "alt")
        (total_pass, total_th, total_pass_th) = compareValues(rescan_result, gemma_result, 1e-2, 2.0)
        rescan_pass_rate[run_count, :] = [i, (total_pass_th/total_th)] 
        testing_result[i,:] = [total_pass, total_th, total_pass_th, sigma2, h2]
        println("total_th: $total_th")
        println("rescan time: $rescan_time")
    end

    # display(testing_result[i,:])
    # println("passrate is $(cv[1]/size(julia_gemma)[1])")
    # display(lod)
    # display(julia_result)
    # display(gemma_result)
    julia_gemma[:,1]=rescan_result
    julia_gemma[:,2]=gemma_result
end

# display(testing_result)
# display(rescan_pass_rate[1:20,:])
# display(julia_gemma)
# writeToFile(testing_result,"./result/testing_result.txt")
