include("../src/scan.jl")
include("../src/kinship.jl")
include("../src/readData.jl")
include("../src/wls.jl")

using DelimitedFiles
using LinearAlgebra
using Optim
using Distributions
# using FaSTLMM

using BenchmarkTools

pheno_file = "../data/bxdData/traits.csv"
pheno = readBXDpheno(pheno_file)
geno_file = "../data/bxdData/geno_prob.csv"
geno = readGenoProb(geno_file)
k = calcKinship(geno)

## genome scan
time = @btime lod = scan(reshape(pheno[:,1], :, 1), geno, k, true)
## genome scan permutation
#@btime scan(reshape(pheno[:,1], :, 1), geno, k, 1024,1,true);

## transform LOD to -log10(p) (univariate)
julia_result = -log.(10,(ccdf.(Chisq(1),2*log(10)*lod)));
julia_result = julia_result[1:2:end]


#################################################################
#                              gemma                            #
#################################################################

## Converting data sets to format usable by gemma: 

geno_output_file = "../data/bxdData/bxd_geno_for_gemma.txt"
transform_bxd_geno_to_gemma(geno_file, geno_output_file);

pheno_output_file = "../data/bxdData/pheno_for_gemma.txt"
transform_bxd_pheno_to_gemma(pheno_file,pheno_output_file);

## run gemma:

gemma_bin = "../software/gemma-0.98.1-linux-static"
# Run this command in terminal to get kinship matrix from gemma. 
run(`$gemma_bin -g $geno_output_file -p $pheno_output_file -gk -no-check`)
# Run this command in terminal to get gemma result, scan_result is the output file.  
run(`$gemma_bin -g $geno_output_file -p $pheno_output_file -k ../software/output/result.cXX.txt -lmm 2 -o scan_result -no-check`)

##for gemma ouput (LRT) :-log10(p) transformation
#gemma=readdlm("lrt_chr1.assoc.txt";header=true)
#lrtp=gemma[1][:,end]
#-log.(10,lrtp)
gemma_scan = readdlm("../software/output/scan_result.assoc.txt";header=true)
lrtp=gemma_scan[1][:,end]
gemma_result = -log.(10,lrtp)


#################################################################
#                              compare                          #
#################################################################

@test julia_result ≈ gemma_result atol=1e-4
println("Compare result: $(isapprox(julia_result, gemma_result, atol=1e-4)")