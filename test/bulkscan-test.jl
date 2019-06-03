using BenchmarkTools
using UnicodePlots

include("../src/bulkscan.jl")
include("../src/readData.jl")

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

pheno_file = "../data/bxdData/traits.csv"
pheno = readBXDpheno(pheno_file)
geno_file = "../data/bxdData/geno_prob.csv"
geno = readGenoProb(geno_file)

geno_output_file = "../data/bxdData/bxd_geno_for_gemma.txt"
transform_bxd_geno_to_gemma(geno_file, geno_output_file);

## run gemma:
# Run this command in terminal to get kinship matrix from gemma. 
run(`gemma -g $geno_output_file -p ../data/bxdData/pheno_for_gemma.txt -gk -no-check`)

# getting kinship matrix from gemma 
K = convert(Array{Float64,2},readdlm("./output/result.cXX.txt", '\t'));

intercept = ones(size(pheno,1),1);
(y0,X0,lambda) = rotateData(pheno,[intercept geno],K);

h2 = bulkesth2(y0,reshape(X0[:,1],:,1),lambda);
histogram(h2)
