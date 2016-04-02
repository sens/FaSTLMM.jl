####################################################################
# Pylmm code for comparing results and performance
####################################################################
# Modified version of Karl Broman's code
####################################################################

# import modeules
from pylmm import lmm
import numpy
import time

# load the datasets
kinship = numpy.genfromtxt('../data/kinship.csv', delimiter=',')
pheno =   numpy.genfromtxt('../data/pheno.csv', delimiter=',')
covar =   numpy.genfromtxt('../data/covar.csv', delimiter=',')

# remove first row and column to conform to pylmm expectations
pheno = pheno[1:pheno.shape[0],1:pheno.shape[1]]
kinship = kinship[1:kinship.shape[0],1:kinship.shape[1]]
covar = covar[1:covar.shape[0],1:covar.shape[1]]

####################################################################
# results
####################################################################

# open file for writing
f = open('pylmm_results.csv', 'w')

# function to print results
def print_result(result, index, method):
    result_arr = [result[0], result[1][0,0], result[1][1,0], result[2][0,0], result[3]]
    for j in range(len(result_arr)):
        result_arr[j] = str(result_arr[j])
    f.write(','.join([str(index), method] + result_arr) + '\n')

# print header row
f.write('index,method,hsq,intercept,sex,sigmasq,loglik\n')

# loop over the columns in pheno
#    fit the LMM by REML or ML and print the results to STDIN
for i in range(pheno.shape[1]):
    result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=True)
    print_result(result, i+1, 'reml')
    result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=False)
    print_result(result, i+1, 'ml')
f.close()

####################################################################
# performance
####################################################################

# number of times to run for benchmarking
n_times = 100

# open file for writing
f = open('pylmm_time.csv', 'w')

# loop over the columns in pheno
#    fit the LMM by REML or ML and print the results to STDIN
for j in range(n_times):
    start_time = time.time()
    for i in range(pheno.shape[1]):
        result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=True)
    stop_time = time.time()
    this_time = stop_time - start_time
    f.write(str(this_time)+'\n')
f.close()
