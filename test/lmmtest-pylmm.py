from pylmm import lmm
import numpy
import time

n_times = 100

# load the datasets
kinship = numpy.genfromtxt('kinship.csv', delimiter=',')
pheno =   numpy.genfromtxt('pheno.csv', delimiter=',')
covar =   numpy.genfromtxt('covar.csv', delimiter=',')

# loop over the columns in pheno
#    fit the LMM by REML or ML and print the results to STDIN
start_time = time.time()
for j in range(n_times):
    for i in range(pheno.shape[1]):
        result = lmm.LMM(pheno[:,[i]], kinship, X0=covar).fit(REML=True)
stop_time = time.time()
ave_time = (stop_time - start_time)/n_times
open('pylmm_time.txt', 'w').write(str(ave_time) + "\n")
