#################################
#Author: Shweta Ramdas
#necessary argument: --infile INFILE: file name with coordinate and expression alue
#necessary argument: --outfileOUTFILE: output file suffix
#Optional arguments:
#--spearman: to compate spearman correlation
#--min_lag : change the minimum lag
#--max_lag : change the maximum lag
#Output: Single file with autocorrelation (pearson) in first column and coordinate in second column
#################################

import sys
import numpy as np
from pandas import Series
from matplotlib import pyplot as plt
from scipy.stats.stats import pearsonr 
from scipy.stats.stats import spearmanr
import datetime
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
import statistics
import argparse

#input file, formatted as having 2 columns: position, and value
parser = argparse.ArgumentParser(description='Calculate autocorrelation...')
parser.add_argument("--infile")
parser.add_argument("--outfile")
parser.add_argument("--min_lag")
parser.add_argument("--max_lag")
parser.add_argument("--spearman",action="store_true")
args = parser.parse_args()


filename = args.infile
outname = args.outfile
txn = pd.read_table(filename,header=None)
txn = txn.drop_duplicates()
txn.columns = ['coord', 'value']

unique_coords = txn.coord.drop_duplicates()

txn2 = txn.groupby(txn.coord)[['value']].median()

#for parallel processing
num_cores = multiprocessing.cpu_count()

#lags to compute
MIN_LAG = 1
MAX_LAG = 60000

if args.min_lag:
	MIN_lAG = int(args.min_lag)
if args.max_lag:
	MAX_LAG = int(args.max_lag)

LENGTH_GENOME = 4641652

#For each lag, this for loop identifies pairs of loci in the genome with that lag, and saves their expression values in a list. Then it computes the correlations between the 2 lists.
def get_pearson_corr(lag):
#	print(lag)
#	print(datetime.datetime.now())
	exp1 = []
	exp2 = []
	for i in txn2.index:
		e1 = txn2.value[i]
		if (i + lag) > LENGTH_GENOME:
#			e2 = txn2.value[i + lag - LENGTH_GENOME]
			newind = i + lag - LENGTH_GENOME
		else:
#			e2 = txn2.value[i + lag]
			newind = i + lag

		if newind in txn2.index:
			e2 = txn2.value[newind]
			exp1.append(e1)
			exp2.append(e2)
	return [lag, pearsonr(exp1,exp2)[0]]

def get_spearman_corr(lag):
#       print(lag)
#       print(datetime.datetime.now())
        exp1 = []
        exp2 = []
        for i in txn2.index:
                e1 = txn2.value[i]
                if (i + lag) > LENGTH_GENOME:
#                       e2 = txn2.value[i + lag - LENGTH_GENOME]
                        newind = i + lag - LENGTH_GENOME
                else:
#                       e2 = txn2.value[i + lag]
                        newind = i + lag

                if newind in txn2.index:
                        e2 = txn2.value[newind]
                        exp1.append(e1)
                        exp2.append(e2)
        return [lag, spearmanr(exp1,exp2)[0]]


print(datetime.datetime.now())
if args.spearman:
	RES = Parallel(n_jobs=num_cores)(delayed(get_spearman_corr)(i) for i in range(MIN_LAG,MAX_LAG))
else:
	RES = Parallel(n_jobs=num_cores)(delayed(get_pearson_corr)(i) for i in range(MIN_LAG,MAX_LAG))

print(datetime.datetime.now())

lags = []
cors = []
for i in range(0, len(RES)):
	lags.append(RES[i][0])
	cors.append(RES[i][1])

#Display auto-correlation plot
#plt.scatter(lags, cors)
#plt.show()

df = pd.DataFrame(lags, cors)
df.to_csv(str(MIN_LAG) + "_" + str(MAX_LAG) + "_" + outname,sep="\t")
