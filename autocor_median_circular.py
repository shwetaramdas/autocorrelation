#################################
Author: Shweta Ramdas
Input Argument 1: file name with coordinate and expression alue
Input Argument 2: output file suffix
Output: Single file with autocorrelation in first column and coordinate in second column
#################################

import sys
import numpy as np
from pandas import Series
from matplotlib import pyplot as plt
from scipy.stats.stats import pearsonr 
import datetime
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
import statistics

#input file, formatted as having 2 columns: position, and value
filename = sys.argv[1]
outname = sys.argv[2]
txn = pd.read_table(filename,header=None)
txn = txn.drop_duplicates()
txn.columns = ['coord', 'value']

unique_coords = txn.coord.drop_duplicates()
#txn2 = txn.set_index('coord')

txn2 = txn.groupby(txn.coord)[['value']].median()

#for i in unique_coords:
#	txn2.value[i] = statistics.median(txn.value[txn.coord == i])

#for parallel processing
num_cores = multiprocessing.cpu_count()

MIN_LAG = 20000
MAX_LAG = 40000

LENGTH_GENOME = 4641652

#For each lag, this for loop identifies pairs of loci in the genome with that lag, and saves their expression values in a list. Then it computes the correlations between the 2 lists.
def get_corr(lag):
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

print(datetime.datetime.now())
RES = Parallel(n_jobs=num_cores)(delayed(get_corr)(i) for i in range(MIN_LAG,MAX_LAG))
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
