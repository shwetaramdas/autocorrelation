import pandas as pd
import numpy as np
from pandas import Series
from matplotlib import pyplot as plt
from scipy.stats.stats import pearsonr 
import datetime

#input file, formatted as having 2 columns: position, and value
filename = sys.argv[1]
txn = pd.read_table(filename,header=None)
txn = txn.drop_duplicates()
txn.columns = ['coord', 'value']
txn = txn.set_index('coord')

MAX_LAG = 1000

LENGTH_GENOME = max(txn.index)

lags = []
cors = []
#For each lag, this for loop identifies pairs of loci in the genome with that lag, and saves their expression values in a list. Then it computes the correlations between the 2 lists.
for lag in range(1, 1000):
	print(lag)
	print(datetime.datetime.now())
	exp1 = []
	exp2 = []
	for i in txn.index:
		e1 = txn.value[i]
		if (i + lag) > LENGTH_GENOME:
			break
		if (i + lag) in txn.index:
			e2 = txn.value[i+lag]
	#		print(e2)
			if  hasattr(e1,"__len__"):
				e1 = e1.iloc[0]
			if  hasattr(e2,"__len__"):
				e2 = e2.iloc[0]
			exp1.append(e1)
			exp2.append(e2)
			
	lags.append(lag)
	cors.append(pearsonr(exp1,exp2)[0])

#Display auto-correlation plot
plt.scatter(lags, cors)
plt.show()
