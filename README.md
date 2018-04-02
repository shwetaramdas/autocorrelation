# autocorrelation
Python script for computing auto-correlations in expression in a genome at different lags

To cite:
[![DOI](https://zenodo.org/badge/111421140.svg)](https://zenodo.org/badge/latestdoi/111421140)



Requires: Python 3, python package statistics, pandas, numpy

###autocor.py

Input: File with 2 columns: first column containing coordinates, and second column containing values at co-ordinate.
Output: A plot showing correlation between pairs of values separated by all values of lag from 0 to MAX_LAG. 

###autocor_median_circular.py

Example Usage: python autocor_median_circular.py --infile input_file_name.txt --outfile filecors.txt --min_lag 1 --max_lag 1000 --spearman

Required Arguments:
--infile
--outfile

Optional arguments:
--min_lag
--max_lag
--spearman

Input file:  File with 2 columns: first column containing coordinates, and second column containing values at co-ordinate. We assume that this is a circular genome, and the program wraps around the end position while calculating correlations for any particular lag value.

Output: File with autocorrelations in first column and lag in second column
