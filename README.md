# autocorrelation
Python script for auto-correlations in expression in a genome at different lags

Requires: Python 3

autocor.py
Input: File with 2 columns: first column containing coordinates, and second column containing values at co-ordinate.
Output: A plot showing correlation between pairs of values separated by all values of lag from 0 to MAX_LAG. 

autocor_median_circular.py
Input:  File with 2 columns: first column containing coordinates, and second column containing values at co-ordinate.
Output: File with autocorrelations in first column and lag in second column
