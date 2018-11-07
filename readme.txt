1. unzip the file 'GSE59739_DataTable.zip' in the save folder with the source codes
2. open the file 'GSE59739_DataTable.txt' using Excel and save it as 'GSE59739_DataTable.xlsx'
3. run 'run_OGFSC_demo.m' to perform a demo analysis 



OGFSC input&output parameters 

Input:

data          - a variable vs samples matrix
nBins         - number of bins, by default 60.
minBinSize    - minimeansDatam bin size, by default 100.
LR_p          - p-value threshold to identify valid linear regression model, by default 0.01.
alpha         - list of different threshold, i.e., the upper boundary of confidence interval, by default [0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999].
TW_threshold  -  the threshold to determine the number of eigenvalues on TW distribution, by default 0.0001.
plot_option   - plot the outputs of ODFSC, by default 0. 



Output:


OGFSC_idx  - the index of OGFSC selected genes
idx_output - the index of genes successfully modelled by MLM
cv2_threshold - the estimated CV2 threshold for each successfully modelled gene

