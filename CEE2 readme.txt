CEE2 programs

This is a set of MATLAB programs that replicates CEE(1999), allows for three different estimators of the VCV matrix,
and conducts a set of Monte Carlos to find coverage rates for chosen IRFs.


CEE2.m is a program that estimates and boostraps IRFs for the CEE model

The data used to estimate and generate IRFs is found in the file, CEEM1Qtrunc.xls in cells B2:H125.

CEE2getxy.m takes a TxQ data matrix called 'data' and sets up the X and Y matrices needed to estimate a VAR(P)

CEE2estim.m is a function that takes Y & X as data inputs and dfa as an indicator of which estimator to use for the VCV matrix, 
  it estimates the VAR(P) by OLS and outputs the coeffiecient estimates, the VCV matrix and the series of residuals.

CEE2irfsolve.m is a function that takes coefficient estimates and the VCV matrix as inputs and computes and outputs the
  structural IRFs.  A lower Cholesky decomposition is used.

CEE2plotirfs.m plots the structural IRF chosen


CEE2cover.m is a program that runs a series of Monte Carlos to determine coverage rates for a chosen IRF

CEE2generate.m creates a new dataset for the Monte Carlo iteration using the estimated reduced form VAR as the data generating
  process

CEE2estboot.m replicates CEE2.m by estimating and bootstrapping the IRFs for a Monte Carlo iteration

CEE2check.m determines if the original IRF lies above, below, or inside the confidence bands for a particular Monte Carlo
  iteration






