This file describes the programs used to investigate estimators other than those associated with CEE (1999)


biasUNIVARIATE.m  performs Monte Carlos for a univariate deterministic OLS regression.  It compares various estimates of the 
  variance of the error terms.

biaBIVARIATE.m  performs a similar exercise for a bivarite SUR and compares the estimates of the 2x2 VCV matrix.

biasAR.m  performs Monte Carlos for a univariate AR regression and compares estimates of the variance.

biasVAR.m  performs Monte Carlos for a bivariate VAR regression and compares estimates of the 2x2 VCV matrix.

bootBQ.m  replicates Blanchard & Quah (1989) and plots IRFs

BQadjusted.wk1 is a data file containing the data for estimation in bootBQ.m

A0solveBQ.m  is a function that is used to solve for the A0 matrix which, in turn, generates the IRFs in bootBQ.m. the
  function imposes long-run restrictions and is used in conjunction with MATLABs fsolve command.