% This set of commands pulls out the desired IRFs from the original OLS
% regression, the averages of the bootstrapped IRFs and the upper and lower
% confidence bands and then plots them.

% The original estimate is bold, the average of the boostraps is a normal
% solid line, and the confidence bands are dashed.

IRF = reshape(IRFmat1(variable,shock,:),K+1,1);
IRFavg = reshape(IRFavgmat(variable,shock,:),K+1,1);
IRFupp = reshape(IRFuppmat(variable,shock,:),K+1,1);
IRFlow = reshape(IRFlowmat(variable,shock,:),K+1,1);
% PLOT IRFS
CEE2irf([IRF IRFavg IRFupp IRFlow])