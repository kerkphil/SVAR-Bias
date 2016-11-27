% This set of commands generates data for Monte Carlo iterations given an
% initial estimate of the VAR coefficients (betaorig) and the A0orig matrix

% generate a series of uncorrelated shocks (as in the structural representation)
% use twice sample size
E = randn(2*T+P,Q);
% generate correlated shocks (like the residuals in the OLS reduced form)
U = E*A0orig;

% iniitalize data
% Here the data are arranged with observations in rows and variables in
% columns.
data = zeros(2*T,Q);
data(1,:) = dataorig(1,:);
for t=2:2*T
    for j=1:min(P,t-1)
        data(t,:) = data(t,:) + data(t-j,:)*Borig(:,:,j);
    end
    data(t,:) = data(t,:) + consorig + U(t,:);
end
data = data(T+1:2*T,:);