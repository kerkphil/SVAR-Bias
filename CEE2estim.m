function [beta,Sig,u] = CEEestim(y,x,dfa)
% This function takes an S by Q dependent variable vector (y) and 
% a S by P*Q+1 vector of lagged variables (x).  The degrees of freedom
% adjustment for the estimate of Sig is given by an indicator variable
% (dfa).  See below for the actual formulas.
% It ouputs the P*Q+1 by Q OLS coefficient matrix (beta), 
% the variance-covariance matrix for the residuals (Sig),
% and the P*Q+1 by Q matrix of residuals (u)

global P Q S

beta = (x'*x)^(-1)*x'*y;
u = y - x*beta;
if dfa == 1
    Sig = u'*u./S; %MLE estimate
elseif dfa == 2
    Sig = u'*u./(S-P*Q-1);  %DFA1 estimate (biased for bootstrap)
else
    Sig = S.*(u'*u)./(S-P*Q-1)^2;  %DFA1 estimate (unbiased for bootstrap)
end