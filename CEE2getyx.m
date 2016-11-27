% This set of coommands takes a T by Q dataset and creates two matrices
% suitable for OLS estimation:  a T-P by Q y matrix and a T-P by P*Q+1 x
% matrix

%put current value into y vector
y = data(P+1:T,:);
% delcare x vector: Q*P no constants, Q*P+1 includes constants
x = ones(T-P,Q*P+1);
%put lagged values into x matrix
for p=1:P
    x(:,Q*(p-1)+1:Q*p) = data(P+1-p:T-p,:);
end