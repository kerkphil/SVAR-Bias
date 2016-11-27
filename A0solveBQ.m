function y = A0solveBQ(A0)
global Sig B C A Acum K
% This function takes A0 as a direct input and the estimated B(L)
% polynomial & Sig as global inputs.  It outputs the difference between
% A0'A0 and Sig.  It is used in the main program bootBQ.m in an fsolve
% command to numerically solve for the value of A0.

% A is a 2x2 matrix
% B is a 2x2xP array of estimated coefficients
% Sig is a 2x2 matrix of reduced form variance/covariances

[row,col,P]=size(B);
% declare C matrices
C = zeros(2,2,K);
C(:,:,1) = [1 0; 0 1];
for k=2:K
    for p=2:P
        if k-p+1>0
            C(:,:,k) = C(:,:,k) - B(:,:,p)*C(:,:,k-p+1);
        end
    end
end
% declare A matrices
A = zeros(2,2,K);
Acum = A;
%calculate A(L)
for k=1:K
    A(:,:,k) = A0*C(:,:,k);
end
%sum up A(L) for effects on GDP
Acum(:,:,1) = A(:,:,1);
for k=2:K
    Acum(:,:,k) = Acum(:,:,k-1) + A(:,:,k);
end
% solve for LR cumulative effect as sum goes to infinity
B1 = zeros(2,2);
for p=1:P
    B1 = B1 + B(:,:,p);
end
ALR = (B1^-1)*A0;

%check for consistent direction of IRFs
% sum the effects on GDP
crit1 = Acum(1,1,7);
crit2 = Acum(2,1,K);
%crit2 = Acum(2,1,K);
% for k=2:K
%     %crit1 = crit1 + Acum(1,1,k);
%     crit2 = crit2 + Acum(2,1,k);
% end
% change sign if necessary
if crit1 < 0
    A0(1,1) = -1*A0(1,1);
    A0(1,2) = -1*A0(1,2);
    for k=1:K
        A(1,1,k) = -1*A(1,1,k);
        A(1,2,k) = -1*A(1,2,k);
        Acum(1,1,k) = -1*Acum(1,1,k);
        Acum(1,2,k) = -1*Acum(1,2,k);
    end
end
if crit2 < 0
    A0(2,1) = -1*A0(2,1);
    A0(2,2) = -1*A0(2,2);
    for k=1:K
        A(2,1,k) = -1*A(2,1,k);
        A(2,2,k) = -1*A(2,2,k);
        Acum(2,1,k) = -1*Acum(2,1,k);
        Acum(2,2,k) = -1*Acum(2,2,k);
    end
end

dev = A0'*A0 - Sig;
y = [dev(1,1) dev(1,2) dev(2,2) ALR(1,1)];