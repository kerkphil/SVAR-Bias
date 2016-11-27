function A = CEEirfsolve(beta,Sig)
global P Q K

% This functions takes the P*Q+1 by Q coefficient matrix (beta) and the Q
% by Q variance-covariance matrix (Sig) as inputs.  And outputs a Q by Q by
% K array of the A(L) lag polynomial coefficients.  

%Y(t)=B(L)Y(t)+u(t) is the reduced form
%Y(t)=C(L)u(t) is the MA representation
%Y(t)=A(L)e(t) is the structural VAR
% Y(t),u(t) & e(t) are column vectors
% beta comes from the OLS regression where Y(t) is a row vector

%Formulas are from notes by David E. Spencer

% find B(L)
B = zeros(Q,Q,P);
for p=1:P
    B(:,:,p) = beta(Q*(p-1)+1:Q*p,:)';
end
% find C(L)
C = zeros(Q,Q,K+1);
C(:,:,1) = eye(Q);
for k=1:K
    for i=1:min(P,k)
        C(:,:,1+k) = C(:,:,k+1) + B(:,:,i)*C(:,:,k+1-i);
    end
end
% find A0
A0 = chol(Sig,'lower');
% find A(L)
A = zeros(Q,Q,K+1);
A(:,:,1) = A0;
for k=1:K
    A(:,:,k+1) = C(:,:,k+1)*A0;
end