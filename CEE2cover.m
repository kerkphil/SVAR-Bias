% This program finds coverage rates for chosen IRFs from estimates of 
% a QxQ reduced form VAR of order P
% The particular application is to CEE

clear

global P Q K S

%DECLARE PARAMETERS
P = 4; % P is order of AR process to be estimated
K = 16; % K in number of IRF lags to be reported
N = 200; % N is number of bootstrap iterations
M = 1000; % M is number of Monte Carlos
confid = .05; % set confidence level
lower = floor(confid/2*N)+1;  % the integer value of lower confidence level in N
upper = ceil((1-confid/2)*N); % the integer value of upper confidence level in N
Origdf = 2; % use OLS estimator to find DGP
Startmethod = 1;  %1=data observations, 2=random draws
shock = 4;  %shock to plot (4)
variable = 1; %variable to plot (1)

%READ DATA FROM EXTERNAL XLS FILE
data = xlsread('C:MATLABwork\CEE\CEEM1Qtrunc2.xls','B2:H125');
% find sample size (T) and number of variables (Q) from sample
[T,Q]=size(data)
S = T-P; % S is the usable sample size
% GET y & x MATRICES
CEE2getyx
%ESTIMATE REDUCED FORM VAR
[beta,Sig,u] = CEE2estim(y,x,Origdf);
%SOLVE FOR IRFS
IRFmat = CEE2irfsolve(beta,Sig);
%SAVE RESULTS
dataorig = data;
betaorig = beta;
Sigorig = Sig;
% Sigorig(1,2) = 0;
% Sigorig(2,1) = 0;
IRForigmat = IRFmat;
IRForig = reshape(IRForigmat(variable,shock,:),K+1,1);

%FIND VALUES NEEDED TO GENERATE MONTE CARLO DATASETS
% find B(L)
Borig = zeros(Q,Q,P);
for p=1:P
    % no need to transpose here since we will be using this to generate
    % data with observations as row vectors.
    Borig(:,:,p) = betaorig(Q*(p-1)+1:Q*p,:);
end
% find regression constants
consorig = betaorig(P*Q+1,:);
% find A0 matrix
A0orig = chol(Sigorig,'lower')';

% choose estimators for monte carlos
Origdf = 1;  %1=MLE, 2=DFA1
Bootdf = 1;  %1=MLE, 2=DFA1, 3=DFA3  see CEE2estim.m for formulas

%MONTE CARLO LOOP
Coveragemc = zeros(K+1,1);
Abovemc = zeros(K+1,1);
Belowmc = zeros(K+1,1);
for m=1:M
    m
    %GENERATE NEW DATA
    CEE2generate
    %RUN REGRESSIONS & BOOTSTRAPS ON NEW DATA
    CEE2estboot
    %CHECK FOR COVERAGE & UPDATE MONTE CARLO AVERAGES
    CEE2check
    for k=1:K+1
        Coveragemc(k,1) = (m-1)*Coveragemc(k,1)/m + Coverage(k,1)/m;
        Abovemc(k,1) = (m-1)*Abovemc(k,1)/m + Above(k,1)/m;
        Belowmc(k,1) = (m-1)*Belowmc(k,1)/m + Below(k,1)/m;
    end
end
%REPORT
[Coveragemc Abovemc Belowmc]
figure;
plot(Coveragemc);
axis([0 K+1 0 1])