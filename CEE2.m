% This program estimates a QxQ reduced form VAR of order P
% It also solves for a SVAR respresentation of the estimated process
% And it bootstraps confidence bands about these structural IRFs
% The particular application is to CEE
clear

global P Q K S

%DECLARE PARAMETERS
P = 4; % P is order of AR process to be estimated
K = 16; % K in number of IRF lags to be reported
N = 1000; % N is number of bootstrap iterations
confid = .05; % set confidence level
lower = floor(confid/2*N)+1;  % the integer value of lower confidence level in N
upper = ceil((1-confid/2)*N); % the integer value of upper confidence level in N
Origdf = 2;  %1=MLE, 2=DFA1
Bootdf = 3;  %1=MLE, 2=DFA1, 3=DFA3  see CEE2estim.m for formulas
Startmethod = 1;  %1=data observations, 2=random draws
shock = 4;  %shock to plot (4)
variable = 1; %variable to plot (1)


%READ DATA FROM EXTERNAL XLS FILE
data = xlsread('C:MATLABwork\CEE\CEEM1Qtrunc.xls','B2:H125');
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
data1 = data;
beta1 = beta;
Sig1 = Sig;
IRFmat1 = IRFmat;

%BOOTSTRAP CONFIDENCE BANDS
% save first P values in order to construct bootstrap data
datastart = data(1:P,:);
% find constants
cons = beta1(P*Q+1,:);
% initialize the bootstrap holding matrix, coeffeicients in columns,
% bootstrap iterations in rows
boot1 = zeros(N,Q*Q*(K+1));
for n=1:N
    %CONSTRUCT BOOTSTRAP DATA
    draw = unidrnd(S,S,1);  %choose random draws with replacement for errors
    uhat = u(draw,:);  %create uhat series
    %construct artificial data series
    data = zeros(T,Q); %initialize the series with zeros for all observations
    if Startmethod == 2
        draw2 = unidrnd(T,P,1);  %choose random draws with replacement for initial values
        datastart = data1(draw2,:);  % create initial values
    end
    data(1:P,:) = datastart(1:P,:); %read starting y's into first P observations
    % construct last T obervations using estimated model
    for t=P+1:T
        data(t,:) = uhat(t-P,:) + cons; %add in randomly drawn error term & constant
        for p=1:P
            data(t,:) = data(t,:) + data(t-p,:)*beta1(Q*(p-1)+1:Q*p,:); %add in coefficients times lags of data
        end
    end
    %GET y & x MATRICES
    CEE2getyx
    %ESTIMATE REDUCED FORM VAR
    [beta,Sig] = CEE2estim(y,x,Bootdf);
    %SOLVE FOR IRFS
    IRFmat = CEE2irfsolve(beta,Sig);
    % read these into holding matrix
    boot(n,:) = reshape(IRFmat,1,Q*Q*(K+1));  
end

% ANALYZE THE BOOTSTRAP RESULTS
bootavg = mean(boot);  %get average value of bootstrap A's
bootvar = var(boot); %get variances of bootstrap A's
% get avg IRFs from boot1avg
IRFavgmat = reshape(bootavg,Q,Q,K+1);
% DO PERCENTILES
% sort each column
bootsort = sort(boot);
% get upper bounds
IRFuppmat = reshape(bootsort(upper,:),Q,Q,K+1);
% get lower bounds
IRFlowmat = reshape(bootsort(lower,:),Q,Q,K+1);

% read selected IRFs from 3 dimensional arrays
CEE2plotirfs