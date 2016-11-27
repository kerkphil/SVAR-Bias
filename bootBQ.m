% This program bootstraps a bivariate reduced form bivariate VAR of order P
% It also bootstraps the IRFs to obtain confidence bands.
% The identification scheme and data are the same as Blanchard & Quah

% This program calls the subroutine A0solveBQ.m
clear

global Sig B C A Acum K

N=200; %number of bootstrap iterations
P=8; % P is order of AR process to be estimated
K=40; %length of IRFs to display
confid = .05; % set confidence level
lower = floor(confid/2*N)+1;
upper = ceil((1-confid/2)*N);

% READ IN DATA
%x1=wk1read('C:MATLABwork\bootstrap\BQunadjusted.wk1');
x1=wk1read('C:MATLABwork\bootstrap\BQadjusted.wk1');
%x1=wk1read('C:MATLABwork\bootstrap\BQadjustedlong.wk1');

%x1 = x1(75:150,:);

[T,r1]=size(x1);
% %Plot data
% figure;
% subplot(2,1,1); plot(x1(:,1))
% title('growth of GDP')
% subplot(2,1,2); plot(x1(:,2))
% title('unemployment rate')

% ESTIMATE VAR
%save first P values in order to construct bootstrap data
ystart = x1(1:P,:);
%put current value into y vector
y = x1(P+1:T,:);
% delcare x vector: 2*P no constants, 2*P+1 includes constants
x = ones(T-P,2*P+1);
%put lagged values into x matrix
for p=1:P
    x(:,2*p-1:2*p) = x1(P+1-p:T-p,:);
end
T = T-P;  %recalculate sample size
%calculate AR coefficients
%regression 1 coeffs are in column 1, regression 2 are in column 2
XX = (x'*x)^(-1);
GAM = x*XX*x';
beta = XX*x'*y;  
eps = y - x*beta;  %calculate errors
%construct var-covar matrix
Sigest = eps'*eps./(T-2*P-1);
%Sigest = eps'*eps./T; %MLE
Sig = Sigest;

% FIND MA REPRESENTATION
% put elements of beta into approriate 2x2 matrices
B = zeros(2,2,P+1);
B(:,:,1) = [1 0; 0 1];
for p=2:P+1
    B(:,:,p)=-beta(2*p-3:2*p-2,:);
end
Best = B;
options = optimset('Display','off','TolX', 1e-8,'TolFun', 1e-10);
[A0est,fest,exitest] = fsolve(@A0solveBQ,Sigest,options);
% get IRFS A and Acum
IRFy1 = zeros(K,1);
IRFy2 = IRFy1;
IRFu1 = IRFy1;
IRFu2 = IRFy1;
for k=1:K
    IRFy1(k) = Acum(1,1,k);
    IRFu1(k) = A(1,2,k);
    IRFy2(k) = Acum(2,1,k);
    IRFu2(k) = A(2,2,k);
end
Aorig = A;
A0orig=A0est;
IRFy1orig = IRFy1;
IRFu1orig = IRFu1;
IRFy2orig = IRFy2;
IRFu2orig = IRFu2;

% PERFORM BOOTSTRAP
% for matrices below, bootstrap iterations are in rows
boot1 = zeros(N,2*2*K);  %initialize matrix to store right bootstrap A coeffs
boot1w = zeros(N,2*2*K);  %initialize matrix to store wrong bootstrap A coeffs
% run bootstraps
for n=1:N
    draw = unidrnd(T,T,1);  %chose random draws with replacement for errors
    epshat = eps(draw,:);  %create uhat series
    %construct artificial data series
    x1new = zeros(T+P,2); %initialize the series with zeros for all observations
    draw2 = unidrnd(T+P,P,1);  %chose random draws with replacement for initial values
    ystartnew = x1(draw2,:);  % creat initial values
    for t=1:P
        x1new(t,:) = ystartnew(t,:); %read randomly drawn y's into first P observations
    end
    % construct last T obervations using estimated model
    for t=P+1:T+P
        x1new(t,:) = epshat(t-P,:); %add in randomly drawn error term
        for p=1:P
            x1new(t,:) = x1new(t,:)+x1new(t-p,:)*beta(2*p-1:2*p,:); %add in coefficients times lags of x1
        end
    end
    %put current value into y vector
    ynew = x1new(P+1:T,:);
    % delcare x vector: 2*P no constants, 2*P+1 includes constants
    xnew = ones(T-P,2*P+1);
    %put lagged values into x matrix
    for p=1:P
        xnew(:,2*p-1:2*p) = x1new(P+1-p:T-p,:);
    end
    % estimate using OLS
    XXnew = (xnew'*xnew)^(-1);
    betanew = XXnew*xnew'*ynew; %calculate and store the estimated betas
    epsnew = ynew - xnew*betanew;  %calculate errors
    %construct var-covar matrix
    Sigright = epsnew'*epsnew./((T-2*P-1)^2/T); %right way
    Sigwrong = epsnew'*epsnew/(T-2*P-1); %wrong way
    %SigMLE = epsnew'*epsnew/T; %wrong way - MLE
    
    Sig = Sigright;
    % put elements of beta into approriate 2x2 matrices
    B = zeros(2,2,P+1);
    B(:,:,1) = [1 0; 0 1];
    for p=2:P+1
        B(:,:,p)=-betanew(2*p-3:2*p-2,:);
    end
    % solve for A coffeicients by getting A0 from A0solve
    [A0new,fnew,exitnew] = fsolve(@A0solveBQ,A0est,options);
    for k=1:K
        boot1(n,4*k-3:4*k) = [Acum(1,1,k) A(1,2,k) Acum(2,1,k) A(2,2,k)];
    end
    
    Sig = Sigwrong;
    % put elements of beta into approriate 2x2 matrices
    B = zeros(2,2,P+1);
    B(:,:,1) = [1 0; 0 1];
    for p=2:P+1
        B(:,:,p)=-betanew(2*p-3:2*p-2,:);
    end
    % solve for A coffeicients by getting A0 from A0solve
    [A0new,fnew,exitnew] = fsolve(@A0solveBQ,A0est,options);
    for k=1:K
        boot1w(n,4*k-3:4*k) = [Acum(1,1,k) A(1,2,k) Acum(2,1,k) A(2,2,k)];
    end

end

% ANALYZE THE RIGHT BOOTSTRAP RESULTS
boot1avg = mean(boot1);  %get average value of bootstrap A's
boot1var = var(boot1); %get variances of bootstrap A's

% get avg IRFs from boot1avg
IRFy1a = zeros(K,1);
IRFy2a = zeros(K,1);
IRFu1a = zeros(K,1);
IRFu2a = zeros(K,1);
for k=1:K
    IRFy1a(k) = boot1avg(1,4*k-3);
    IRFu1a(k) = boot1avg(1,4*k-2);
    IRFy2a(k) = boot1avg(1,4*k-1);
    IRFu2a(k) = boot1avg(1,4*k);
end

% DO PERCENTILES
% sort each column
bootsort = sort(boot1);
% get lower bounds
IRFy1l = zeros(K,1);
IRFy2l = zeros(K,1);
IRFu1l = zeros(K,1);
IRFu2l = zeros(K,1);
for k=1:K
    IRFy1l(k) = bootsort(lower,4*k-3);
    IRFu1l(k) = bootsort(lower,4*k-2);
    IRFy2l(k) = bootsort(lower,4*k-1);
    IRFu2l(k) = bootsort(lower,4*k);
end
% get upper bounds
IRFy1u = zeros(K,1);
IRFy2u = zeros(K,1);
IRFu1u = zeros(K,1);
IRFu2u = zeros(K,1);
for k=1:K
    IRFy1u(k) = bootsort(upper,4*k-3);
    IRFu1u(k) = bootsort(upper,4*k-2);
    IRFy2u(k) = bootsort(upper,4*k-1);
    IRFu2u(k) = bootsort(upper,4*k);
end

%Blanchard Quah bounds
IRFy1bqu = zeros(K,1);
IRFy2bqu = zeros(K,1);
IRFu1bqu = zeros(K,1);
IRFu2bqu = zeros(K,1);
IRFy1bql = zeros(K,1);
IRFy2bql = zeros(K,1);
IRFu1bql = zeros(K,1);
IRFu2bql = zeros(K,1);
Ny1u = zeros(K,1);
Nu1u = zeros(K,1);
Ny2u = zeros(K,1);
Nu2u = zeros(K,1);
for k=1:K
    for n=1:N
        if bootsort(n,4*k-3)>boot1avg(1,4*k-3)
            IRFy1bqu(k) = IRFy1bqu(k) + (bootsort(n,4*k-3)-IRFy1a(k))^2;
            Ny1u(k) = Ny1u(k) +1;
        else
            IRFy1bql(k) = IRFy1bql(k) + (bootsort(n,4*k-3)-IRFy1a(k))^2;
        end
        if bootsort(n,4*k-2)>boot1avg(1,4*k-2)
            IRFu1bqu(k) = IRFu1bqu(k) + (bootsort(n,4*k-2)-IRFu1a(k))^2;
            Nu1u(k) = Nu1u(k) +1;
        else
            IRFu1bql(k) = IRFu1bql(k) + (bootsort(n,4*k-2)-IRFu1a(k))^2;
        end
        if bootsort(n,4*k-1)>boot1avg(1,4*k-1)
            IRFy2bqu(k) = IRFy2bqu(k) + (bootsort(n,4*k-1)-IRFy2a(k))^2;
            Ny2u(k) = Ny2u(k) +1;
        else
            IRFy2bql(k) = IRFy2bql(k) + (bootsort(n,4*k-1)-IRFy2a(k))^2;
        end
        if bootsort(n,4*k)>boot1avg(1,4*k)
            IRFu2bqu(k) = IRFu2bqu(k) + (bootsort(n,4*k)-IRFu2a(k))^2;
            Nu2u(k) = Nu2u(k) +1;
        else
            IRFu2bql(k) = IRFu2bql(k) + (bootsort(n,4*k)-IRFu2a(k))^2;
        end
    end
end
Ny1l = N*ones(K,1) - Ny1u;
Nu1l = N*ones(K,1) - Nu1u;
Ny2l = N*ones(K,1) - Ny2u;
Nu2l = N*ones(K,1) - Nu2u;
IRFy1bqu = IRFy1a+(IRFy1bqu./Ny1u).^.5;
IRFu1bqu = IRFu1a+(IRFu1bqu./Nu1u).^.5;
IRFy2bqu = IRFy2a+(IRFy2bqu./Ny2u).^.5;
IRFu2bqu = IRFu2a+(IRFu2bqu./Nu2u).^.5;
IRFy1bql = IRFy1a-(IRFy1bql./Ny1l).^.5;
IRFu1bql = IRFu1a-(IRFu1bql./Nu1l).^.5;
IRFy2bql = IRFy2a-(IRFy2bql./Ny2l).^.5;
IRFu2bql = IRFu2a-(IRFu2bql./Nu2l).^.5;

% ANALYZE THE WRONG BOOTSTRAP RESULTS
boot1wavg = mean(boot1w);  %get average value of bootstrap A's
boot1wvar = var(boot1w); %get variances of bootstrap A's

% get avg IRFs from boot1wavg
IRFy1aw = zeros(K,1);
IRFy2aw = IRFy1aw;
IRFu1aw = IRFy1aw;
IRFu2aw = IRFy1aw;
for k=1:K
    IRFy1aw(k) = boot1wavg(1,4*k-3);
    IRFu1aw(k) = boot1wavg(1,4*k-2);
    IRFy2aw(k) = boot1wavg(1,4*k-1);
    IRFu2aw(k) = boot1wavg(1,4*k);
end

% DO PERCENTILES
% sort each column
bootsort = sort(boot1w);
% get lower bounds
IRFy1lw = zeros(K,1);
IRFy2lw = IRFy1lw;
IRFu1lw = IRFy1lw;
IRFu2lw = IRFy1lw;
for k=1:K
    IRFy1lw(k) = bootsort(lower,4*k-3);
    IRFu1lw(k) = bootsort(lower,4*k-2);
    IRFy2lw(k) = bootsort(lower,4*k-1);
    IRFu2lw(k) = bootsort(lower,4*k);
end
% get upper bounds
IRFy1uw = zeros(K,1);
IRFy2uw = IRFy1uw;
IRFu1uw = IRFy1uw;
IRFu2uw = IRFy1uw;
for k=1:K
    IRFy1uw(k) = bootsort(upper,4*k-3);
    IRFu1uw(k) = bootsort(upper,4*k-2);
    IRFy2uw(k) = bootsort(upper,4*k-1);
    IRFu2uw(k) = bootsort(upper,4*k);
end

%Blanchard Quah bounds
IRFy1bquw = zeros(K,1);
IRFy2bquw = zeros(K,1);
IRFu1bquw = zeros(K,1);
IRFu2bquw = zeros(K,1);
IRFy1bqlw = zeros(K,1);
IRFy2bqlw = zeros(K,1);
IRFu1bqlw = zeros(K,1);
IRFu2bqlw = zeros(K,1);
Ny1u = zeros(K,1);
Nu1u = zeros(K,1);
Ny2u = zeros(K,1);
Nu2u = zeros(K,1);
for k=1:K
    for n=1:N
        if bootsort(n,4*k-3)>boot1wavg(1,4*k-3)
            IRFy1bquw(k) = IRFy1bquw(k) + (bootsort(n,4*k-3)-IRFy1aw(k))^2;
            Ny1u(k) = Ny1u(k) +1;
        else
            IRFy1bqlw(k) = IRFy1bqlw(k) + (bootsort(n,4*k-3)-IRFy1aw(k))^2;
        end
        if bootsort(n,4*k-2)>boot1wavg(1,4*k-2)
            IRFu1bquw(k) = IRFu1bquw(k) + (bootsort(n,4*k-2)-IRFu1aw(k))^2;
            Nu1u(k) = Nu1u(k) +1;
        else
            IRFu1bqlw(k) = IRFu1bqlw(k) + (bootsort(n,4*k-2)-IRFu1aw(k))^2;
        end
        if bootsort(n,4*k-1)>boot1wavg(1,4*k-1)
            IRFy2bquw(k) = IRFy2bquw(k) + (bootsort(n,4*k-1)-IRFy2aw(k))^2;
            Ny2u(k) = Ny2u(k) +1;
        else
            IRFy2bqlw(k) = IRFy2bqlw(k) + (bootsort(n,4*k-1)-IRFy2aw(k))^2;
        end
        if bootsort(n,4*k)>boot1wavg(1,4*k)
            IRFu2bquw(k) = IRFu2bquw(k) + (bootsort(n,4*k)-IRFu2aw(k))^2;
            Nu2u(k) = Nu2u(k) +1;
        else
            IRFu2bqlw(k) = IRFu2bqlw(k) + (bootsort(n,4*k)-IRFu2aw(k))^2;
        end
    end
end
Ny1l = N*ones(K,1) - Ny1u;
Nu1l = N*ones(K,1) - Nu1u;
Ny2l = N*ones(K,1) - Ny2u;
Nu2l = N*ones(K,1) - Nu2u;
IRFy1bquw = IRFy1aw+(IRFy1bquw./Ny1u).^.5;
IRFu1bquw = IRFu1aw+(IRFu1bquw./Nu1u).^.5;
IRFy2bquw = IRFy2aw+(IRFy2bquw./Ny2u).^.5;
IRFu2bquw = IRFu2aw+(IRFu2bquw./Nu2u).^.5;
IRFy1bqlw = IRFy1aw-(IRFy1bqlw./Ny1l).^.5;
IRFu1bqlw = IRFu1aw-(IRFu1bqlw./Nu1l).^.5;
IRFy2bqlw = IRFy2aw-(IRFy2bqlw./Ny2l).^.5;
IRFu2bqlw = IRFu2aw-(IRFu2bqlw./Nu2l).^.5;


        
% plot IRFs
y1=[IRFy1 IRFy1l IRFy1a IRFy1u IRFy1lw IRFy1aw IRFy1uw];
y2=[IRFu1 IRFu1l IRFu1a IRFu1u IRFu1lw IRFu1aw IRFu1uw];
y3=[IRFy2 IRFy2l IRFy2a IRFy2u IRFy2lw IRFy2aw IRFy2uw];
y4=[IRFu2 IRFu2l IRFu2a IRFu2u IRFu2lw IRFu2aw IRFu2uw];
%compare1(y1,y2,y3,y4)

y1=[IRFy1 IRFy1bql IRFy1a IRFy1bqu IRFy1bqlw IRFy1aw IRFy1bquw];
y2=[IRFu1 IRFu1bql IRFu1a IRFu1bqu IRFu1bqlw IRFu1aw IRFu1bquw];
y3=[IRFy2 IRFy2bql IRFy2a IRFy2bqu IRFy2bqlw IRFy2aw IRFy2bquw];
y4=[IRFu2 IRFu2bql IRFu2a IRFu2bqu IRFu2bqlw IRFu2aw IRFu2bquw];
%compare1(y1,y2,y3,y4)

xlswrite('bootBQirf.xls',y1,'y1');
xlswrite('bootBQirf.xls',y2,'u1');
xlswrite('bootBQirf.xls',y3,'y2');
xlswrite('bootBQirf.xls',y4,'u2');
xlswrite('bootBQirf.xls',beta,'beta');
xlswrite('bootBQirf.xls',A0est,'A0');