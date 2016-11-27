%This set of commands checks if the original IRF lies withing the 
% bootstrapped confidence bands.

Coverage = zeros(K+1,1);  %lies within the bands
Above = zeros(K+1,1); %lies above the upper band
Below = zeros(K+1,1); %lies below the lower band
for k=1:K+1
    if (IRForig(k)<=IRFupp(k)) && (IRForig(k)>=IRFlow(k))
        Coverage(k,1)=1;
    end 
    if (IRForig(k)>IRFupp(k))
        Above(k,1)=1;
    end
    if (IRForig(k)<IRFlow(k))
        Below(k,1)=1;
    end
end
