function logpri = logprior(propParam,paramRange,nparam)
% DESCRIPTION:
% This function computes the log prior of the proposed parameters
% 
% INPUT:
% - *propParam* is the new set of parameter values (size nparam)
% - *paramRange* defines the min and max allowed parameter values
% - *nparam* is the number of parameters to use
% 
% OUTPUT:
% - *logpri* value of the resulting logprior
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017
%% Compute log prior (sum of log of each prior distrib)
    logpri=0;
    if nparam ~= 7
        for i = 1:nparam
            logpri = logpri + log(unifpdf(propParam(i),paramRange(i,1),paramRange(i,2)));
        end
    elseif nparam == 7 % Need to treat the obs age uncertainty separately
        for i = 1:nparam-1
            logpri = logpri + log(unifpdf(propParam(i),paramRange(i,1),paramRange(i,2)));
        end
        logpri = logpri + log(normpdf(propParam(7),paramRange(7,1),paramRange(7,2)));
    end
        
end

