function cost = loglikelihoodTWTT(Param,lp,TWTTpik,sigma_ns)
% DESCRIPTION:
% This function computes the log likelihood value of the age profile
% derived using proposed parameter values compared to observations
% 
% INPUT:
% 
% OUTPUT:
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 29 March 2017

%% get piks in same units
%TWTTpik = pik*10^6;
 nparam=length(Param(:,1));

%% 
    vIce = Param(nparam-lp-1,1);
%     dFirn = Param(nparam-lp,1);
    Dr = Param(length(Param)-lp+1:length(Param),1);

%% Simulate TWTT of Dr
    TWTTm = Dr_TWTT(Dr,vIce,sigma_ns/100);
    %TWTTpik'
    
%% Find variance of depth error
    varTWTT = (sigma_ns(2:end)/100).^2; %skip the surface sigma, i = 1
    
%% Compute loglikelihood        
    cost = sum((TWTTm-TWTTpik').^2./(2*varTWTT)');

end