function loglike = loglikelihood(param,H,D,z,obsAge1950,accumFlag)
% DESCRIPTION:
% This function computes the log likelihood value of the age profile
% derived using proposed parameter values compared to observations
% 
% INPUT:
% - *param* is array of proposed parameters (size nparam)
% - *H* is length of the core (m)
% - *D* is observed chronology depth (m)
% - *z* is height above the bed (m)
% - *obsAge1950* is observed chronology (years before 1950)
% - *accumflag* specifies the accumulation function to use, as defined in setParams.m
% 
% OUTPUT:
% - *loglike* value of the resulting loglikelihood
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017
%% Simulate age profile from proposed parameters
    age= byrdModels(param,z,H,accumFlag);
    
%% Reverse reference frame of observed depths to be elevation above bed
    ZZ=H-D; 
    
%% Index to find depths where both obs simulated ages exist for comparison
    [~,ind]=intersect(z,round(ZZ)); % Model evaluated at 1m intervals, so round volcanic z to compare
    ind=flipud(ind);
    
%% Find variance of observed age
    if length(param) == 7
        %sigAge = param(7)*obsAge1950;
        sigAge=param(7);
    else
        sigAge = 0.03*(obsAge1950);
    end
    varAge=sigAge.^2;
    
%% Compute loglikelihood        
     loglike = -(sum((age(ind)-obsAge1950).^2.*(varAge)/2));%/length(obsAge));
     %loglike = -(sum((age(ind)-obsAge).^2)/(2*varAge))'/length(obsAge);
end

