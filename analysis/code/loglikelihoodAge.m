function [cost,S,regularization] = loglikelihoodAge(param,H,D,z,obsAge1950,accumFlag,lp)
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
    age= byrdModels2(param,z,H,accumFlag,lp);
    
%% Reverse reference frame of observed depths to be elevation above bed
    ZZ=H-D; 
    
%% Index to find depths where both obs simulated ages exist for comparison
    % Model evaluated at 1m intervals, so round volcanic z to compare
    [~,ind]=intersect(z,round(ZZ)); 
    ind=flipud(ind);
    
%% Variance of observed age (to be modulated by S)
    sigAge = 0.01*(obsAge1950);
    varAge=sigAge.^2;
    
%% Compute S
    % from Charles's GMDD paper
    % S = Ga(ke/2+alpha, E(m)+beta)
    
    %Regularization
    % Calculate the variance with a moving average for the mean
    % var = 1/(N-1) sum(|Ai - mu)|^2)
    N = length(param)-lp-3;
    var = 1/(N-1)*sum((param(2:length(param)-lp-3) - smooth(param(2:length(param)-lp-3),3)).^2);
    Reg = var;
    %Reg = var(param(2:length(param)-lp-3,1));
    %accRefReg = 0.000204421610432737;
    accRefReg = 1.5989e-04;
     regularization = Reg/accRefReg;

     E_m = sum((age(ind)-obsAge1950).^2./varAge);
     ke = lp; %number of horizons (assumes independence)
     alpha = 1;
     beta = 1;

     S = gamrnd(ke/2+alpha, 1./(E_m+beta));
     
%% Evaluate loglikelihood
     cost = E_m/2 + regularization/2;
     
end

