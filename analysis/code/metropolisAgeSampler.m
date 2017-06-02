function [Param,Age,numsteps,posterior,burnin]=metropolisAgeSampler(core,D,z,H,obsAge1950,accumFlag,seed)
% DESCRIPTION:
% This functions runs a metropolis sampler to compute an ensemble of ice
% flow parameters which agree to within uncertainty with an observed ice 
% core chronology 
% 
% INPUT:
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *D* is observed chronology depth (m)
% - *z* is height above the bed (m)
% - *H* is length of the core (m)
% - *obsAge1950* is observed chronology (years before 1950)
% - *accumflag* specifies the accumulation function to use, as defined in setParams.m
% - *seed* is an arbitrary seed for random number generators (for
%   reproducibility across runs)
%
% OUTPUT:
% - *nAccept* is the number of Metropolis iterations with new accepted params
% - *Param* is ensemble of values for each parameter (nparam x totsteps)
% - *Age* is an H by totsteps matrix of ages for each meter depth in the core
% - *numsteps* is number of Metropolis samples after burnin
% - *posterior* contains the posterior value for each metropolis iteration (1 x totsteps)
% - *burnin* is number of iterations to ignore from the start of the Metropolis sampling
% 
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

disp('Running the metropolis algorithm')

%% Set algorithm parameters
    rng(seed);
    burnin   = 5000;  
    numsteps = 50000;

    %this one requires a LOT of iterations
    if accumFlag == 100
          burnin=50000;
          numsteps = 75000;
    end
            
    totsteps = numsteps + burnin;
    
    posterior = nan(totsteps, 1);
    Age = nan(H,totsteps);
   
%% Initialize posterior 
    [paramRange,nparam] = setParams(accumFlag,H);
    param=paramRange(:,1)+rand(nparam,1).*(paramRange(:,2)-paramRange(:,1));
    Param(:,1) = proposeParams(nparam,paramRange,param);
    
    loglike = loglikelihood(Param(:,1),H,D,z,obsAge1950,accumFlag);
    
    prevStep = loglike;
    nAccept = 0;

%% Run Metropolis algorithm
    for i = 1:(totsteps)
        %propose new parameters
        seed=seed+1;
        propParam = proposeParams(nparam,paramRange,Param(:,i));
        
        %calculate proposed posterior
        loglike=loglikelihood(propParam,H,D,z,obsAge1950,accumFlag);
        prop = loglike;
        
        %accept/reject proposal
        rho = exp(prop-prevStep);
        u = rand;
        if rho > u %accept
            nAccept=nAccept+1;
            if i < totsteps
                Param(:,i+1) = propParam; 
            end
            prevStep = prop;
            posterior(i) = exp(prop);
        else       %reject 
            if i < totsteps
                Param(:,i+1) = Param(:,i);
            end
            posterior(i) = exp(prevStep);
        end
        
        % compute age profile from this iteration's parameters
        Age(:,i) = getAge(core,Param(:,i),z,H,accumFlag);
    end   
    fprintf('Metropolis acceptance rate: %0.4f%%\n',nAccept/numsteps*100.)
end