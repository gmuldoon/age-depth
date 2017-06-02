function [Param,numsteps,cost,burnin,S,reg]=metropolisAgeSampler2(D,z,H,obsAge1950,accumFlag,seed,paramRange,nparam,pik)
% DESCRIPTION:
% This functions runs a metropolis sampler to compute an ensemble of ice
% flow parameters which agree to within uncertainty with an observed ice 
% core chronology 
% 
% INPUT:
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
    numsteps = 10000;
            
    totsteps = numsteps + burnin;
    
    cost = nan(totsteps, 2);
    Param = nan(nparam,totsteps);
    lp = length(pik);
    sigmaTWTT = estimate_sigmaTWTT(seed);
   
%% Initial proposal
    % propose parameter values
    prevParam(1:nparam-lp,1)=paramRange(1:nparam-lp,1)+rand(nparam-lp,1).*(paramRange(1:nparam-lp,2)-paramRange(1:nparam-lp,1));
    %prevParam(2:5,1) = normrnd(paramRange(2:5,1),paramRange(2:5,2));
    %Param(1:nparam-lp,1) = proposeParams2(nparam,paramRange,param,lp);
    % 
    %prevParam(nparam-lp+1:nparam,1) = [500,600,900,1300,1700]';
    %prevParam(nparam-lp+1:nparam,1) = [400 800 1000 1400 1800];
    v_est = 1.68*10^8;
    pik_sec = pik*10^-6;
    depth_est = v_est*pik_sec/2;
    prevParam(nparam-lp+1:nparam,1) = depth_est;         
                
    %calculate an initial proposal
    [costAge,~,~] = loglikelihoodAge(prevParam,H,D,z,obsAge1950,accumFlag,lp);
    costTWTT = loglikelihoodTWTT(prevParam,lp,pik,sigmaTWTT);
    
    nAccept = 0;
    S = nan(totsteps,1);

%% Run Metropolis algorithm
    for i = 1:(totsteps)
    %for i = 502:1500
        if rem(i,5000) == 0 
             disp(i) 
        end
        %propose new parameters
        seed=seed+1;
        
        %if not at initial step, propose from previous accepted params
        if i~=1
            prevParam = Param(:,i-1);
        end
        propParam = proposeParams2(nparam,paramRange,prevParam,lp,H);
        
        %calculate proposal cost
        [costAge_prop,S_prop,reg_prop]  = loglikelihoodAge(propParam,H,D,z,obsAge1950,accumFlag,lp);
        costTWTT_prop = loglikelihoodTWTT(propParam,lp,pik,sigmaTWTT);
        
        %accept/reject proposal    
        
        if (costAge_prop < costAge) || (exp(-S_prop*(costAge_prop - costAge)) > rand) %accept Age params
            if (costTWTT_prop < costTWTT) || (exp(-(costTWTT_prop - costTWTT)) > rand) % accept TWTT
                % Accept all params if both likelihoods suggest to
                if i > burnin
                    nAccept=nAccept+1;
                end
                Param(:,i) = propParam;
                S(i) = S_prop;
                costTWTT = costTWTT_prop;
                costAge = costAge_prop;
                regularization = reg_prop;
                %cost(i,:) = [costTWTT_prop costAge_prop];
            else %reject all
               Param(:,i) = prevParam; 
               if i > 1
                    S(i) = S(i-1);
               end
            end
        else       %reject 
            Param(:,i) = prevParam;
            if i > 1
                S(i) = S(i-1);
           end
        end
        cost(i,1) = costTWTT;
        cost(i,2) = costAge;
        reg(i,1) = regularization;
    end   
    fprintf('Metropolis acceptance rate: %0.4f%%\n',nAccept/numsteps*100.)
end