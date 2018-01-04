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
    burnin   = 50000;  
    numsteps = 200000;
            
    totsteps = numsteps + burnin;
    
    cost = nan(totsteps, 2);
    Param = nan(nparam,totsteps);
    lp = length(pik);
    %sigmaTWTT = estimate_sigmaTWTT(seed);
    [sigma_ns,sigma_meters,SNR_dB] = radar_power(lp);
   
%% Initial proposal
    % propose parameter values
    prevParam(1:nparam-lp,1)=paramRange(1:nparam-lp,1)+rand(nparam-lp,1).*(paramRange(1:nparam-lp,2)-paramRange(1:nparam-lp,1));
    %set reasonable estimates of accumulation to start

%      prevParam(1,1) = 0.3981;
%      prevParam(2,1) = 0.1575; %deepest part
%      prevParam(3,1) = 0.1606;
%      prevParam(4,1) = 0.1013;
%      prevParam(5,1) = 0.1016;
%      prevParam(6,1) = 0.1386;
%      prevParam(7,1) = 0.1560;
%      prevParam(8,1) = 0.1477;
%      prevParam(9,1) = 0.1394; 
%      prevParam(10,1)= 0.1323; 
%     prevParam(11,1)= 0.1531; %shallowest part

%     prevParam(12,1)= 0.3119;
%     prevParam(13,1) = 168745989;
    
    regularization = nan;
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
    costTWTT = loglikelihoodTWTT(prevParam,lp,pik,sigma_ns);
    
    nAccept = 0;
    S = nan(totsteps,1);
    S(1) = nan;
%% Run Metropolis algorithm
    reg = 0;
    %for i = 1:15
    for i = 1:(totsteps)
    %for i = 502:1500
        if rem(i,10000) == 0 
             disp(i) 
        end
        %propose new parameters
        seed=seed+1;
        
        %if not at initial step, propose from previous accepted params
        if i~=1
            prevParam = Param(:,i-1);
        end
        propParam = proposeParams2(nparam,paramRange,prevParam,lp,H);
        %propParam(2:11,1) = smooth(propParam(2:11,1),3);
        
        % To compute error budget, set propParam to optimal value for each
        % param at a time. 
%            propParam(13,1) = 168881631.929264; %vice
        %propParam(14,1) = 27.5;
%           propParam(1,1) = 0.611;
%            propParam(12,1) = 0.390597; %h
        
%         propParam(2,1)  = 0.142275;
%         propParam(3,1)  = 0.118138;
%         propParam(4,1)  = 0.08495;
%         propParam(5,1)  = 0.0765647;
%         propParam(6,1)  = 0.118395;
%         propParam(7,1)  = 0.14147;
%         propParam(8,1)  = 0.1376075;
%         propParam(9,1)  = 0.11692;
%         propParam(10,1) = 0.121137;
%         propParam(11,1) = 0.14808;

        
        %calculate proposal cost
        [costAge_prop,S_prop,reg_prop]  = loglikelihoodAge(propParam,H,D,z,obsAge1950,accumFlag,lp);
        
%         S_prop = 0.03769;
        costTWTT_prop = loglikelihoodTWTT(propParam,lp,pik,sigma_ns);
        
        %accept/reject proposal
        random_number1 = rand;
        random_number2 = rand;
%         if i < 10            
%             i
%             S_prop, S
%             costAge_prop
%             costAge
%             costTWTT_prop
%             costTWTT
%             reg_prop
%             regularization
%             random_number1
%             random_number2
%             exp(-S_prop*(costAge_prop - costAge))
%             exp(-(costTWTT_prop - costTWTT))
%         end

        if (costTWTT_prop < costTWTT) || (exp(-(costTWTT_prop - costTWTT)) > random_number2) % accept TWTT
            if (costAge_prop < costAge) || (exp(-S_prop*(costAge_prop - costAge)) > random_number1) %accept Age params           
                % Accept all params if both likelihoods suggest to
                if i > burnin
                    nAccept=nAccept+1;
                end
%                 disp('Accept')
                
                Param(:,i) = propParam;
                S(i) = S_prop;
                costTWTT = costTWTT_prop;
                costAge = costAge_prop;
                regularization = reg_prop;
                %cost(i,:) = [costTWTT_prop costAge_prop];
            else %reject all
                %disp('Reject')
               Param(:,i) = prevParam;
               if i > 1
                    S(i) = S(i-1);
               end
            end
        else       %reject 
            %disp('Reject')
            Param(:,i) = prevParam;
            if i == 1
                regularization = 0;
            else
                regularization = reg(i-1);
            end
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