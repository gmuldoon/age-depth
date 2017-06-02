function [paramDistrib,burnin,numsteps,ageDistrib,obsAge1950,D,H,z,pikDepthSamp,pikDepthStats,pikAge,pikAgeStats,posterior]=calcAgeDepth(accumFlag,datFlag,core,plotting)
% DESCRIPTION:
% This function performs the steps of the Byrd chronology analysis:
% 1. Set up the environment
% 2. Load data including observed chronologies and layer TWTTs.
% 3. Run Metropolis algorithm to derive flow parameters
% 4. Derive depth from layer TWTT
% 5. Sample layer depth uncertainty
% 6. Calculated layer age given layer depth
% 7. Calculate convenient stats on the result
%
% NOTE:
% plotAgeDepth.m calls this function and plots the results. That script
% makes it possible to loop through this analysis for different
% configurations
%
% INPUT:
% - *accumFlag* specifies the accumulation function to use, as defined in setParams.m
% - *datFlag* specifies the observed chronology at Byrd to use 
% - *core* is ice core to derive chronology for (Byrd or WD)
%
%
% OUTPUT:
% - *paramDistrib* is a matrix with ensemble of values for each parameter (nparam x totsteps)
% (*totsteps* in numsteps + burnin)
% - *burnin* is number of iterations to ignore from the start of the Metropolis sampling
% - *numsteps* is number of Metropolis samples after burnin
% - *ageDistrib* is an H by totsteps matrix of ages for each meter depth in the core
% - *obsAge1950* is observed chronology (years before 1950)
% - *D* is observed chronology depth (m)
% - *H* is length of the core (m)
% - *z* is height above the bed (m)
% - *pikDepthSamp* is is the lp-by-numsteps ensemble of depths for each horizon
% - *pikDepthStats* contains the TWTT, mean/median/stddev of horizon depth
% - *pikAge* is the lp-by-numsteps ensemble of ages for each horizon
% - *pikAgeStats* contains the TWTT, mean/median/stddev of horizon age
% - *posterior* contains the posterior value for each metropolis iteration (1 x totsteps)
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

    tic
%% 1. Set up environment
    addpath('../Matlab_functions');
    set(0,'DefaultFigureWindowStyle','docked')
    seed=123;      % seed for reproducibility

%% 2. Load data
    [H,D,z,dFirn,pik,obsAge1950,obsAge1950Unc]=loadCorePikData(core,datFlag);

%% 3. Calculate parameter sets which give reasonable age profiles
    [paramDistrib,ageDistrib,numsteps,posterior,burnin]=metropolisAgeSampler(core,D,z,H,obsAge1950,accumFlag,seed);

%% 4. Transform TWTT to depth with uncertainty from radar
    [pikDepth,lp,nsamp]=twtt2depth(pik,H,core,plotting,dFirn,z);

%% 5. Sample depth of each horizon to within uncertainty
    pikDepthSamp = sampleDepth(numsteps,nsamp,pikDepth,lp);

%% 6. Assign age to sampled depths
    pikAge = calcPikAge(pikDepthSamp,ageDistrib,z,core,H,obsAge1950,obsAge1950Unc);

%% 7. Calculate mean/median/std of age and depth for each radar horizon
    [pikDepthStats,pikAgeStats,fullAgeStats] = getAgeDepthStats(core,pik,pikDepthSamp,pikAge,ageDistrib,H,z,obsAge1950,obsAge1950Unc);

%% End
    toc
end
