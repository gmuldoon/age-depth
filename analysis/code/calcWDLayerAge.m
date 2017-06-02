function [WDpikAge,WDpikDepthSamp,WDpikDepthStats,WDpikAgeStats]=calcWDLayerAge(numsteps,plotting)
% DESCRIPTION:
% This function calculates the ages of horizons at WD
%
% INPUT:
% - *plotting* is a flag that makes plots related to the depth conversion if true
% 
% OUTPUT:
% - *WDpikDepth* is an lp-by-nsamp matrix of depths of each horizon with unc
% - *WDpikAge* is an lp-by-nsamp matrix of ages of each horizon with unc
% - *WDpikDepthStats* contains TWTT, mean/median/stddev of each WD horizon depth
% - *WDpikAgeStats* contains TWTT, mean/median/stddev of each WD horizon age
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017
%
%% Setup
datFlag = nan;  %don't need this for WD
core = 'WD';

%% Load data for wd
[H,~,z,dFirn,pik,obsAge1950,obsAge1950Unc]=loadCorePikData(core,datFlag);

%% Convert TWTT to depth for piks
[WDpikDepth,lp,nsamp]=twtt2depth(pik,H,core,plotting,dFirn,z);

%% Sample depth uncertainty
WDpikDepthSamp = sampleDepth(numsteps,nsamp,WDpikDepth,lp);

%% Assign age to sampled depths
ageDistrib=nan; %because there is none for WD
WDpikAge = calcPikAge(WDpikDepthSamp,ageDistrib,z,core,H,obsAge1950,obsAge1950Unc);
%WDpikAge=flipud(WDpikAge);

%% get stats on age and depth for easy plotting
ageDistrib=nan; %don't need this for wd
[WDpikDepthStats,WDpikAgeStats,~] = getAgeDepthStats(core,pik,WDpikDepthSamp,WDpikAge,ageDistrib,H,z,obsAge1950,obsAge1950Unc);

