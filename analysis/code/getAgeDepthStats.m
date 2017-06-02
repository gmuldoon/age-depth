function [pikDepthStd,pikAgeStd,fullAgeStd] = getAgeDepthStats(core,pik,pikDepth,pikAge,ageDistrib,H,z,obsAge1950,obsAge1950Unc)
% DESCRIPTION:
% This file computes basic statitics (mean, median, stddev) of the ensemble
% of age and depth for each horizon and for the whole ice column. These are
% useful for evaluation and plotting.
%
% INPUT:
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *pik*
% - *pikDepth* is 
% - *pikAge* is the lp-by-numsteps ensemble of ages for each horizon
% - *ageDistrib* is an H by totsteps matrix of ages for each meter depth in the core
% - *H* is length of the core (m)
% - *z* is height above the bed (m)
% - *obsAge1950* is observed chronology (years before 1950)
% - *obsAge1950Unc* is the age uncertaint for the WD core (or nan for Byrd)
%
% OUTPUT:
% - *pikDepthStd* contains TWTT, mean/median/stddev of each horizon depth
% - *pikAgeStd* contains TWTT, mean/median/stddev of each horizon age
% - *fullAgeStd* contains z, mean/median/stddev of ice column age

% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 16 Feb 2017

disp('Computing ensemble stats')
%% Initialize some things
    lp = length(pikDepth(:,1));
    pikDepthStd=nan(lp,4);
    pikAgeStd = nan(lp,4);
    
%% Compute stats for horizons only
    for i = 1:lp
        %stats on horizon depth
        pikDepthStd(i,1) = pik(i);                  % TWTT               
        pikDepthStd(i,2) = mean(pikDepth(i,:));     % mean depth         
        pikDepthStd(i,3) = median(pikDepth(i,:));   % median depth
        pikDepthStd(i,4) = std(pikDepth(i,:));      % stddev in depth

        %stats on horizon age
        pikAgeStd(i,1) = pik(i);                    % TWTT
        pikAgeStd(i,2) = mean(pikAge(i,:));         % mean age
        pikAgeStd(i,3) = median(pikAge(i,:));       % median age
        pikAgeStd(i,4) = std(pikAge(i,:));          % stddev in age
    end
    format long g
%% Compute stats for full profile (every 1m) in ice column
    fullAgeStd = nan(H,4);
    if strcmp(core,'WD')
        fullAgeStd(:,1) = z;            % z profile
        fullAgeStd(:,2) = obsAge1950;   % mean age
        fullAgeStd(:,3) = obsAge1950;   % median age
        fullAgeStd(:,4) = obsAge1950Unc;% std in age
    elseif strcmp(core,'Byrd')
        fullAgeStd(:,1) = z;                    % z profile
        fullAgeStd(:,2) = mean(ageDistrib,2);   % mean age
        fullAgeStd(:,3) = median(ageDistrib,2); % median age
        fullAgeStd(:,4) = std(ageDistrib,0,2);  % 0 is normal weighting (normalize by N-1)
    end

end

