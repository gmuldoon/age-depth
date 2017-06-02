function [H,D,z,pik,obsAge1950]=loadCorePikData2(core,datFlag)
% DESCRIPTION:
% This function loads data necessary for the analysis in calcAgeDepth.m
%
% NOTE:
% It loads:
% - Byrd and WAIS Divide ice core characteristics
% - Byrd and WAIS Divide ice core observed chronologies
% - Radar horizon TWTT at each ice core
% - Uncertainty in chronology for WD
%
% INPUT:
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *datFlag* specifies the observed chronology at Byrd to use 
%  
% OUTPUT:
% - *H* is length of the core (m)
% - *D* is observed chronology depth (m)
% - *z* is height above the bed (m)
% - *dFirn* is depth of the firn layer (m)
% - *obsAge1950* is observed chronology (years before 1950)
% - *obsAge1950Unc* is uncertainty on observed chronology
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

    fprintf('Loading age depth data at the %s ice core\n',core)
%% Load/define Byrd ice core data and info
    if strcmp(core,'Byrd')
        H = 2164;    % ice thickness at Byrd [m]
        dFirn=64;    % depth of the firn layer

        %TWTT for traced radar horizons
        filename ='../data/byrd_TWTT_5.txt';
        ncol=1;
        pik=readFloats(filename,ncol)';

        if strcmp(datFlag,'volcanic')   % Blunier et al. 2001 volcanic Byrd chronology back to 50 ka
            disp('Using volcanic data')

            filename='../data/byrd_volcage.txt';
            ncol=2;
            byrdObs=readFloats(filename,ncol);
            D1968=byrdObs(:,1);
            obsAge = byrdObs(:,2);

        elseif strcmp(datFlag,'oxygen') % From methane-based synchronization with Greenland (relative ages)
            disp('Using oxygen data')

            filename='../data/ByrdCH4.csv';
            ncol=2;
            byrdObs=readFloats(filename,ncol);
            D1968=byrdObs(:,2);
            obsAge = byrdObs(:,1);
        end

        % add accum between drilling and obs/WD drilling
        accum = 0.11;
        D = D1968+(2013-1968)*accum; 

        % set ages relative to 1950 like WD core chronology
        obsAge1950=obsAge-(1968-1950); 

        %Only use those depths which are different when rounded
        [D,ind] = unique(ceil(D));
        obsAge1950=obsAge1950(ind);
        obsAge1950Unc = nan(size(obsAge1950)); %We have unc for WD, not for Byrd
    end
%% Load/define WAIS Divide (WD) ice core data and info
    if strcmp(core,'WD')
        dFirn=159;      % depth of firn layer

        %TWTT for traced radar horizons
        filename='../data/WD_TWTT_5.txt';
        ncol=1;
        pik=readFloats(filename,ncol)';

        %WD ice core chronology (from Christo Buizert, personal email)
        filename='../data/WD_chronologyunc.csv';
        ncol=3;          %depth, age, 2sig age unc
        WD=readFloats(filename,ncol);

        %interpolate WD chronology to be at 1m intervals
        [D,H,obsAge1950,obsAge1950Unc] = interpWDobs(WD);

        %Only use those depths which are different when rounded
        %[D,~] = unique(ceil(D));
        %obsAge1950=obsAge1950(ind);
        %obsAge1950Unc = obsAge1950Unc(ind); %We have unc for WD, not for Byrd    
    end
%% Height above bed param, pt every 1m
    z = (1:H)';      
end
