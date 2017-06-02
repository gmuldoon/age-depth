function age = getAge(core,param,z,H,accumFlag)
% DESCRIPTION:
% This functions computes the age profile from a set of proposed
% parameters. For Byrd, it is based on the Schwander et al. model except for
% accumFlag = 4 which used the Morland et al. flow model
% 
% INPUT:
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *param* is ensemble of values for each parameter (size: nparam)
% - *z* is height above the bed (m)
% - *H* is length of the core (m)
% - *accumFlag* specifies the accumulation function to use, as defined in setParams.m
%
% OUTPUT:
% - *age* is the profile of simulated ages (size: H)
% 
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

%% Evaluate the ice flow model for a particular core
    if strcmp(core,'WD')
        age= wdModels(param,z,H,accumFlag);
    elseif strcmp(core,'Byrd')
        age= byrdModels(param,z,H,accumFlag);
    end
end

