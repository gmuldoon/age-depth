function pikAge = calcPikAge(pikDepthSamp,ageDistrib,z,core,H,obsAge1950,obsAge1950Unc)
% DESCRIPTION:
%
% NOTE:
% INPUT:
% - *pikDepth* is the lp-by-numsteps ensemble of depths for each horizon
% - *ageDistrib* is an H-by-totsteps matrix of ages for each meter depth in the core
% - *z* is height above the bed (m)
% - *numsteps* is number of Metropolis samples after burnin
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *H* is length of the core (m)
% - *obsAge1950* is observed chronology (years before 1950)
% - *obsAge1950Unc* is observed chronology (years before 1950)
%
% OUTPUT:
% - *pikAge* is the lp-by-numsteps ensemble of ages for each horizon
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

disp('Computing horizon ages with uncertainty')
%% Define some things
    pikZ = H - pikDepthSamp; %flipping the coordinate axis
    pikAge = nan(size(pikDepthSamp)); % initializing matrix for horizon ages
    m = length(pikDepthSamp(1,:)); % numsteps 
%% Compute age of horizons given a sampled ice flow model   
    for n=1:m %numsteps for Byrd, nsamp for WD
        if strcmp(core,'Byrd')
            % sample an age model (ie a set of flow parameters)
            samp_age = randi([1 m]); 

            % evaluate model at each pikZ for that model to find corresponding age
            [~,ind] = intersect(z,round(pikZ(:,n)));
            pikAge(:,n) = ageDistrib(ind,samp_age);
            
            %flip age bc ordered by z which puts deeper ages first
            pikAge(:,n) = flipud(pikAge(:,n));
            
        elseif strcmp(core,'WD')
            % evaluate model at each pikZ for that model to find corresponding age
            [~,ind] = intersect(z,round(pikZ(:,n)));
            
            %assume WD ages are observed w/ gaussian errors 
            pikAge(:,n) = normrnd(obsAge1950(ind),obsAge1950Unc(ind));
        end
    end
end

