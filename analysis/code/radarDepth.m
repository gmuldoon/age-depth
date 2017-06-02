function [zf,depthCorrected,pikCorrected,vIce,pwErr] = radarDepth(pik,core,H)
% DESCRIPTION:
% Calculates the depth corresponding to radar reflection times for Byrd
%
% INPUT:
% - *pik* is the TWTT at the core to several interpreted horizons
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *H* is length of the core (m)
%
% OUTPUT:
% - *zf* is firn correction throughout ice column
% - *depthCorrected* is firn-corrected depth vector for entire ice column
% - *pikCorrected* is firn-corrected depth at radar horizons
% - *vIce* is random sample velocity in ice in range from Fujita et al 2000
% - *pwErr* is the error due to pulse-width effects at each horizon
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

%% Definitions
    lp=length(pik);    
    density_ice = 9.17*10^5;                  % ice density [g/m^3]
    K = 0.85;                                 % [m^3 Mg^-1]
    c = 2.98925*10^8;                         % speed of light in air [m/s] 

    vIce = randi([1.68*10^8 1.695*10^8],1,1); % random sample velocity in ice in range from Fujita et al 2000
    n_ice = c/vIce;                           % refractive index of ice

    if strcmp(core,'Byrd')                   
        pik_surf = 2.5;                       % surface TWTT, subject to the same unc as any reflector
    elseif strcmp(core,'WD')
        pik_surf = 2.5;
    else
        disp('Core value is invalid. must be Byrd or WD')
    end
    
     % Initiliaze params to help with efficiency
     depth = nan(lp,1);            % ice depth to radar horizons
     dens_int = nan(H,1);          % density integrand of firn correction
     zf = nan(H,1);                % firn correction at all depths
     depthCorrected=zeros(H,1);    % corrected depths throughout ice column 
     pikCorrected = nan(lp,1);     % corrected radar horizon depths
                 
%% Basic, everyday calculation of depth to horizons(not accounting for firn)
    for i=1:lp
        rtime=pik(1,i)-pik_surf;              % subtracting surf reflector bc need air vel above this
        time=(rtime)*10^-6/2;                 % depth=0.5*distance traveled, pik=[us]
        depth(i) = vIce*time;                 % Physics 101
    end
    
%% Load density data
    %Negative depths indicate above the surface at the time of drilling. 
    %large number is used at bottom of Byrd just to extend it far enough
    if strcmp(core,'Byrd') 
        fid = fopen('../data/byrd_densitydepth.txt');
        dat=cell2mat(textscan(fid,'%f %f'));
        density_depth=dat(:,1)';
        densityobs=dat(:,2)';
        fclose(fid);
        
        density_depth = density_depth';  % depth to each density measurement
        densityobs = densityobs';        % firn density measurements
        
        accum=0.11;
        time_correction = (2013 - 1968)*accum; % For depth comparison to WD/modern, add in recent accum
        
    elseif strcmp(core,'WD') 
        fid = fopen('../data/wd_densitydepth_every1000.txt');
        dat=cell2mat(textscan(fid,'%f %f %f %f')); % 4th column here is thickness of core sample
        density_depth=dat(:,1)';         % depth at which density was measured
        densityobs=dat(:,2)';            % density measurements 
        dens_sigma=dat(:,3)';            % uncertainty on density measurements
        fclose(fid);
        
        density_depth = density_depth';  % depth to each density measurement
        densityobs = densityobs';        % firn density measurements
        
        time_correction = 0.;            % For depth comparison to WD/modern, add in recent accum
    else
        disp('Core selection is out of range. Not sure how you managed that.')
    end
    
%% If using the WD data, create a distribution of firn corrections based on uncertainty in density (rho)
    if strcmp(core,'WD') 
        ldens = length(dens_sigma);
        dens_profile_samp = nan(ldens,1);
        for i=1:ldens
            dens_profile_err = normrnd(0,dens_sigma(i),1);                 % for each depth in the ice column, sample an error given observed 1sig error in density
            dens_profile_samp(i,1) = densityobs(i) + dens_profile_err;     % errors can be negative, centered on 0 so just add to observed to get sample of rho 
        end
        
        densityobs = dens_profile_samp';    % save this sampled density to the density param used in the interpolation below
    end
    
%% Interpolate the density profiles       
    depths=(1:H)';       % full depth vector over which to interpolate
    density_interp = interp1(density_depth', densityobs', depths(:,1));    % Interpolating density over entire core
   
%% Compute firn correction for whole column and then for vector of the picked depths only
     % Firn correction equation:
     % zf = (K/n_ice)* integral[(density_ice - density(z)) dz]

    %Integrand of firn correction at initial step 
    dens_int(1) = (density_ice - density_interp(1)*100^3)*(depths(1,1)); % for full column depth

    %Integrand of firn correction at all depths
    for i=2:H
        dens_int(i) = dens_int(i-1,1) + (density_ice - density_interp(i,1)*100^3)*... 
            (depths(i,1)-depths(i-1,1));   % multiplication converts units of densityobs/interp to match dens ice
    end

    %Firn correction throughout ice column
    for z=1:H
        zf(z,1) = K/n_ice*dens_int(z,1)/(10^6);  %10^6 converts to meters.
    end

    % Firn-corrected depth vector for entire ice column
    for i=1:H
        depthCorrected(i) = depths(i) + zf(i);  % Add firn correction bc lower density firn puffs up the ice col thickness
    end
    
%% Interpolate zf at depth of horizons 
    % Should be same for all layers bc all below firn column, but this way 
    % zf will be correct for each layer even if that changes.
    % Each element of zf corresponds to a depth in the ice column.
    
    %Interpolate zf at depths of horizons
    zfPik=interp1(depths,zf,depth(:,1));

    %Firn-corrected depth at radar horizons
    for j=1:lp
        pikCorrected(j,1) = depth(j,1) + zfPik(j,1);
    end
    
%% Incorporate uncertainty from the phase width sampling of horizons
    % Data for pik1 is digitized at 20ns. conservatively, have accuracy of
    % sampling this phase to within lambda/2. 
    pik_pw_err=nan(lp,1);
    pwErr = nan(lp,1);
    
    sigma_pw = vIce / (2 * (1/20e-9));
    surf_pw_err = normrnd(0,sigma_pw);
    for i=1:lp
        pik_pw_err(i) = normrnd(0,sigma_pw);
        pwErr(i) = sqrt(surf_pw_err^2+pik_pw_err(i)^2);
        pikCorrected(j,1) = pikCorrected(j,1)+pwErr(i);
    end
    
%% Account for additional snowball between 1968/core drilling and obs
    pikCorrected = pikCorrected + time_correction;
    
end


