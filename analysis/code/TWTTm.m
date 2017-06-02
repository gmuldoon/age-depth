function TWTT = TWTTm(Param)
%% Read proposed Dr
Dr = Param(end-lp:end,1);

%% Load density data
    %Negative depths indicate above the surface at the time of drilling. 
    %large number is used at bottom of Byrd just to extend it far enough
    
        fid = fopen('../data/byrd_densitydepth.txt');
        dat=cell2mat(textscan(fid,'%f %f'));
        density_depth=dat(:,1)';
        densityobs=dat(:,2)';
        fclose(fid);
        
        density_depth = density_depth';  % depth to each density measurement
        densityobs = densityobs';        % firn density measurements
        
        % Correct time difference between ice core drilling, radar collection
        time_correction = (2013 - 1968)*Param(5,1); 
               
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
        
%% Interpolate zf at depth of horizons 
    % Should be same for all layers bc all below firn column, but this way 
    % zf will be correct for each layer even if that changes.
    % Each element of zf corresponds to a depth in the ice column.
    Dr_firncorrected = nan(lp,1);
    
    %Interpolate zf at values of Dr
    zfPik=interp1(depths,zf,Dr);

    %Firn-corrected depth at radar horizons
    for j=1:lp
        Dr_firncorrected(j,1) = depth(j,1) + zfPik(j,1);
    end
    
%% Account for time between 1968 and data collection
    Dr_timecorrected = Dr_firncorrected - time_correction;
    
 %% Basic equation   
    vIce = 1.68*10^8; %m/s
    TWTT = 2*Dr_timecorrected*vIce;
    
end