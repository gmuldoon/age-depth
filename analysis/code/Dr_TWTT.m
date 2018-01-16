function TWTTm_microsec = Dr_TWTT(Dr,vIce,sigma_ns,sigma_dens)

    % Change normal ice column to one without firn, 
    % so can assume constant velocity in glacial ice
    K = 0.85;                                 % [m^3 Mg^-1]
    c = 2.98925*10^8;                         % speed of light in air [m/s] 
    vIce = randi([1.68*10^8 1.695*10^8],1,1); % random sample velocity in ice in range from Fujita et al 2000
    n_ice = c/vIce; 
    
    
    dFirn = 8.6621; %from radardepth
    sigma_firn =  K/n_ice*sigma_dens/(10^6); %assumes all reflectors are below firn layer
    dFirn_sample = dFirn + sigma_firn;
    
    Dice(:,1) = Dr(:,1) - dFirn_sample;
    
    % Sample ePW
    % Data for pik1 is digitized at 20ns. conservatively, have accuracy of
    % sampling this phase to within lambda/2. 
    ePW = nan(length(Dr),1);
%     sigma_pw = 14e-9; %this is seconds
    %sigma_pw = 1e-9;
     sigma_pw = sigma_ns*1e-9;
%     sigma_pw = ones(5,1)*14e-9;
    pik_pw_err=nan(length(Dr),1);
    surf_pw_err = normrnd(0,sigma_pw(1));
    for j=1:length(Dr)
        pik_pw_err(j) = normrnd(0,sigma_pw(j+1));
        ePW(j) = sqrt(surf_pw_err^2+pik_pw_err(j)^2);
    end
        
    % Classic TWTT
    TWTTice(:,1) = 2*Dice(:,1)/vIce + ePW(:,1); 
    
    % All surfaces flattened to 2.5 us; add in air time to compare to radar
    TWTTsurf = 2.5*10^-6; % seconds
    TWTTm(:,1) = TWTTice(:,1)+TWTTsurf;
    
    TWTTm_microsec = TWTTm*10^6;
end