function TWTTm_microsec = Dr_TWTT(Dr,vIce,dFirn)

    % Change normal ice column to one without firn, 
    % so can assume constant velocity in glacial ice
    Dice(:,1) = Dr(:,1) - dFirn;
    
    % Sample ePW
    % Data for pik1 is digitized at 20ns. conservatively, have accuracy of
    % sampling this phase to within lambda/2. 
    ePW = nan(length(Dr),1);
    sigma_pw = 14e-9; %this is seconds
    %sigma_pw = 1e-9;
    pik_pw_err=nan(length(Dr),1);
    surf_pw_err = normrnd(0,sigma_pw);
    for j=1:length(Dr)
        pik_pw_err(j) = normrnd(0,sigma_pw);
        ePW(j) = sqrt(surf_pw_err^2+pik_pw_err(j)^2);
    end
        
    % Classic TWTT
    TWTTice(:,1) = 2*Dice(:,1)/vIce + ePW(:,1); 
    
    % All surfaces flattened to 2.5 us; add in air time to compare to radar
    TWTTsurf = 2.5*10^-6; % seconds
    TWTTm(:,1) = TWTTice(:,1)+TWTTsurf;
    
    TWTTm_microsec = TWTTm*10^6;
end