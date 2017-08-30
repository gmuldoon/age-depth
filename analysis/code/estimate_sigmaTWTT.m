function sigmaTWTT = estimate_sigmaTWTT(seed)
    rng(seed);
    %some convenient estimates of Dr
    Dr = [500,600,900,1300,1700]';
    nDr = length(Dr);
    %% Sample uncertainty in vIce, ePW, and dFirn; Calculate resulting TWTT  
    nsamp = 1000;
    TWTTm_microsec = nan(nsamp,nDr);
    for i = 1:nsamp    
        
        % Sample vIce
        vIce = randi([1.68*10^8 1.695*10^8],1,1); % random sample velocity in ice in range from Fujita et al 2000
        
        % Sample firn correction 
        dFirn = randi([40 100],1,1); %guesstimated range [m]

        % Model TWTT for those sampled values & with noise from pulse width
        TWTTm_microsec(i,:) = Dr_TWTT(Dr,vIce,dFirn);
    end
    
    TWTTm_mean = mean(TWTTm_microsec);
    
    
    perfect_model_discrep = (TWTTm_microsec - TWTTm_mean).^2;
    costTWTT_star = sum(perfect_model_discrep,2)/2.; %assumes the sigma here is 1
    
    %% According to Charles's GMDD paper:
    % 1/sigmaTWTT = sqrt(ke/2)/sigma_costTWTT_star
    
    ke = nDr/2;
    sigmaTWTT = std(costTWTT_star)/sqrt(ke/2);
end