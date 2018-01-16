 function [sigma_ns,sigma_meters,SNR_dB] = radar_power(lp)
    %% Settings
    set(0,'defaulttextinterpreter','latex')
    addpath('../../../../Matlab_functions/distinguishable_colors');

    %% Import power data
    dat = importdata('/Users/gail/Downloads/pik1.dat');

    TWTT_microsec = dat(:,1);
    pwr_dBm = dat(:,2);

    %% Align LM surface with the power surface so can overplot reflectors
    [~,idx] = max(pwr_dBm);
    pwr_surf = TWTT_microsec(idx); %put in LM units
    LM_surf = 250/100;
    surf_diff = pwr_surf - LM_surf;


    %% Plot
    y=[0,-100];
    %   x = pik'; % not at the exact same location as pwr diagram, so don't use
    x = [ 781,958,1304,1755]/100;

%     figure(100)
%     clf
%     plot(TWTT_microsec,pwr_dBm); hold on
%     plot([4.04,4.04],[0,-100],'k-.','LineWidth',2); hold on %surface
%     plot([29.18,29.18],[0,-100],'k-.','LineWidth',2); hold on %bed
%     colors=distinguishable_colors(lp);
%     colors(3,:) = [0 0.5 0]; %make the green darker so it's no blinding
%     for i = 1:lp  
%         plot([x(i),x(i)]-surf_diff,[0,-100],'Color',colors(i,:),'LineWidth',4); hold on
%     end
%      xlim([0 35])
%      xlabel('TWTT ($\mu$s)')
%      ylabel('Power (dB)')


    %% Computing SNR
    % from plot
    signal = [-1.504,-17,-28.69,-37.18,-51.44];
    noise = [-25.24,-27.41,-36.82,-48.81,-72.66];

    SNR_dB = signal-noise;
    SNR = 10.^(SNR_dB/10);

    %% Computing range resolution
    % bandwidth = 15e6;
    v_ice = 1.685*10^8;
    % delta_r  = v_ice/(2*bandwidth);
    delta_r = 8.4; %from young et al. 2011 because not actually as good as theoretical

    %% Compute uncertainty in precision
    sigma_meters = delta_r ./ sqrt(SNR);
    sigma_ns = sigma_meters/v_ice*1e9;

    %sigmaTWTT = sigma_ns/100;

 end