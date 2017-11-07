dat = importdata('/Users/gail/Downloads/pik1.dat');

TWTT_microsec = dat(:,1);
pwr_dBm = dat(:,2);


% figure(100)
% clf
% plot(pwr_dBm)

pwr_sort = sort(pwr_dBm,'descend');
%Looking for noise level
%Examples from different peaks:
%[-65.33,-70.88]
%[-60.38,-62.28]
%[-44.37,-45.08]
%[-77.81,-76.1]
%[-61.24,-58.19]
%[-54.65,-60.53]
%[-54.87,-48.93]
%[-76.58,-85.6]
%[-59.28,-72.27]
%[-65.33,-70.88]
noise = mean([-65.33,-70.88]);

%Looking for signal (peak)
%find average of top 100 values
%-53.76
%-48.12
%-36.39
%-64.89
%-46.95
%-43.15
%-37.31
%-68.88
%-51.64
%-54.2
signal = -54.2;

%Computing SNR
SNR_dB = signal-noise;

%Computing range resolution
% bandwidth = 15e6;
% v_ice = 1.685*10^8;
% delta_r  = v_ice/(2*bandwidth);
delta_r = 8.4; %from young et al. 2011 because not actually as good as theoretical

%Compute uncertainty in precision
sigma_meters = delta_r / sqrt(SNR_dB);
sigma_ns = sigma_meters/v_ice*1e9;


%sigma_ns = 
%13.1622
%13.7160
%17.2674
%14.3521
%13.950
%13.1189
%13.0512
%14.2666
%13.2596
%13.3688

sigma_ns_list = [13.1622, 13.7160, 17.2674,14.3521, 13.950, 13.1189, ...
    13.0512, 14.2666, 13.2596, 13.3688];

%mean_sigma_ns = mean(sigma_ns_list);
% ~ 14 ns