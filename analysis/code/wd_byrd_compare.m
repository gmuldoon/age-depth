addpath('../Matlab_functions/herrorbar');

%depth from histogram
byrd_depth = [152.8, 379.9, 466.8, 574.9, 631.4, 722.5, 866, 1140.6];
byrd_depthstd = [0.4, 1, 1.2, 1.5, 1.7, 1.9, 2.3, 6.9];

%from mean of pikDepthUnc
wd_depth = [586.28, 1107.79, 1284.15, 1390.48, 1561.78, 1708.6, 1891.73, 2417.46];
%from depth hist
wd_depthstd = [1.5,2.84,3.29, 3.57, 4.01, 4.38, 4.85, 6.2];


%byrd age from chronology (ka)
byrd_age = [1120, 2980, 3770, 4830, 5500, 6710, 8820, 13440];
byrd_agestd = [30, 70, 80, 90, 90, 90, 120, 160];

%wd age from chronology for these particular piks

wd_age = [2386, 4906, 5908, 6595, 7856, 9116, 10908, 17489];
wd_agestd = [4.7, 19.2, 25, 29.1, 36.5, 43.9, 54.5, 174.9];

figure(20)
clf
%plot(byrd_age, wd_age); hold on
herrorbar(byrd_age,wd_age,2*wd_agestd); hold on
errorbar(byrd_age,wd_age,2*byrd_agestd); hold on
refline(1,0); hold on
xlabel('WAIS Divide age (ka)')
ylabel('Byrd age (ka)')
legend('WD Age error','Byrd age error', 'Identity line')

