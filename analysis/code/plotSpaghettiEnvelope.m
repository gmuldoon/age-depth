
function [p,lab,pikDepthUnc,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(ageDistrib,obsAge1950,D,H,accumFlag,z,burnin,paramDistrib,datFlag,plotData,Ar,Zr)
% DESCRIPTION:
% This function make a line plot of age-depth with an envelope the width of
% 2std on either side of the mean.
%
% INPUT:
% - *ageDistrib* is an H by totsteps matrix of ages for each meter depth in the core
% - *obsAge1950* is observed chronology (years before 1950)
% - *D* is observed chronology depth (m)
% - *H* is length of the core (m)
% - *accumFlag* specifies the accumulation function to use, as defined in setParams.m
% - *z* is height above the bed (m)
% - *burnin* is number of iterations to ignore from the start of the Metropolis sampling
% - *paramDistrib* is a nparam-by-totsteps ensemble of parameter values from the Metropolis algorithm
% - *datFlag* specifies the observed chronology at Byrd to use 
% 
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017
    addpath('/Users/gail/Documents/Research/Matlab_functions/boundedline/');
    figure(10)

    %% Assign colors to each accumFlag so they can be distinguished if plotted together
    if accumFlag == 7
        shadecolor='r--';
        lab = sprintf('Simulated %s record',string(datFlag));
    end

    %% Plot the age-depth envelope
    %p = boundedline(mean(ageDistrib(:,end-1500:end),2)/1000,flipud(z),0.03*mean(ageDistrib(:,end-1500:end),2)/1000,'orientation', 'horiz','alpha','b','transparency', 0.1); hold on

    p = boundedline(mean(ageDistrib(:,burnin:end),2)/1000,flipud(z)/H,2*std(ageDistrib(:,burnin-1500:end),1,2)/1000,'orientation', 'horiz','alpha',shadecolor,'transparency', 0.1); hold on
    

    %% Plot data over top
    if plotData
        p2 = plot(obsAge1950/1000,D/H,'ko','LineWidth',2,'MarkerSize',7); hold on
        lab = sprintf('Observed %s ages',datFlag);
    end
    
    %% Plot reflector age/depth
    p3 = plot(mean(Ar/1000,2),mean((H-Zr)/H,2),'b^','LineWidth',2,'MarkerSize',7); hold on
    
    %% Plot WAIS Divide chronology
    filename='../data/WD_chronologyunc.csv';
    ncol=3;          %depth, age, 2sig age unc
    WD=readFloats(filename,ncol);
    [~,Hwd,wdAge1950,wdAge1950Unc] = interpWDobs(WD);

    p4 = plot(wdAge1950/1000,(1:Hwd)/Hwd,'g-'); hold on
    
    
    %% Plot reflectors at WD
    %TWTT for traced radar horizons
    zwd = 1:Hwd;
    %TWTT for traced radar horizons
    filename='../data/WD_TWTT_5.txt';
    pik=readFloats(filename,1)';
    [pikDepthUnc,~,~]=twtt2depth(pik,3404,'WD',0,159,zwd);
    WD_depth = mean(pikDepthUnc,2);
    wd_age = wdAge1950(round(WD_depth));

    plot(wd_age/1000,WD_depth/Hwd,'b^','LineWidth',2,'MarkerSize',7); hold on
    %% Set plot params
    set(gca,'YDir','reverse');
    xlabel('Age (ka)','Fontsize',10)
    ylabel('Normalized depth','Fontsize',10)
    set(gca, 'FontSize', 10)
    axis([0 65 0 1])
    pbaspect([1 1 1])
    
    legend([p,p4,p2,p3],'Estimated Byrd chronology','WD chronology','Volcanic chronology','Observed Reflectors')

end

