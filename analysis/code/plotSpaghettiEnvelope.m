
function [p,lab] = plotSpaghettiEnvelope(ageDistrib,obsAge1950,D,H,accumFlag,z,burnin,paramDistrib,datFlag,plotData)
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

    p = boundedline(mean(ageDistrib(:,end-1500:end),2)/1000,flipud(z),2*std(ageDistrib(:,end-1500:end),1,2)/1000,'orientation', 'horiz','alpha',shadecolor,'transparency', 0.1); hold on
    

    %% Plot data over top
    if plotData
        plot(obsAge1950/1000,D,'k.','LineWidth',2,'MarkerSize',14); hold on
        lab = sprintf('Observed %s ages',datFlag);
    end
    %% Set plot params
    set(gca,'YDir','reverse');
    xlabel('Age (ka)','Fontsize',14)
    ylabel('Depth (m)','Fontsize',14)
    set(gca, 'FontSize', 16)
    axis([0 65 0 H])

end

