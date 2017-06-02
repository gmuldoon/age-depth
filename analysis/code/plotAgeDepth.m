% NOTE: 
% plotAgeDepth is the main program file for the Byrd layers code.
%
% DESCRIPTION:
% This script does a Bayesian analysis of the Byrd ice core chronology
% and plots the results. 
%
% USAGE:
% Set the flags below  (accum, datIn, core) as appropriate. Then run.
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017
%
%% Set up the environment
% clear
close all
addpath('../Matlab_functions/shadedErrorBar');
addpath('../Matlab_functions/herrorbar');
addpath('../Matlab_functions/boundedline');

%% Set flags 
accum=7;                    % type of accum solver. (see setParams.m for options)
datIn=string('volcanic');   % observational chronology to constrain analysis
core = 'Byrd';              % ice core to derive chronology for (WD or Byrd)
depthPlotting = 1;          % whether or not to plot depth calculation

%% Do the analysis and plot
legendSym = [];
legendCnt = 0;             % number of entries in the legend
for n = 1:length(accum)         % choose accum flag (if more than one)

     accumFlag = accum(n);
    for nn = 1:length(datIn)    % choose constraining data (if more than one)
        datFlag = datIn(nn);

        % calculate age-depth for Byrd
        [paramDistrib,burnin,numsteps,ageDistrib,obsAge1950,D,H,z,...
            pikDepth,pikDepthStats,pikAge,pikAgeStats,posterior]...
            = calcAgeDepth(accumFlag,datFlag,core,depthPlotting);

        % plot histograms of age and depth for Byrd
        plotAgeDepthHisto(pikDepth,pikDepthStats,core,pikAge,pikAgeStats,burnin,datFlag)
        
        % plot convergence of model parameters
        plotConvergence(paramDistrib,core,burnin,accumFlag,datFlag)

        % plot age-depth profile for Byrd
        [p,lab] = plotSpaghettiEnvelope(ageDistrib,obsAge1950,D,H,accumFlag,z,burnin,paramDistrib,datFlag,0); hold on
        legendLab{legendCnt+1} = lab; 
        legendCnt = legendCnt+1; %increment legend entries
        legendSym= [legendSym,p];
        if strcmp('oxygen',datFlag)
            D_o2 = D;
            obsAge1950_o2 = obsAge1950;
            paramDistrib_o2 = paramDistrib;
        elseif strcmp('volcanic',datFlag)
            D_volc = D;
            obsAge1950_volc = obsAge1950;
            paramDistrib_volc = paramDistrib;
        end
        
        % Overplot observed data age-depth profile
        [p,lab] = plotSpaghetti(ageDistrib,numsteps,burnin,obsAge1950,D,H,accumFlag,z,datFlag); hold on
        legendLab{legendCnt+1} = lab; 
        legendCnt = legendCnt+1; %increment legend entries
        legendSym= [legendSym,p];        
    end
end

% Plot obs data on age-depth line plot

%    [p,lab] = plotSpaghettiEnvelope(ageDistrib,obsAge1950,D,H,accumFlag,z,burnin,paramDistrib,datFlag,1);
%    legendLab{nn+1}
%p1=plot(obsAge1950_volc/1000,D_volc,'k.','LineWidth',2,'MarkerSize',14); hold on
%legendSym= [legendSym,p1];
legend(legendSym,legendLab)
legend('boxoff')

%% plot histograms of WD piks
%clear D, obsAge1950
[WDpikAge,WDpikDepth,WDpikDepthStats,WDpikAgeStats]=calcWDLayerAge(numsteps,depthPlotting);
plotAgeDepthHisto(WDpikDepth,WDpikDepthStats,'WD',WDpikAge,WDpikAgeStats,burnin,datFlag)

%%
% Plot the results
% % plotAgeDepthError(ageDistrib,z,obsAge1950,D,pikAgeStats,pikDepthStats,accumFlag,H,fullAgeStats,burnin)
% plotAgeDepthHisto(pikDepthUnc,pikDepthStats,core,pikAge,pikAgeStats,burnin,datFlag)
% plotCostHisto(posterior)

%if accumFlag == 100
%plotSpaghetti(ageDistrib,numsteps,obsAge1950,D,H,accumFlag,z,burnin)
%plotAgeDepthError(ageDistrib,z,obsAge1950,D,pikAgeStats,pikDepthStats,accumFlag,H,fullAgeStats,burnin)

