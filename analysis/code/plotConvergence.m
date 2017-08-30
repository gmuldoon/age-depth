function  plotConvergence(paramDistrib,core,burnin,accumFlag,datFlag,paramRange)
% DESCRIPTION:
% This function makes convergence plots of the ice flow parameters selected
% by the Metropolis sampler
%
% INPUT:
% - *paramDistrib* is a nparam-by-totsteps ensemble of parameter values from the Metropolis algorithm
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *burnin* is number of iterations to ignore from the start of the Metropolis sampling
% - *accumFlag* specifies the accumulation function to use, as defined in setParams.m
% - *datFlag* specifies the observed chronology at Byrd to use 
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 16 Feb 2017
%
    set(0,'defaulttextinterpreter','latex')
%% Set up a loop index
    niter = length(paramDistrib(end,burnin:end));
    n = 1:niter';
    nparam = length(paramDistrib(:,1));
    lp=4;
    paramDistrib = paramDistrib(:,burnin:end);
%% Create convergence plots for each parameter 
    figure(20)
    clf
    
    subplot(4,3,1:2)
    plot(paramDistrib(1,:),'LineWidth',3)
    ylabel('q')
    set(gca,'FontSize',15)
    set(gca,'XTickLabel','')
    xlim([0 niter])
    
    subplot(4,3,3)
    hist(gca,paramDistrib(1,:))
    xlabel('q')
    set(gca,'FontSize',15)
    
    subplot(4,3,4:5)
    plot(paramDistrib(nparam-lp-2,:),'LineWidth',3)
    ylabel('h')
    set(gca,'FontSize',15)
    set(gca,'XTickLabel','')
    xlim([0 niter])
    
    subplot(4,3,6)
    hist(gca,paramDistrib(nparam-lp-2,:))
    xlabel('h')
    set(gca,'FontSize',15)
    
    subplot(4,3,7:8)
    plot(paramDistrib(nparam-lp-1,:),'LineWidth',3)
    ylabel('$v_{ice}$ (m/s)')
    set(gca,'FontSize',15)
    set(gca,'XTickLabel','')
    xlim([0 niter])
    
    subplot(4,3,9)
    hist(gca,paramDistrib(nparam-lp-1,:))
    xlabel('$v_{ice}$ (m/s)')
    set(gca,'FontSize',15)
    
    subplot(4,3,10:11)
    plot(paramDistrib(nparam-lp,:),'LineWidth',3)
    ylabel('$d_{firn}$ (m)')
    xlabel('Metropolis iteration')
    set(gca,'FontSize',15)
    xlim([0 niter])
    
    subplot(4,3,12)
    hist(gca,paramDistrib(nparam-lp,:))
    xlabel('$d_{firn}$ (m)')
    set(gca,'FontSize',15)
    
    print('../figures/convergence1','-dpng')  
 %%   
    figure(21);
    clf
    nthaccum = 11;
    set(0,'defaulttextinterpreter','latex')
    acc_depths = fliplr([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2164]);
    for i = 3:3:30
        %plot per iteration
        subplot(10,3,i-2:i-1);
        plot(paramDistrib(nthaccum,:),'LineWidth',3); hold on
        ylim([0 0.5])
        xlim([0 niter])
        ylabel('$\dot{a}_{d}$')
        set(gca, 'FontSize', 14,'yaxislocation','left')
        if i == 30
            xlabel('Metropolis iteration')
            set(gca,'XTickLabel',0:50000:niter,'XTick',0:50000:niter)
        else
            set(gca,'XTickLabel','')
        end
        % abuse the legend to show the depth range
        p = plot(NaN,NaN,'w');
        legend(p,sprintf('%i m < d < %i m',acc_depths(nthaccum),acc_depths(nthaccum-1)),'Location','northwest');
        legend('boxoff')

        %plot histogram
        subplot(10,3,i)
        hist(gca,paramDistrib(nthaccum,:))
        xlim([0 0.25])   
        set(gca, 'FontSize', 14,'YTickLabel','')
        if i == 30
            xlabel('Accumulation rate (m/a)')
        else
            set(gca,'XTickLabel','')
        end
        % go to next parameter value
        nthaccum = nthaccum-1;
    end
    
    print('../figures/convergence2','-dpng')

%%    
    
%     addpath('../../../../Matlab_functions/distinguishable_colors');
    colors=distinguishable_colors(4);
    colors(3,:) = [0 0.5 0]; %make the green darker so it's no blinding
 
    figure(22)
    clf   
    plot(paramDistrib(15,:),'Color',colors(1,:),'Linewidth',3); hold on
    plot(paramDistrib(16,:),'Color',colors(2,:),'Linewidth',3); hold on
    plot(paramDistrib(17,:),'Color',colors(3,:),'Linewidth',3); hold on
    plot(paramDistrib(18,:),'Color',colors(4,:),'Linewidth',3); hold on
    ylim([0 2164])
    xlim([0 niter])
    xlabel('Metropolis iteration')
    ylabel('Depth')
    set(gca,'YDir','reverse','Fontsize',14);

    print('../figures/convergence3','-dpng')
    
    

end

