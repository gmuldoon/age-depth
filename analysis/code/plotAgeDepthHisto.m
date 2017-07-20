function plotAgeDepthHisto2(pikDepth,pikDepthStats,core,pikAge,pikAgeStats,burnin,datFlag)
% DESCRIPTION:
% This function plot age and depth histograms for horizons
%
% INPUT:
% - *pikDepth* is an lp-by-nsamp matrix of depths of each horizon with unc
% - *pikDepthStats* contains TWTT, mean/median/stddev of each WD horizon depth
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *pikAge* is an lp-by-nsamp matrix of ages of each horizon with unc
% - *pikAgeStats* contains TWTT, mean/median/stddev of each WD horizon age
% - *burnin* is number of iterations to ignore from the start of the Metropolis sampling
% - *datFlag* specifies the observed chronology at Byrd to use 
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 16 Feb 2017
%
    %% Set up some colors and open a new figure
    colors=distinguishable_colors(length(pikAge(:,1)));
    figure
    clf
    Nmax = nan(length(pikAge(:,1)),1);
    
    %% Plot pik depth histogram (with variable depth max on plot)
    subplot(2,1,1)
    for i=1:length(pikAge(:,1))
        % set plot characteristics so that each pik is a different color
        f=histogram(gca,(pikDepth(i,:)));hold on
        f.FaceColor=colors(i,:);hold on
        f.EdgeColor=colors(i,:);hold on
        %f.BinWidth=2;hold on
        f.FaceAlpha=0.25; hold on
        f.EdgeAlpha=0.5; hold on

        % Label mean and stddev for each histo so they hopefully don't overlap
        if rem(i,2)
            text(pikDepthStats(i,2),max(f.Values)+randi([0,10])*7,strcat(num2str(round(pikDepthStats(i,2),1)),' \pm ',...
            num2str(round(pikDepthStats(i,4),1)),' m'),'Fontsize',14,'Color',colors(i,:),'FontWeight','bold')
        else
            text(pikDepthStats(i,2),max(f.Values)+randi([-10,0])*9,strcat(num2str(round(pikDepthStats(i,2),1)),' \pm ',...
            num2str(round(pikDepthStats(i,4),1)),' m'),'Fontsize',14,'Color',colors(i,:),'FontWeight','bold')
        end
        Nmax(i)=max(f.Values);
    end

    % add a title and other axis labels
    ax1=set(gca);
    if strcmp(core,'Byrd')
        title(sprintf('Depth Distrib for selected %s horizons using %s data',core,string(datFlag)))
    else
        title(sprintf('Depth Distrib for selected %s horizons',core))
    end
    xlabel('Depth (m)')
    ylabel('N')
    X1.YLim=[0 max(Nmax)+40];

    %% Plot pik age histogram (with fixed age max on plot)
    subplot(2,1,2)
    for i=1:length(pikAge(:,1))
        % set plot characteristics so that each pik is a different color
        f=histogram(pikAge(i,burnin:end)/1000); hold on
        f.FaceColor=colors(i,:);hold on
        f.EdgeColor=colors(i,:);hold on
        f.FaceAlpha=0.25; hold on
        f.EdgeAlpha=0.5; hold on
        %f.BinWidth=0.1; hold on
        
        % Label mean and stddev for each histo so they hopefully don't overlap
        if rem(i,2)
            text(pikAgeStats(i,3)/1000-0.05,max(f.Values)+randi([-10,0])*200,strcat(num2str(round(pikAgeStats(i,2)/1000,2)),...
            ' \pm ',num2str(round(pikAgeStats(i,4)/1000,2)),' ka'),'Fontsize',14,'Color',colors(i,:),'FontWeight','bold')
        else
            text(pikAgeStats(i,3)/1000-0.05,max(f.Values)+randi([0,10])*220,strcat(num2str(round(pikAgeStats(i,2)/1000,2)),...
            ' \pm ',num2str(round(pikAgeStats(i,4)/1000,2)),' ka'),'Fontsize',14,'Color',colors(i,:),'FontWeight','bold')
        end
         Nmax(i,1)=max(f.Values);
    end

    % add a title and other axis labels
    ax2=set(gca);
    if strcmp(core,'Byrd')
        title(sprintf('Age Distrib for selected %s horizons using %s data',core,string(datFlag)))
    else
        title(sprintf('Age Distrib for selected %s horizons',core))
    end
    xlabel('Age (ka)')
    ylabel('N')
    ylim([0 max(Nmax)+75]);hold on
    xlim([0 65])

end

