function plotAgeDepthHisto2(Dr,core,Ar,burnin,datFlag,pikDepthUncWD,wdAge1950,wdAge1950Unc)
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
    set(0,'defaulttextinterpreter','latex')
    addpath('../../../../Matlab_functions/distinguishable_colors');
    %% Set up some colors and open a new figure
    colors=distinguishable_colors(length(Ar(:,1)));
    colors(3,:) = [0 0.5 0]; %make the green darker so it's no blinding
    figure(11)
    clf
    Nmax = nan(length(Ar(:,1)),1);

%     %% Plot pik depth histogram (with variable depth max on plot)

    subplot(2,1,1)

    for i=1:length(Ar(:,1))
        % set plot characteristics so that each pik is a different color
        f=histogram(gca,(Dr(i,burnin:end)));hold on
        f.FaceColor=colors(i,:);hold on
        f.EdgeColor=colors(i,:);hold on
        f.BinWidth=1;hold on
        %f.FaceAlpha=0.25; hold on
        %f.EdgeAlpha=0.75; hold on

        % Label mean and stddev for each histo so they hopefully don't overlap
        if rem(i,2)
            text(mean(Dr(i,burnin:end))*1.02,max(f.Values)*0.93,strcat(num2str(round(mean(Dr(i,burnin:end)),1)),' $\pm$~ ',...
            num2str(round(std(Dr(i,burnin:end)),1)),' m'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        else
            text(mean(Dr(i,burnin:end))*1.02,max(f.Values)*0.93,strcat(num2str(round(mean(Dr(i,burnin:end)),1)),' $\pm$~ ',...
            num2str(round(std(Dr(i,burnin:end)),1)),' m'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        end
%         text(4000,mean(Dr(i,burnin:end))+50,strcat(num2str(round(mean(Ar(i,burnin:end))/1000,2)),...
%             ' $\pm~$ ',num2str(round(std(Ar(i,burnin:end))/1000,3),'%2.2f'),' ka'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        Nmax(i)=max(f.Values);
    end

    % add a title and other axis labels
    ax1=set(gca);
%     if strcmp(core,'Byrd')
%         title(sprintf('Depth Distrib for selected %s horizons',core))
%     else
%         title(sprintf('Depth Distrib for selected %s horizons',core))
%     end
    xlabel('Depth (m)','FontSize',15)
    ylabel('N')
    X1.YLim=[0 max(Nmax)+40];
    set(gca,'FontSize',15)
    xlim([200 1700])
    text(0.01, 0.9, 'A) Byrd Reflector Depth','FontSize',15,'FontWeight','bold','Units','Normalized')
%     view(90, -90)


    %% Plot pik age histogram (with fixed age max on plot)
    
    %Read WD stuff 
    filename='../data/WD_chronologyunc.csv';
    ncol=3;          %depth, age, 2sig age unc
    WD=readFloats(filename,ncol);
    [~,Hwd,wdAge1950,wdAge1950Unc] = interpWDobs(WD);
    
    zwd = 1:Hwd;
    filename='../data/WD_TWTT_4.txt';
    pik=readFloats(filename,1)';
    [pikDepthUnc,~,~]=twtt2depth(pik,3404,'WD',0,159,zwd);
    WD_depth = mean(pikDepthUnc,2);
    wd_age = wdAge1950(round(WD_depth));
    
    subplot(2,1,2)
    for i=1:length(Ar(:,1))
        % set plot characteristics so that each pik is a different color
        f=histogram(Ar(i,burnin:end)/1000); hold on
        f.FaceColor=colors(i,:);hold on
        f.EdgeColor=colors(i,:);hold on
        %f.FaceAlpha=0.75; hold on
        %f.EdgeAlpha=0.75; hold on
        f.BinWidth=0.01; hold on
        Nmax(i,1)=max(f.Values);
        
        % Label mean and stddev for each histo so they hopefully don't overlap
        if rem(i,2)
            text(mean(Ar(i,burnin:end))/1000*1.05,max(f.Values)*0.93,strcat(num2str(round(mean(Ar(i,burnin:end))/1000,2)),...
            ' $\pm~$ ',num2str(round(std(Ar(i,burnin:end))/1000,3),'%2.2f'),' ka'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        else
            text(mean(Ar(i,burnin:end))/1000*1.05,max(f.Values)*0.93,strcat(num2str(round(mean(Ar(i,burnin:end))/1000,2)),...
            ' $\pm~$ ',num2str(round(std(Ar(i,burnin:end))/1000,3),'%2.2f'),' ka'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        end
        text(wd_age(i)/1000-0.5,max(f.Values)+50,'*','Fontsize',25,'Color','m','FontWeight','bold')
        
%         text(1000,mean(Ar(i,burnin:end))/1000-2,strcat(num2str(round(mean(Ar(i,burnin:end))/1000,2)),...
%             ' $\pm~$ ',num2str(round(std(Ar(i,burnin:end))/1000,3),'%2.2f'),' ka'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
    end

    % add a title and other axis labels
    ax2=set(gca);
%     if strcmp(core,'Byrd')
%         title(sprintf('Age Distrib for selected %s horizons',core))
%     else
%         title(sprintf('Age Distrib for selected %s horizons',core))
%     end
    xlabel('Age (ka)','FontSize',15)
    ylabel('N')
%     ylim([0 max(Nmax)+75]);hold on
    xlim([0 30])
    set(gca,'FontSize',15)
    text(0.01,0.9, 'B) Byrd Reflector Age','FontSize',15,'FontWeight','bold','Units','Normalized')
    
    text(0.5,0.79,'* ','Color','m','Units','Normalized','FontSize',25)
    text(0.52,0.82,'Age at WAIS Divide','Units','Normalized','FontWeight','bold','FontSize',15);
    rectangle('Position',[14.7 1550 4.5 200])
%     view(90, -90)
    %% Plot WD age distributions
%     WD_depth = mean(pikDepthUncWD,2);
%     wd_age = wdAge1950(round(WD_depth),:);
%     wd_age_unc = wdAge1950Unc(round(WD_depth),:);
    
%     subplot(3,1,3)
%     for i=1:length(wd_age_unc(:,1))
%         x = [0:0.005:65];
%         norm = normpdf(x,wd_age(i)/1000,wd_age_unc(i)/1000);
%         f = area(x,norm); hold on
%         % set plot characteristics so that each pik is a different color
% %         f=histogram(wd_age_unc(i,burnin:end)/1000); hold on
%           f.FaceColor=colors(i,:);hold on
%           f.EdgeColor=colors(i,:);hold on
% %         %f.FaceAlpha=0.75; hold on
% %         %f.EdgeAlpha=0.75; hold on
% %         f.BinWidth=0.01; hold on
% %         
% %         % Label mean and stddev for each histo so they hopefully don't overlap
%         if rem(i,2)
%             text(wd_age(i)/1000+0.25,max(f.YData)+2,strcat(num2str(round(wd_age(i)/1000,2)),...
%             ' $\pm~$ ',num2str(round(2*wd_age_unc(i)/1000,2)),' ka'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
%         else
%             text(wd_age(i)/1000+0.25,max(f.YData)+2,strcat(num2str(round(wd_age(i)/1000,2)),...
%             ' $\pm~$ ',num2str(round(2*wd_age_unc(i)/1000,2)),' ka'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
%         end
% %          Nmax(i,1)=max(f.Values);
%     end
% 
%     % add a title and other axis labels
% %     title('Age Distrib for selected WAIS Divide horizons')
%     xlabel('Age (ka)','FontSize',15)
%     ylabel('N')
% %     ylim([0 max(Nmax)+75]);hold on
%      xlim([0 30])
%      set(gca,'FontSize',15)
%      text(0.01, 0.8, 'C) WAIS Divide Age','FontSize',15,'FontWeight','bold','Units','Normalized')
     
     print('../figures/agedepthhisto','-dpng')

end

