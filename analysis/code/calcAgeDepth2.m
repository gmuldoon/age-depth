%function [Ar, Param]=calcAgeDepth2(accumFlag,datFlag,core,plotting)
tic;
accumFlag=7;                    % type of accum solver. (see setParams.m for options)
datFlag=string('volcanic');   % observational chronology to constrain analysis
core = 'Byrd';              % ice core to derive chronology for (WD or Byrd)
depthPlotting = 1;          % whether or not to plot depth calculation

%% 1. Set up environment
    addpath('/Users/gail/Documents/Research/Matlab_functions');
    addpath('/Users/gail/Documents/Research/Matlab_functions/export_fig');
    set(0,'DefaultFigureWindowStyle','docked')
    seed=123;      % seed for reproducibility

%% 2. Load data
    [H,D,z,pik,obsAge1950]=loadCorePikData2(core,datFlag);
    
%% 3. Set priors for QOI
    [paramRange,nparam] = setParams2(accumFlag,H,pik);
    
%% 4. Evaluate proposed parameters with likelihoods on Aic, TWTTr
    [Param,numsteps,cost,burnin,S,reg]=metropolisAgeSampler2(D,z,H,obsAge1950,accumFlag,seed,paramRange,nparam,pik);
    
%% 5. Compute Ar profile from accepted parameter values
    age = nan(H,numsteps+burnin);
    Zr = H-Param(nparam-length(pik)+1:nparam,:); %change to elevation from depth
    Ar = nan(size(Zr));
    lp = length(pik);
    
    for i = 1:burnin+numsteps
        % Full age-depth profile
        age(:,i) = byrdModels2(Param(:,i),z,H,accumFlag,lp);
        
        % Age of radar horizons
        for j = 1:length(pik)
            Ar(j,i) = byrdModels2(Param(:,i),z,H,accumFlag,lp,Zr(j,i));
        end
    end
    
%% THE REST IS PLOTTING RESULTS    
%%  Plot convergence of parameters
% plotConvergence(Param,core,burnin,accumFlag,datFlag,paramRange)    

%%  Make spaghetti plot
figure(10)
clf
[~,~,pikDepthUncWD,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(age,obsAge1950,D,H,accumFlag,z,burnin,Param,datFlag,1,Ar,Zr,1);
% print('../figures/spaghetti','-dpng','-r1000')
% savefig('../figures/spaghetti.fig')
%% Plot age-depth histograms
plotAgeDepthHisto2(Param(nparam-lp+1:nparam,:),core,Ar,burnin,datFlag,pikDepthUncWD,wdAge1950,wdAge1950Unc)
savefig('../figures/agedepthhisto.fig')
print('../figures/agedepthhisto','-dpng','-r1000')
 %% Regularization plot
set(0,'defaulttextinterpreter','latex')
figure(12)
plot(reg(burnin:end),'LineWidth',3);
% hist(reg(burnin:end))
xlabel('Metropolis iteration','FontSize', 15)
ylabel('$r$')
set(gca, 'FontSize', 15)
%xlim([0 numsteps])
% print('../figures/regularization','-dpng','-r1000')

%% Cost plots
set(0,'defaulttextinterpreter','latex')
figure(13)
clf
subplot(2,1,1)
plot(cost(burnin:end,1),'r','LineWidth',3);hold on % TWTT
set(gca,'XTickLabel','')
ylabel('$cost_{TWTT}$','FontSize', 15)
xlim([0 numsteps])
set(gca, 'FontSize', 15)

subplot(2,1,2)
plot(cost(burnin:end,2),'b','LineWidth',3);hold on  % Age
xlabel('Metropolis iteration','FontSize', 15)
ylabel('$cost_{Age}$','FontSize', 15)
xlim([0 numsteps])
set(gca, 'FontSize', 15)
% print('../figures/cost','-dpng','-r1000')

%% Accumulation depth plot
           
acc_depths = fliplr([ 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]);
figure(14)
clf

for i = burnin:1000:burnin+numsteps
    %subplot(2,1,1)
     plot(acc_depths,Param(2:nparam-lp-3,i)); hold on
     smoothed = smooth(Param(2:nparam-lp-3,i),3);
     plot(acc_depths,smoothed, '-.','LineWidth',3); hold on

end

set(gca, 'FontSize', 12)
% print('../figures/accumdepth','-dpng','-r1000')

%% Accumulation paper figure
%sort accums by age cost
[cost_sorted,cost_sorted_I] = sort(cost(:,2),1); 
Param_sorted = Param(:,cost_sorted_I);

acc_depths = fliplr([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]);
naccum = length(acc_depths);
AccumSortedFull = zeros(H,burnin+numsteps);
for i = 1:naccum
    if i == 1 % bottom of the ice
        for j = H:-1:acc_depths(1)
            AccumSortedFull(j,:) = Param_sorted(i+1,:);
        end
    elseif i == naccum % top of the ice
        for j = acc_depths(i):-1:1
            AccumSortedFull(j,:) = Param_sorted(i+2,:);
        end
        for j = acc_depths(i-1):-1:acc_depths(i)
            AccumSortedFull(j,:) = Param_sorted(i+1,:);
        end
    else
        for j = acc_depths(i-1):-1:acc_depths(i)
            AccumSortedFull(j,:) = Param_sorted(i+1,:); 
        end
    end
end
set(0,'defaulttextinterpreter','latex')
figure(15)
clf
n=20000;
cm=colormap(parula((burnin+numsteps)/n+1));
cnt = 0;
for i = burnin:n:burnin+numsteps
    cnt = cnt+1;
    plot(smooth(AccumSortedFull(:,i),50),'color',cm(cnt,:),'LineWidth',3); hold on
    %accVar(i) = var(Param(2:nparam-lp-3,i));    
end

xlabel('Depth (m)','FontSize',18)
ylabel('$\dot{a}$ (m/a)','FontSize',18)
xlim([0 2164]);
h=colorbar;
set(get(h,'title'),'string','Normalized Cost');
set(gca,'FontSize',18);

% Add tick marks for Age

ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'));
set(ax2, 'XTick', get(ax1, 'XTick'));
AgeTickLabels= {'0' '1.4' '3.3' '5.4' '7.7' '10.7' '14.8' '21.8' '30.9' '42.1' '63.7'};
DummyYLabels = {''};
set(ax2, 'XTickLabel', AgeTickLabels,'YTickLabel',DummyYLabels);
xlabel('Age (ka)');
set(gca,'FontSize',18);

% savefig('../figures/accumdepthSorted.fig')
% print('../figures/accumdepthSorted','-dpng','-r1000')


%% Degrees of freedom plots
% agediv4 = age;
% Paramdiv4 = Param;
% Ardiv4 = Ar;
% Zrdiv4 = Zr;
% Drdiv4 = Paramdiv4(end-lp:end,:);
% [~,~,pikDepthUncWD,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(...
%     agediv4,obsAge1950,D,H,accumFlag,z,burnin,Paramdiv4,datFlag,1,Ardiv4,Zrdiv4,4);

% agediv2 = age;
% Paramdiv2 = Param;
% Ardiv2 = Ar;
% Zrdiv2 = Zr;
% Drdiv2 = Paramdiv2(end-lp:end,:);

%[~,~,pikDepthUncWD,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(...
%   agediv2,obsAge1950,D,H,accumFlag,z,burnin,Paramdiv2,datFlag,1,Ardiv2,Zrdiv2,2);

% agediv1 = age;
% Paramdiv1 = Param;
% Ardiv1 = Ar;
% Zrdiv1 = Zr;
% Drdiv1 = Paramdiv1(end-lp:end,:);
% [~,~,pikDepthUncWD,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(...
%     agediv1,obsAge1950,D,H,accumFlag,z,burnin,Paramdiv1,datFlag,1,Ardiv1,Zrdiv1,1);

% agediv5 = age;
% Paramdiv5 = Param;
% Ardiv5 = Ar;
% Zrdiv5 = Zr;
% Drdiv5 = Paramdiv5(end-lp:end,:);
% [~,~,pikDepthUncWD,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(...
%     agediv5,obsAge1950,D,H,accumFlag,z,burnin,Paramdiv5,datFlag,1,Ardiv5,Zrdiv5,5);

% agediv10 = age;
% Paramdiv10 = Param;
% Ardiv10 = Ar;
% Zrdiv10 = Zr;
% Drdiv10 = Paramdiv10(end-lp:end,:);

%%
% set(0,'defaulttextinterpreter','latex')
% figure(30);
% clf
% % text(0.25,0,'Depth (m)')
% % text(0.75,0, 'Age (ka)')
% % 
% subplot(4,1,1)
% [ax1] = plotAgeDepthHisto2_ke(Drdiv1,core,Ardiv1,burnin,1);
% 
% subplot(4,1,2)
% %subplot(4,2,4)
% [ax2] = plotAgeDepthHisto2_ke(Drdiv2,core,Ardiv2,burnin,2);hold on
% 
% subplot(4,1,3)
% [ax3] = plotAgeDepthHisto2_ke(Drdiv4,core,Ardiv4,burnin,3);
% % 
% subplot(4,1,4)
% [ax4] = plotAgeDepthHisto2_ke(Drdiv5,core,Ardiv5,burnin,4);

% print('../figures/keCompare','-dpng')

%% Correlation plots
set(0,'defaulttextinterpreter','latex')
% Write out to file 
%acc_depths = {'$< 200 m$', '$< 400 m$', '$< 600 m$', '$< 800 m$ ', '$< 1000 m$ ', '$< 1200 m$', '$< 1400 m$', '$< 1600 m$', '$< 1800 m$', '$> 1800 m$'}';
% labels = {'\boldmath{$q$}','\boldmath{$\dot{a}_{<\,1400\, m}$}', '\boldmath{$\dot{a}_{<\,1600\,m}$}', '\boldmath{$\dot{a}_{<\,1800\, m}$}', '\boldmath{$\dot{a}_{>\, 1800\, m}$}'}';
%n = length(labels);
n=11;
tmp = Param(1:11,:);


figure(16)
clf
%[Scatter,SAx,BigAx,Histo,HAx]= plotmatrix_withr(Param(12:-1:1,burnin:1000:end)'); % 11 is the shallowest
[Scatter,SAx,BigAx,Histo,HAx]= plotmatrix_withr(tmp(1:end,burnin:100:end)'); 

for i = 1:n
    %SAx(i,1).YLabel.String=string(labels(i));
    SAx(i,1).YLabel.Rotation=0;
    SAx(i,1).YLabel.FontSize=20;
    SAx(i,1).YLabel.Units='Normalized';
    SAx(i,1).YLabel.Position=[-0.3 0.5000 0];
    for j = 1:n
        if j > 1 && i > 1  
%             SAx(i,j).XLim = [0.0 0.25];
%             SAx(i,j).YLim = [0.0 0.25]; 
        end
        
        %SAx(end,j).XLabel.String = string(labels(j));
        SAx(end,j).XLabel.FontWeight = 'bold';
        SAx(end,j).XLabel.FontSize=20;
    end
    if i > 1
%         HAx(i).XLim = [0.0 0.25];
        Histo(i).BinWidth = 0.005;
    end
end

 print('../figures/accumCorrelation','-dpng','-r1000')
% print('-depsc2', '-loose', '../figures/accumCorrelation.eps');
% saveas(gcf,'../figures/accumCorrelation.png')
%  export_fig ../figures/accumCorrelation.png
%% Depth correlation
set(0,'defaulttextinterpreter','latex')
% Write out to file 
%acc_depths = {'$< 200 m$', '$< 400 m$', '$< 600 m$', '$< 800 m$ ', '$< 1000 m$ ', '$< 1200 m$', '$< 1400 m$', '$< 1600 m$', '$< 1800 m$', '$> 1800 m$'}';

%labels={'\boldmath{$d_{firn}$}','\boldmath{${D_1}$}','\boldmath{$D_2$}','\boldmath{$D_3$}','\boldmath{$D_4$}'};
labels={'\boldmath{${D_1}$}','\boldmath{$D_2$}','\boldmath{$D_3$}','\boldmath{$D_4$}'};
figure(17)
clf
[Scatter,SAx,BigAx,Histo,HAx] = plotmatrix_withr(Param(15:nparam,burnin:100:end)');

for i = 1:length(labels)
    SAx(i,1).YLabel.String=string(labels(i));
    SAx(i,1).YLabel.Rotation=0;
    SAx(i,1).YLabel.FontSize=14;
    SAx(i,1).YLabel.Units='Normalized';
    SAx(i,1).YLabel.Position=[-0.25 0.5000 0];
    for j = 1:length(labels)
%         SAx(i,j).XLim = [0 0.25];
%         SAx(i,j).YLim = [0 0.25]; 
        
        SAx(end,j).XLabel.String = string(labels(j));
        SAx(end,j).XLabel.FontSize=14;
    end
%     HAx(i).XLim = [0 0.25];
%     Histo(i).BinWidth = 0.01;
end
set(gca,'FontSize',15);
print('../figures/depthCorrelation','-dpng','-r1000')

%% Everything correlation
set(0,'defaulttextinterpreter','latex')
% Write out to file 
%acc_depths = {'$< 200 m$', '$< 400 m$', '$< 600 m$', '$< 800 m$ ', '$< 1000 m$ ', '$< 1200 m$', '$< 1400 m$', '$< 1600 m$', '$< 1800 m$', '$> 1800 m$'}';

tmp2 = Param;
tmp2(19,:) = S;

figure(18)
clf
[Scatter,SAx,BigAx,Histo,HAx] = plotmatrix_withr(tmp2(1:end,burnin:1000:end-1)');

% for i = 2:naccum+1
%     SAx(i,1).YLabel.String=string(acc_depths(i));
%     SAx(i,1).YLabel.Rotation=0;
%     SAx(i,1).YLabel.Units='Normalized';
%     SAx(i,1).YLabel.Position=[-0.5 0.5000 0];
%     for j = 2:naccum+1
%         SAx(i,j).XLim = [0 0.25];
%         SAx(i,j).YLim = [0 0.25]; 
%         
%         SAx(end,j).XLabel.String = string(acc_depths(j));
%         SAx(end,j).XLabel.FontWeight = 'bold';
%     end
% %     HAx(i).XLim = [0 0.25];
% %     Histo(i).BinWidth = 0.01;
% end

%% ice flow correlation
set(0,'defaulttextinterpreter','latex')
% Write out to file 
%acc_depths = {'$< 200 m$', '$< 400 m$', '$< 600 m$', '$< 800 m$ ', '$< 1000 m$ ', '$< 1200 m$', '$< 1400 m$', '$< 1600 m$', '$< 1800 m$', '$> 1800 m$'}';

tmp3 = Param([1,12,13,14],:);
tmp3(5,:) = S;

figure(19)
clf
[Scatter,SAx,BigAx,Histo,HAx] = plotmatrix_withr(tmp3(1:end,burnin:100:end-1)');



%% Quantify uncertainty in the difference in ages between two layers

% Take the difference between each layer for each iteration
L = nan(lp-1,numsteps+burnin);
for i = 1:lp-1
    L(i,:) = Param(nparam-lp+i+1,:) - Param(nparam-lp+i,:);
end


% Plot the distribution of each layer (any smaller?)
figure(50)
clf
set(0,'defaulttextinterpreter','latex')
colors=distinguishable_colors(length(L(:,1)));
colors(3,:) = [0 0.5 0]; %make the green darker so it's not blinding
    for i=1:lp-1
        % set plot characteristics so that each pik is a different color
        f=histogram(L(i,burnin:end)); hold on
        f.FaceColor=colors(i,:);hold on
        f.EdgeColor=colors(i,:);hold on
        %f.FaceAlpha=0.05; hold on
        f.EdgeAlpha=0.25; hold on
        %f.BinWidth=0.01; hold on
        
        % Label mean and stddev for each histo so they hopefully don't overlap
        if rem(i,2)
            text(mean(L(i,burnin:end)),max(f.Values)+1000*0.5,strcat(num2str(round(mean(L(i,burnin:end)),2)),...
            ' $\pm~$ ',num2str(round(2*std(L(i,burnin:end)),2)),' m'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        else
            text(mean(L(i,burnin:end)),max(f.Values)+1000*.75,strcat(num2str(round(mean(L(i,burnin:end)),2)),...
            ' $\pm~$ ',num2str(round(2*std(L(i,burnin:end)),2)),' m'),'Fontsize',15,'Color',colors(i,:),'FontWeight','bold')
        end
         Nmax(i,1)=max(f.Values);
    end

%% Compute mean and standard deviation of all params

for i = 1:nparam
    param_stats(i,1) = median(Param(i,:));
    param_stats(i,2) = std(Param(i,:));
end

param_stats(nparam+1,1) = mean(S);
param_stats(nparam+1,2) = std(S);
param_stats(nparam+2,1) = mean(reg);
param_stats(nparam+2,2) = std(reg);

% param_stats


%% Plotting parameter histograms
param_S = Param(11:-1:2,:); %exclude depths
param_S(11,:) = Param(1,:); %exclude depths
param_S(12:14,:) = Param(12:14,:); %exclude depths
param_S(15,:) = S;
figure(40)
clf
set(gca,'FontSize',18)
set(gca,'XTickLabel','')
set(0,'defaulttextinterpreter','latex')
for i=1:length(param_S(:,1))
    subplot(3,5,i)
    hist(gca,param_S(i,:),30); hold on
    set(gca,'FontSize',15)
    if i ==1
        xlabel('\boldmath{$\dot{a}_{0-200}$ (m/a)}')
    elseif i == 2
        xlabel('\boldmath{$\dot{a}_{200-400}$ (m/a)}')
    elseif i == 3
        xlabel('\boldmath{$\dot{a}_{400-600}$ (m/a)}')
    elseif i == 4
        xlabel('\boldmath{$\dot{a}_{600-800}$ (m/a)}')
    elseif i == 5
        xlabel('\boldmath{$\dot{a}_{800-1000}$ (m/a)}')
    elseif i == 6
        xlabel('\boldmath{$\dot{a}_{1000-1200}$ (m/a)}')
    elseif i == 7
        xlabel('\boldmath{$\dot{a}_{1200-1400}$ (m/a)}')
    elseif i == 8
        xlabel('\boldmath{$\dot{a}_{1400-1600}$ (m/a)}')
    elseif i == 9
        xlabel('\boldmath{$\dot{a}_{1600-1800}$ (m/a)}')
    elseif i == 10
        xlabel('\boldmath{$\dot{a}_{1800-2000}$ (m/a)}')
    elseif i == 11
        xlabel('\boldmath{$q$}')
    elseif i == 12
        xlabel('\boldmath{$h$}')
    elseif i == 13
        xlabel('\boldmath{$v_{ice}$}')
    elseif i == 14
        xlabel('\boldmath{$\epsilon_{firn}$}')
    elseif i == 15
        xlabel('\boldmath{$S$}')       
    end
    
    if i == 1 ||i == 6 || i == 11
        ylabel('N')
    end
%     
%     if i < 11
%         xlim([0.06, 0.16])
%     end
end
savefig('../figures/paramHist.fig')
print('../figures/paramHist','-dpng','-r1000')
%%
    
toc;