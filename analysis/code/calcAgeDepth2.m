%function [Ar, Param]=calcAgeDepth2(accumFlag,datFlag,core,plotting)
tic;
accumFlag=7;                    % type of accum solver. (see setParams.m for options)
datFlag=string('volcanic');   % observational chronology to constrain analysis
core = 'Byrd';              % ice core to derive chronology for (WD or Byrd)
depthPlotting = 1;          % whether or not to plot depth calculation

%% 1. Set up environment
    addpath('/Users/gail/Documents/Research/Matlab_functions');
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
plotConvergence(Param,core,burnin,accumFlag,datFlag,paramRange)    

%%  Make spaghetti plot
figure(10)
clf
[~,~,pikDepthUncWD,wdAge1950,wdAge1950Unc] = plotSpaghettiEnvelope(age,obsAge1950,D,H,accumFlag,z,burnin,Param,datFlag,1,Ar,Zr,1);
print('../figures/spaghetti','-dpng')

%% Plot age-depth histograms
plotAgeDepthHisto2(Param(nparam-lp+1:nparam,:),core,Ar,burnin,datFlag,pikDepthUncWD,wdAge1950,wdAge1950Unc)

%% Regularization plot
set(0,'defaulttextinterpreter','latex')
figure(12)
plot(reg(burnin:end),'LineWidth',3);
xlabel('Metropolis iteration','FontSize', 15)
ylabel('$r$')
set(gca, 'FontSize', 15)
xlim([0 numsteps])
print('../figures/regularization','-dpng')

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
print('../figures/cost','-dpng')

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
print('../figures/accumdepth','-dpng')

%% Accumulation paper figure
%sort accums by age cost
[cost_sorted,cost_sorted_I] = sort(cost(:,2),1); 
Param_sorted = Param(:,cost_sorted_I);

acc_depths = fliplr([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]);
naccum = length(acc_depths);
AccumSortedFull = zeros(H,burnin+numsteps);
for i = 1:naccum
    if i == 1
        for j = H:-1:acc_depths(1)
            AccumSortedFull(j,:) = Param_sorted(i+1,:);
        end
    elseif i == naccum
        for j = acc_depths(i-1):-1:1
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
cm=colormap(parula((burnin+numsteps)/n));
cnt = 0;
for i = burnin:n:burnin+numsteps
    cnt = cnt+1;
    plot(smooth(AccumSortedFull(:,i),50),'color',cm(cnt,:),'LineWidth',3); hold on
    %accVar(i) = var(Param(2:nparam-lp-3,i));    
end

xlabel('Depth (m)','FontSize',15)
ylabel('$\dot{a}$ (m/a)','FontSize',15)
xlim([0 2164]);
h=colorbar;
set(get(h,'title'),'string','Normalized Cost');
set(gca,'FontSize',15);

print('../figures/accumdepthSorted','-dpng')


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
%    agediv2,obsAge1950,D,H,accumFlag,z,burnin,Paramdiv2,datFlag,1,Ardiv2,Zrdiv2,2);

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
% figure(17);
% clf
% text(0.25,0,'Depth (m)')
% text(0.75,0, 'Age (ka)')
% 
% subplot(4,2,1)
% subplot(4,2,2)
% [ax1,ax2] = plotAgeDepthHisto2_ke(Drdiv1,core,Ardiv1,burnin,1);
% 
% subplot(4,2,3)
% subplot(4,2,4)
% [ax3,ax4] = plotAgeDepthHisto2_ke(Drdiv2,core,Ardiv2,burnin,3);
% 
% subplot(4,2,5)
% subplot(4,2,6)
% [ax5,ax6] = plotAgeDepthHisto2_ke(Drdiv4,core,Ardiv4,burnin,5);
% 
% subplot(4,2,7)
% subplot(4,2,8)
% [ax7,ax8] = plotAgeDepthHisto2_ke(Drdiv5,core,Ardiv5,burnin,7);

%print('../figures/keCompare','-dpng')

%% Correlation plots
set(0,'defaulttextinterpreter','latex')
% Write out to file 
acc_depths = {'$< 200 m$', '$< 400 m$', '$< 600 m$', '$< 800 m$ ', '$< 1000 m$ ', '$< 1200 m$', '$< 1400 m$', '$< 1600 m$', '$< 1800 m$', '$> 1800 m$'}';

figure(16)
clf
[Scatter,SAx,BigAx,Histo,HAx] = plotmatrix_withr(Param(11:-1:2,:)'); % 11 is the shallowest

for i = 1:naccum
    SAx(i,1).YLabel.String=string(acc_depths(i));
    for j = 1:naccum
        SAx(i,j).XLim = [0 0.25];
        SAx(i,j).YLim = [0 0.25]; 
        
        SAx(end,j).XLabel.String = string(acc_depths(j));
        SAx(end,j).XLabel.FontWeight = 'bold';
    end
    HAx(i).XLim = [0 0.25];
    Histo(i).BinWidth = 0.01;
end

print('../figures/accumCorrelation','-dpng')
%%



toc;