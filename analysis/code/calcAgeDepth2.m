%function [Ar, Param]=calcAgeDepth2(accumFlag,datFlag,core,plotting)
tic;
clear
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
        
        % Age of radar horizo1ns
        for j = 1:length(pik)
            Ar(j,i) = byrdModels2(Param(:,i),z,H,accumFlag,lp,Zr(j,i));
        end
    end
    
    
%%  The rest is plotting to test results

%plotAgeDepthHisto(Param(nparam-lp+1:nparam,:),pikDepthStats,core,pikAge,pikAgeStats,burnin,datFlag)  
plotConvergence(Param,core,burnin,accumFlag,datFlag,paramRange)    
%end
% %%
% figure(1)
% clf
% lp = length(pik);
% subplot(2,1,1)
% histogram(Param(nparam-lp+1,:),'BinWidth',10,'FaceColor','k'); hold on
% histogram(Param(nparam-lp+2,:),'BinWidth',10,'FaceColor','b'); hold on
% histogram(Param(nparam-lp+3,:),'BinWidth',10,'FaceColor','g'); hold on
% histogram(Param(nparam-lp+4,:),'BinWidth',10,'FaceColor','r'); hold on
% histogram(Param(nparam-lp+5,:),'BinWidth',10,'FaceColor','m'); hold on
% xlim([0 H])
% xlabel('Depth (m)')
% 
% subplot(2,1,2)
% histogram(Ar(1,:)/1000,'BinWidth',0.25,'FaceColor','k'); hold on
% histogram(Ar(2,:)/1000,'Binwidth',0.25,'FaceColor','b'); hold on
% histogram(Ar(3,:)/1000,'BinWidth',0.25,'FaceColor','g'); hold on
% histogram(Ar(4,:)/1000,'BinWidth',0.25,'FaceColor','r'); hold on
% histogram(Ar(5,:)/1000,'BinWidth',0.25,'FaceColor','m'); hold on
% xlim([0 100])
% xlabel('Age (ka)')

%%
figure(11)
clf
plotAgeDepthHisto2(Param(nparam-lp+1:nparam,:),core,Ar,burnin,datFlag)
  %%  
figure(10)
clf
plotSpaghettiEnvelope(age,obsAge1950,D,H,accumFlag,z,burnin,Param,datFlag,1);

% %% 
% figure(3)
% clf
% % for i = numsteps+burnin-100:numsteps+burnin
% %     plot(age(:,i)/1000,flipud(z)); hold on
% % end
% plot(mean(age(:,end-1500:end),2)/1000,flipud(z)); hold on
% plot(obsAge1950/1000,D,'k.','LineWidth',2,'MarkerSize',14);
% %plot(mean(age,2)/1000-std(age,1,2)/1000,flipud(z));
% xlim([0 100])
% ylim([0 H])
% set(gca,'YDir','reverse');
% %%
% figure(4)
% clf
% % subplot(2,1,1)
% plot(S)
% ylabel('S')
% xlabel('N_{iter}')
% 
% subplot(2,1,2)
% plot(sqrt(1./S))
% xlabel('\sigma_{age} = sqrt(1/S)')
%% Regularization
figure(12)
plot(reg(burnin:end));
xlabel('iteration')
ylabel('regularization parameter')

%% Cost
figure(13)
clf
subplot(2,1,1)
plot(cost(burnin:end,1),'r');hold on % TWTT
ylabel('TWTT cost')
subplot(2,1,2)
plot(cost(burnin:end,2),'b');hold on  % Age
xlabel('iteration')
ylabel('Age cost')

%% Accumulation
% acc_depths = [100,200,300,400,500,600,700,800,900,950,1000,1050,1100,1150,...
%               1200,1250,1300,1350,1400,1450,1500,1550,1575,1600,1625,1650,...
%               1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1950,2000,2050];
%           
acc_depths = fliplr([ 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]);
figure(14)
clf



for i = burnin:5000:burnin+numsteps
    %subplot(2,1,1)
     plot(acc_depths,Param(2:nparam-lp-3,i)); hold on
     smoothed = smooth(Param(2:nparam-lp-3,i),3);
     plot(acc_depths,smoothed, '-.','LineWidth',5); hold on
     
    % plot(smooth(Param(2:length(Param(:,1)-lp-3)),3),'lineWidth',11); hold on
%     subplot(2,1,2)
%     plot(i,var(Param(2:nparam-lp-3,:),0,1)'); hold on
end
% plot(acc_depths,mode(Param(2:length(Param(:,1))-lp-3,:),2),'k','Linewidth',10); hold on
% plot(acc_depths,mean(Param(2:length(Param(:,1))-lp-3,:),2),'b','Linewidth',10); hold on
% plot(acc_depths,median(Param(2:length(Param(:,1))-lp-3,:),2),'r','Linewidth',10); hold on

toc;
