function [pikDepthUnc,lp,nsamp]=twtt2depth(pik,H,core,plotting,dFirn,z)
% DESCRIPTION:
% Transform TWTT to depth. Sample vIce 
%
% INPUT:
% - *pik* is the TWTT at the core to several interpreted horizons
% - *H* is length of the core (m)
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *plotting* is a flag for whether or not to output plots here
% - *dFirn* is depth of the firn layer
% - *z* elevation above the bedrock
%
% OUTPUT:
% - *pikDepthUnc* is distribution of layer depths for each sample with unc in vIce (and rho for WD)
% - *lp* is the number of horizons with TWTT
% - *nsamp* is number of samples of vIce
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

disp('Computing the depth of the ice column with uncertainty')
addpath('../../../../Matlab_functions/distinguishable_colors/');

%% Initialize useful variables
nsamp = 1000;                            % number of samples to get for v_ice and density profile for WD
lp = length(pik);                        % number of layers
zfUnc = nan(H,nsamp);
depthUnc = nan(H,nsamp);                 % distribution of firn corrected depths for the whole ice column with unc in v_ice (and rho for WD)
pikDepthUnc = nan(lp,nsamp);             % distribution of layer depths for each sample with unc in v_ice (and rho for WD)
pikDepthStd = nan(lp,1);                 % standard dev of depth to each pik (incorporates all three errors)
zfStd = nan(H,1);                        % standard deviation of firn correction
pwErrUnc = nan(lp,nsamp);                % distribution of pw errors for each layer
pwErrStd = nan(lp,1);                    % standard dev of pulse wifth error to each layer

%% Sample depth error
for n=1:nsamp
    [zf,depth,pikDepth,vIce,pwErr] = radarDepth(pik,core,H); % each iteration script randomly samples v_ice for both ice cores, also samples rho for WD, phase width sampling of the surface and layers
    zfUnc(:,n) = zf;
    depthUnc(:,n) = depth;
    pikDepthUnc(:,n) = pikDepth;
    vIce(n,1) = vIce;
    pwErrUnc(:,n) = pwErr;
end

for i=1:lp
    pikDepthStd(i) = std(pikDepthUnc(i,:)); % uncertainty in depth from vIce, density, phase width sampling
    pwErrStd(i) = std(pwErrUnc(i,:));       % uncertainty from only the pulse width
end

for i=1:H
    zfStd(i) = std(zfUnc(i,:));             % uncertainty from firn correction as a function of depth
end

%% Make figures related to the depth calculation
colors=distinguishable_colors(lp);
if plotting
    equal={' =  '};
    space={'  '};
% Plot histograms of layer depths with unc from ice velocity, firn correction
    figure
    subplot(1,1,1)
    clf
    Nmax=nan(lp,1);
    for i=1:lp
        h=histogram(gca,(pikDepthUnc(i,:)));hold on
        h.FaceColor=colors(i,:);hold on
        h.EdgeColor=colors(i,:); hold on
        h.EdgeAlpha=0.5; hold on
        h.FaceAlpha=0.15; hold on
        text(mean(pikDepthUnc(i,:)),max(h.Values)+2,strcat('\sigma_{depth}',equal,num2str(round(pikDepthStd(i),2)),...
            'm,',space,'\sigma_{PW}',equal,num2str(round(pwErrStd(i),2)),'m'),'Fontsize',12,'Color',colors(i,:),'FontWeight','bold')
        Nmax(i)=max(h.Values);
    end
    set(gca, 'Xdir', 'reverse')
    set(gca,'view',[90 -90]) %rotate appearance of axes
    if strcmp(core,'Byrd')
        title('Distrib of layer depths from unc in ice velocity and density for Byrd')
    elseif strcmp(core,'WD')
        title('Distrib of layer depths from unc in ice velocity and density for WD')
    else
        title('Not sure what this is')
    end
    xlabel('Depth (m)')
    ylabel('Frequency')
    set(gca, 'FontSize', 14)
    %xmax=max(pikDepthUnc(lp,:));
    xlim([0 H])
    ymax=max(Nmax(i))+80;
    ylim([0 ymax])
    text(dFirn/2,ymax/2,strcat('\sigma_{firn}',equal,num2str(round(zfStd(H),2)),'m'),'Fontsize',12,'Color','k','FontWeight','bold')
    patch([0 dFirn dFirn 0], [0 0 max(ylim)*[1 1]], [0.8 0.8 0.8],'FaceAlpha',0.3,'EdgeAlpha',0)
    plot([dFirn dFirn],ylim,'k--','LineWidth',2)
    %set(gca,'view',[0 0])
    
    print('../figures/firncorrection','-dpng')
    
% Plot ensemble of firn corrections and std of distrib as a func of depth (only applied to WD)
    if strcmp(core,'WD')
        figure
        clf
        for i=1:50:nsamp
            line(zfUnc(:,i),z,'Color','r'); hold on
        end
            
        text(min(zfUnc(end,:))-1,dFirn+3,strcat('\mu_{firn}',equal,num2str(round(mean(zfUnc(end,:)),2,'significant')),' m'),...
            'Color','r','FontSize',14,'Fontweight','bold')
        text(min(zfUnc(end,:))-1,dFirn+7,strcat('\sigma_{firn}',equal,num2str(round(zfStd(H),2,'significant')),' m'),...
            'Color','k','FontSize',14,'Fontweight','bold')
        text(0.5,dFirn-2, strcat('Firn layer'),'Color','k','FontSize',14,'Fontweight','bold')
        text(0.5,dFirn+2, strcat('Ice'),'Color','k','FontSize',14,'Fontweight','bold')
        ax1 = gca; % current axes
        ax1.XColor = 'r';
        ax1.YColor = 'r';
        ax1.YLim=([0 dFirn+20]); 
        %xmax=mean(zf_unc(end,:))*1.1;
        %ax1.XLim=([0 xmax]);
        ax1.YDir='reverse';
        ax1.FontSize=14;
        ax1_pos = ax1.Position; % position of first axes
        ylabel('depth (m)','Fontsize',14)
        xlabel('firn correction distribution (m)','Fontsize',14)
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'YDir','reverse',...
            'YLim',([0 dFirn+20]),...
            'FontSize',14);
        line(zfStd,z,'Parent',ax2,'Color','k'); hold on
        xlabel('std of firn correction distribution (m)','Fontsize',14)
        ylabel('depth (m)','Fontsize',14)
        set(gca, 'FontSize', 14)
        plot(ax2.XLim,[dFirn dFirn],'k--','LineWidth',3)
        patch([0 max(ax2.XLim) max(ax2.XLim) 0], [0 0 dFirn dFirn], [0.8 0.8 0.8],'FaceAlpha',0.3,'EdgeAlpha',0)         
    end
end
