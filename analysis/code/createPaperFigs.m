%% createPaperFigs.m
%This script takes the flag and core as input and creates figures for the
%layers paper. This is necessary because the agedepth_main.m script runs
%the analysis for one core and one flag at a time. I'm interested in
%creating figures with multiple cores and flags included.

clear
addpath('../Matlab_functions');

global figCount colors
figCount=0;
%==========================================================================
%Select flag based on the type of accum inversion to do
% 0: Byrd ice core, flag will be automatically set to 0 if byrd also selected
%    At Byrd, layer thicknesses are used as proxy
% 1: invert for a constant accumulation rate with depth
% 2: assume an accumulation rate varying with depth according to the Morse
% et al 2002 functional form
% 3: invert for a depth-varying accumulation rate (Not YET IMPLEMENTED)
%==========================================================================
flag = [1;2;3;4];

core = {'byrd'; 'wd'};

%% Call the function to get output for the desired flags and cores.
for i = 1:size(core)
    Core = core(i);
    if strcmp(Core,'wd')
       %for j = 1:size(flag)
       for j = 1:1
           Flag = flag(j);
            [z,D,pikDepthUnc,depthStd,obsAge1950,fullageStd,ageStd,deAges,H]=ageDepthMain(Flag,Core);
%            if Flag == 1
%                disp('WD1')
%                WD1 = {z;D;pik_depth_unc;depth_sd;obs_age1950;fullage_sd;age_sd;de_ages;H};
%            elseif Flag == 2
%                disp('WD2')
%                WD2 = {z;D;pik_depth_unc;depth_sd;obs_age1950;fullage_sd;age_sd;de_ages;H};
%            elseif Flag == 3
%                disp('WD3')
%                WD3 = {z;D;pik_depth_unc;depth_sd;obs_age1950;fullage_sd;age_sd;de_ages;H};
%             elseif Flag == 4
%                disp('WD4')
%                WD4 = {z;D;pik_depth_unc;depth_sd;obs_age1950;fullage_sd;age_sd;de_ages;H};
%            end
        end
         WD = {z;D;pikDepthUnc;depthStd;obsAge1950;fullageStd;ageStd;deAges;H};

    elseif strcmp(Core,'byrd')
        Flag = 0;
        [z,D,pikDepthUnc,depthStd,obsAge1950,fullageStd,ageStd,deAges,H]=ageDepthMain(Flag,Core);
        Byrd = {z;D;pikDepthUnc;depthStd;obsAge1950;fullageStd;ageStd;deAges;H};
    end
end

%% Make plots as desired


% % Age-depth profiles of Byrd core with different accumulation functions
%figCount = figCount+1;
% figure(figCount)
% clf
% plot(WD1{1},WD1{6}(:,3)/1000,'k'); hold on
% plot(WD2{1},WD2{6}(:,3)/1000,'r'); hold on
% plot(WD3{1},WD3{6}(:,3)/1000,'b'); hold on
% plot(WD4{1},WD4{6}(:,3)/1000,'g'); hold on
% 
% legend({'constant accum','Morse et al. accum','Byrd-like accum', 'WD obs accum'},'Fontsize',16)
% 
% shadedErrorBar(WD1{1},WD1{6}(:,3)/1000,2*WD1{6}(:,4)/1000,'k--',1); hold on
% shadedErrorBar(WD1{1},WD1{6}(:,3)/1000,WD1{6}(:,4)/1000,{'k--','markeredgecolor','k'},1);hold on
% errorbar(WD1{4}(:,3),WD1{7}(:,3)/1000,WD1{7}(:,4)*2/1000,'k.'); hold on
% herrorbar(WD1{4}(:,3),WD1{7}(:,3)/1000,WD1{4}(:,4)*2,'k.'); hold on
% 
% shadedErrorBar(WD2{1},WD2{6}(:,3)/1000,2*WD2{6}(:,4)/1000,'r--',1); hold on
% shadedErrorBar(WD2{1},WD2{6}(:,3)/1000,WD2{6}(:,4)/1000,{'r--','markeredgecolor','r'},1);hold on
% errorbar(WD2{4}(:,3),WD2{7}(:,3)/1000,WD2{7}(:,4)*2/1000,'r.'); hold on
% herrorbar(WD2{4}(:,3),WD2{7}(:,3)/1000,WD2{4}(:,4)*2,'r.'); hold on
% 
% shadedErrorBar(WD3{1},WD3{6}(:,3)/1000,2*WD3{6}(:,4)/1000,'b--',1); hold on
% shadedErrorBar(WD3{1},WD3{6}(:,3)/1000,WD3{6}(:,4)/1000,{'b-.','markeredgecolor','b'},1);hold on
% errorbar(WD3{4}(:,3),WD3{7}(:,3)/1000,WD3{7}(:,4)*2/1000,'b.'); hold on
% herrorbar(WD3{4}(:,3),WD3{7}(:,3)/1000,WD3{4}(:,4)*2,'b.'); hold on
% 
% shadedErrorBar(WD4{1},WD4{6}(:,3)/1000,2*WD4{6}(:,4)/1000,'g--',1); hold on
% shadedErrorBar(WD4{1},WD4{6}(:,3)/1000,WD4{6}(:,4)/1000,{'g-.','markeredgecolor','g'},1);hold on
% errorbar(WD4{4}(:,3),WD4{7}(:,3)/1000,WD4{7}(:,4)*2/1000,'g.'); hold on
% herrorbar(WD4{4}(:,3),WD4{7}(:,3)/1000,WD4{4}(:,4)*2,'g.'); hold on
% 
% plot(WD1{2},WD1{5}/1000,'kv','LineWidth',5); hold on
% 
% axis([0 WD1{9} 0 75 ])
% set(gca,'Xdir','reverse')
% ylabel('Age (ka)','Fontsize',14)
% xlabel('Depth (m)','Fontsize',14)
% set(gca,'view',[90 -90])

%%
% Age-depth profile of WD and Byrd ice cores (using lowest cost accum)
figCount = figCount+1;
figure(figCount)
clf

shadedErrorBar(Byrd{1},Byrd{6}(:,3)/1000,2*Byrd{6}(:,4)/1000,'r--',1); hold on
shadedErrorBar(Byrd{1},Byrd{6}(:,3)/1000,Byrd{6}(:,4)/1000,{'r-.','markeredgecolor','r'},1);hold on
errorbar(Byrd{4}(:,3),Byrd{7}(:,3)/1000,Byrd{7}(:,4)*2/1000,'r'); hold on
herrorbar(Byrd{4}(:,3),Byrd{7}(:,3)/1000,Byrd{4}(:,4)*2/1000,'r'); hold on

shadedErrorBar(WD{6}(:,1),WD{6}(:,2)/1000,WD{6}(:,4)/1000,'k',1); hold on
shadedErrorBar(WD{6}(:,1),WD{6}(:,2)/1000,WD{6}(:,4)/1000/2,'k',1); hold on

axis([0 WD{6}(end,1) 0 75 ])
set(gca,'Xdir','reverse')
ylabel('Age (ka)','Fontsize',14)
xlabel('Depth (m)','Fontsize',14)
set(gca,'view',[90 -90])

 %%
% % Age and depth histos for both cores
% %need: for byrd and wd, flag = 3
% %pik_depth_unc, depth_sd, de_ages,age_sd
% figCount = figCount+1;
% figure(figCount)
% clf
% subplot(2,1,1)
% for i=1:lp
%     f=histogram(gca,(pik_depth_unc(i,:)));hold on
%     f.FaceColor=colors(i,:);hold on
%     f.EdgeColor=colors(i,:);hold on
%     f.BinWidth=2;hold on
%     f.FaceAlpha=0.25; hold on
%     f.EdgeAlpha=0.5; hold on
%     text(depth_sd(i,2),max(f.Values)+randi([-10,10])*0.5,strcat(num2str(round(depth_sd(i,2),1)),' \pm ',...
%         num2str(round(depth_sd(i,4),1)),' m'),'Fontsize',12,'Color',colors(i,:),'FontWeight','bold')
%     Nmax(i)=max(f.Values);
% end
% ax1=set(gca);
% %set(gca,'view',[90 -90]) %rotate appearance of axes
% if strcmp(core,'byrd')
%     title('Depth Distrib for selected Byrd layers')
% elseif strcmp(core,'wd')
%     title('Depth Distrib for selected WD layers')
% else
%     title('Not sure what this is')
% end
% xlabel('Depth (m)')
% ylabel('N')
% %xlim([0 H])
% X1.YLim=[0 max(Nmax)+40];
% 
% subplot(2,1,2)
% for i=1:lp 
%     f=histogram(de_ages(i,:)/1000); hold on
%     f.FaceColor=colors(i,:);hold on
%     f.EdgeColor=colors(i,:);hold on
%     f.FaceAlpha=0.25; hold on
%     f.EdgeAlpha=0.5; hold on
%     f.BinWidth=0.1; hold on
%     text(age_sd(i,3)/1000-0.05,max(f.Values)+randi([-10,10])*5,strcat(num2str(round(age_sd(i,2)/1000,2)),...
%         ' \pm ',num2str(round(age_sd(i,4)/1000,2)),' ka'),'Fontsize',12,'Color',colors(i,:),'FontWeight','bold')
%     Nmax(i)=max(f.Values);
% end 
% %set(gca, 'Xdir', 'reverse')
% %set(gca,'view',[90 -90]) %rotate appearance of axes
% if strcmp(core,'byrd')
%     title('Age Distrib for selected Byrd layers')
% elseif strcmp(core,'wd')
%     title('Age Distrib for selected WD layers')
% else
%     title('Not sure what this is')
% end
% ax2=set(gca);
% 
% xlabel('Age (ka)')
% ylabel('N')
% ylim([0 max(Nmax)+75]);hold on
% %xlim([0 max(de_ages(lp,:)/1000)]);hold on
% xlim([0 65])
% %xlim([0 age_sd(lp,2)/1000+5*age_sd(lp,4)/1000]);hold on
% 
% 

