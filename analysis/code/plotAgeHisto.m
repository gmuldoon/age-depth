function plotAgeDepthHisto(pikDepthUnc,pikDepthStats,core,pikAge,pikAgeStats)

colors=distinguishable_colors(length(pikAge(:,1)));
figure
clf
subplot(2,1,1)
for i=1:lp
    f=histogram(gca,(pikDepthUnc(i,:)));hold on
    f.FaceColor=colors(i,:);hold on
    f.EdgeColor=colors(i,:);hold on
    f.BinWidth=2;hold on
    f.FaceAlpha=0.25; hold on
    f.EdgeAlpha=0.5; hold on
    text(pikDepthStats(i,2),max(f.Values)+randi([-10,10])*0.5,strcat(num2str(round(pikDepthStats(i,2),1)),' \pm ',...
        num2str(round(pikDepthStats(i,4),1)),' m'),'Fontsize',12,'Color',colors(i,:),'FontWeight','bold')
    Nmax(i)=max(f.Values);
end
ax1=set(gca);
%set(gca,'view',[90 -90]) %rotate appearance of axes
if strcmp(core,'byrd')
    title('Depth Distrib for selected Byrd layers')
elseif strcmp(core,'wd')
    title('Depth Distrib for selected WD layers')
else
    title('Not sure what this is')
end
xlabel('Depth (m)')
ylabel('N')
%xlim([0 H])
X1.YLim=[0 max(Nmax)+40];

subplot(2,1,2)
for i=1:lp 
    f=histogram(pikAge(i,:)/1000); hold on
    f.FaceColor=colors(i,:);hold on
    f.EdgeColor=colors(i,:);hold on
    f.FaceAlpha=0.25; hold on
    f.EdgeAlpha=0.5; hold on
    f.BinWidth=0.1; hold on
    text(pikAgeStats(i,3)/1000-0.05,max(f.Values)+randi([-10,10])*5,strcat(num2str(round(pikAgeStats(i,2)/1000,2)),...
        ' \pm ',num2str(round(pikAgeStats(i,4)/1000,2)),' ka'),'Fontsize',12,'Color',colors(i,:),'FontWeight','bold')
    Nmax(i)=max(f.Values);
end 
%set(gca, 'Xdir', 'reverse')
%set(gca,'view',[90 -90]) %rotate appearance of axes
if strcmp(core,'byrd')
    title('Age Distrib for selected Byrd layers')
elseif strcmp(core,'wd')
    title('Age Distrib for selected WD layers')
else
    title('Not sure what this is')
end
ax2=set(gca);

xlabel('Age (ka)')
ylabel('N')
ylim([0 max(Nmax)+75]);hold on
%xlim([0 max(de_ages(lp,:)/1000)]);hold on
%xlim([0 65])
%xlim([0 age_sd(lp,2)/1000+5*age_sd(lp,4)/1000]);hold on


end

