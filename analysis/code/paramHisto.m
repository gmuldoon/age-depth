function paramHisto(paramDistrib,accumFlag,paramRange,core)

%===========================
% Model Parameter Histograms
%===========================


figure
clf
f=histogram(gca,(paramDistrib(1,:)));
text(median(paramDistrib(1,:))-0.005,max(f.Values)+10,strcat(num2str(median(paramDistrib(1,:)))),'Fontsize',12,'Color','black','FontWeight','bold')
    
xlabel('Ratio of bed to surface velocity at WD')
ylabel('N')


if accumFlag == 1 
    figure
    %clf
    hist(paramDistrib(3,:))
    xlabel('Constant accumulation rate in the ice core m/yr)')
    ylabel('N')
    xlim([paramRange(2,1),paramRange(2,2)])
end
if strcmp(core,'wd')
    figure
    %clf
    f=histogram(gca,(paramDistrib(2,:)));
    text(median(paramDistrib(2,:))-0.01,max(f.Values)+10,strcat(num2str(median(paramDistrib(2,:)))),'Fontsize',12,'Color','black','FontWeight','bold')
    xlabel('D-J "kink" as a fraction of depth to bed')
    ylabel('N')
    xlim([paramRange(2,1),paramRange(2,2)])
end

end

