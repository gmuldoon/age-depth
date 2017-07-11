figure(13)
clf
for i = 1:numsteps
    plot(mean(Param(2:nparam-lp-3,:),2)); hold on
 end

%Use variance of the mean as the reference
rednoise = cumsum(randn(1,length(Param(1,:))));
rednoise_normal = rednoise/max(abs(rednoise))*0.1;
var(mean(Param(2:nparam-lp-3,:)+rednoise_normal,2))


plot(mean(Param(2:nparam-lp-3,:),2)); hold on
plot(mean(Param(2:nparam-lp-3,:)+rednoise_normal,2))

plot(flipud(smooth(mean(Param(2:nparam-lp-3,1:500:end),2),5)))

var(flipud(smooth(mean(Param(2:nparam-lp-3,:),2),5)))

var(flipud(mean(smooth(Param(2:nparam-lp-3,:),5),3)))

figure(1)
clf
acc_depths = [ 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000];
for i = burnin:500:numsteps
    plot(acc_depths',flipud(smooth(Param(2:nparam-lp-3,i),3)));hold on
    smoothed(:,i) = flipud(smooth(Param(2:nparam-lp-3,i),3));
end
%plot(acc_depths',mean(smoothed,2),'k','LineWidth',10);hold on