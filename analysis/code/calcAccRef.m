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