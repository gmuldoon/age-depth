function plotAgeDepthError(ageDistrib,z,obsAge1950,D,pikAgeStats,pikDepthStats,accumFlag,H,fullAgeStats,burnin)

    for n=1:H
        [f,age]=ksdensity(ageDistrib(n,burnin:end));
        [~, I] = max(f);
        bestage(n) = age(I);
    end

    addpath('../Matlab_functions/herrorbar');
    figure
    clf
    plot(bestage/1000,flipud(z),'b','LineWidth',1); hold on
    plot(fullAgeStats(:,3)/1000,flipud(z),'k'); hold on
    plot(obsAge1950/1000,D,'k.','LineWidth',5,'Markersize',15); hold on
    errorbar(pikAgeStats(:,3)/1000,pikDepthStats(:,3),pikDepthStats(:,4)*2,'b.'); hold on
    herrorbar(pikAgeStats(:,3)/1000,pikDepthStats(:,3),pikAgeStats(:,4)*2/1000,'r.'); hold on

    set(gca,'YDir','reverse');
    xlabel('Age (ka)','Fontsize',14)
    ylabel('Depth (m)','Fontsize',14)
    %xt = get(gca, 'XTick');yt=get(gca,'YTick');
    set(gca, 'FontSize', 16)
    axis([0 50 0 H])
    legend({'max ksdensity','median simulated age','volcanic data','2 sd age error','2 sd depth error'},'Fontsize',14)
    if accumFlag == 1
        title('Metropolis result: Age-depth profiles for constant accumulation inversion','Fontsize',16)
    elseif accumFlag == 2
        title('Metropolis result: Age-depth profile for accum rate parameterized by Morse et al 02 func','Fontsize',16)
    elseif accumFlag == 3
        title('Metropolis result: Age-depth profile for depth-varying accum inversion','Fontsize',16)
    elseif accumFlag == 0
        title('Metropolis result: Age-depth profile using Byrd ice core','Fontsize',16)
    end
end
