function [p, lab] = plotSpaghetti(ageDistrib,numsteps,burnin,obsAge1950,D,H,accumFlag,z,datFlag)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

totsteps = numsteps+burnin;
%figure
%clf
ColorSet=varycolor(length(ageDistrib(1,:)));

%% Plot data over top
lab = sprintf('Observed %s ages',datFlag);


% for i=burnin:1000:totsteps
%     plot(ageDistrib(:,i)/1000,flipud(z),'Color',ColorSet(i,:)); hold on
% end
p = plot(obsAge1950/1000,D,'k.','LineWidth',2,'MarkerSize',14); hold on
set(gca,'YDir','reverse');
xlabel('Age (ka)','Fontsize',14)
ylabel('Depth (m)','Fontsize',14)
xt = get(gca, 'XTick');
yt=get(gca,'YTick');
set(gca, 'FontSize', 16)
%legend({'simulated age','volcanic data'},'Fontsize',16)
axis([0 65 0 H])
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

