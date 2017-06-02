function plotCostHisto(posterior)

%================
% Cost Histogram
%================
figure
clf
subplot(2,1,1)
histogram(posterior);
%f=histogram(gca,posterior);
%f.BinWidth=1e-19;
title('Cost')
ylabel('N')
xlabel('Cost val')
%text(median(posterior)-250,max(f.Values)+10,strcat(num2str(median(posterior))),'Fontsize',12,'Color','black','FontWeight','bold')
%axis([0,1e6,param_range(5,1),param_range(5,2)])
xlim([0,max(posterior)])
%axis 'auto y'
%axis 'auto x'

subplot(2,1,2)
plot(posterior)
title('Posterior')
xlabel('N_{iter}')
ylabel('Posterior value')


end

