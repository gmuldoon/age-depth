function costHisto(posterior)

%================
% Cost Histogram
%================
figure
clf
f=histogram(gca,posterior);
f.BinWidth=50;
title('Cost')
ylabel('N')
xlabel('Cost val')
text(median(posterior)-250,max(f.Values)+10,strcat(num2str(median(posterior))),'Fontsize',12,'Color','black','FontWeight','bold')
%axis([0,1e6,param_range(5,1),param_range(5,2)])
%xlim([0,3e3])
%axis 'auto y'
%axis 'auto x'

figure
clf
subplot(1,1,1)
plot(posterior)
title('Posterior')
xlabel('N_{iter}')
ylabel('Posterior value')


end

