function  plotConvergence(paramDistrib,core,burnin,accumFlag,datFlag,paramRange)
% DESCRIPTION:
% This function makes convergence plots of the ice flow parameters selected
% by the Metropolis sampler
%
% INPUT:
% - *paramDistrib* is a nparam-by-totsteps ensemble of parameter values from the Metropolis algorithm
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *burnin* is number of iterations to ignore from the start of the Metropolis sampling
% - *accumFlag* specifies the accumulation function to use, as defined in setParams.m
% - *datFlag* specifies the observed chronology at Byrd to use 
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 16 Feb 2017
%
%% Set up a loop index
    niter = length(paramDistrib(end,burnin:end));
    n = 1:niter';
%% Create convergence plots for each parameter
%     figure
%     clf
%     if strcmp(core,'Byrd')
%         suptitle(sprintf('Parameter convergence for %s core with %s data',core,string(datFlag)))
%     else
%         suptitle(sprintf('Parameter convergence for %s core',core))
%     end
%     if accumFlag == 6 || accumFlag == 5 || accumFlag == 7
%         subplot(2,3,1)
%         %plot(n,paramDistrib(1,burnin:end))
%         hist(paramDistrib(1,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('s param')
%         ylim([paramRange(1,1) paramRange(1,2)])
%         
% 
%         %figure(2)
%         %clf
%         subplot(2,3,2)
%         %plot(n,paramDistrib(2,burnin:end))
%         hist(paramDistrib(2,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('acccum 1294m< depth< 2191m')
%         %ylim([paramRange(2,1) paramRange(2,2)])
% 
%         %figure(3)
%         %clf
%         subplot(2,3,3)
%         %plot(n,paramDistrib(3,burnin:end))
%         hist(paramDistrib(3,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('%accum 1024m < depth < 1294 m')
%         %ylim([paramRange(3,1) paramRange(3,2)])
% 
%         %figure(4)
%         %clf
%         subplot(2,3,4)
%         %plot(n,paramDistrib(4,burnin:end))
%         hist(paramDistrib(4,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('accum 150m<depth<1024m')
%         %ylim([paramRange(4,1) paramRange(4,2)])
% 
%         %figure(5)
%         %clf
%         subplot(2,3,5)
%         %plot(n,paramDistrib(5,burnin:end))
%         hist(paramDistrib(5,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('accumulation depth < 150m')
%         %ylim([paramRange(5,1) paramRange(5,2)])

%         if accumFlag == 6 || accumFlag == 7
%             subplot(2,3,6)
%             %plot(n,(paramDistrib(6,burnin:end)))
%             hist((paramDistrib(6,burnin:end)))
%             xlabel('N_{iter}')
%             ylabel('h')
%             ylim([paramRange(6,1) paramRange(6,2)])
%         end
%         if accumFlag == 7
%             figure
%             if strcmp(core,'Byrd')
%                 suptitle(sprintf('Parameter convergence for %s core with %s data',core,string(datFlag)))
%             else
%                 suptitle(sprintf('Parameter convergence for %s core',core))
%             end
%             subplot(2,1,1)
%             %plot(n,paramDistrib(7,burnin:end))
%             hist(paramDistrib(7,burnin:end))
%             xlabel('N_{iter}')
%             ylabel(sprintf('v_{ice}'))
%             ylim([paramRange(7,1) paramRange(7,2)])
%             
%             subplot(2,1,2)
%             %plot(n,paramDistrib(8,burnin:end))
%             hist(paramDistrib(8,burnin:end))
%             xlabel('N_{iter}')
%             ylabel(sprintf('dFirn'))
%             ylim([paramRange(8,1) paramRange(8,2)])
%             
%             for i = 1:5
%                 subplot(2,3,i)
%                 %plot(n,paramDistrib(length(paramRange(:,1))-5+i,burnin:end))
%                 hist(paramDistrib(length(paramRange(:,1))-5+i,burnin:end))
%                 ylim([paramRange(length(paramRange(:,1))-5+i,1) paramRange(length(paramRange(:,1))-5+i,2)])
%                 xlabel('N_{iter}')
%                 ylabel(sprintf('horizon%i', i ))
%             end
            
            
%         end
%     elseif accumFlag == 4
%         figure
%         clf
%         subplot(2,2,1)
%         plot(n,paramDistrib(1,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('acccum 1294m< depth< 2191m')
% 
%         %figure(3)
%         %clf
%         subplot(2,2,2)
%         plot(n,paramDistrib(2,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('%accum 1024m < depth < 1294 m')
% 
%         %figure(4)
%         %clf
%         subplot(2,2,3)
%         plot(n,paramDistrib(3,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('accum 150m<depth<1024m')
% 
%         %figure(5)
%         %clf
%         subplot(2,2,4)
%         plot(n,paramDistrib(4,burnin:end))
%         xlabel('N_{iter}')
%         ylabel('accumulation depth < 150m')
% 
%     elseif accumFlag == 3
%         figure
%         clf
% 
%         subplot(2,2,1)
%         plot(n,(paramDistrib(1,burnin:end)))
%         xlabel('N_{iter}')
%         ylabel('s param')
% 
%         subplot(2,2,2)
%         plot(n,(paramDistrib(2,burnin:end)))
%         xlabel('N_{iter}')
%         ylabel('h')
% 
%         subplot(2,2,3)
%         plot(n,(paramDistrib(3,burnin:end)))
%         xlabel('N_{iter}')
%         ylabel('constant accum')
% 
%     elseif accumFlag == 2
%         figure
%         clf
%         subplot(2,1,1)
%         plot(n,(paramDistrib(1,burnin:end)))
%         xlabel('N_{iter}')
%         ylabel('s param')
% 
%         subplot(2,1,2)
%         plot(n,(paramDistrib(2,burnin:end)))
%         xlabel('N_{iter}')
%         ylabel('s param')
%     else 
%         figure
%         clf
%         fact = factor(length(paramDistrib(:,1))); %determine number of rows and cols
%         if length(fact) == 1
%             fact = factor(length(paramDistrib(:,1))+1);
%         end
%         nrows = fact(end);
%         ncols = ceil(length(paramDistrib(:,1))/nrows);
% 
%         for i = 1:length(paramDistrib(:,1))         %plot the convergence plot
%             subplot(nrows,ncols,i)
%             plot(n,(paramDistrib(i,burnin:end)))
%             xlabel('N_{iter}')
%             ylabel(sprintf('param %d',i))  
%         end
%     end

%     set(gcf,'NextPlot','add');
%     axes;
%     set(gca,'Visible','off');
%     
    nparam = length(paramDistrib(:,1));
    lp=5;
    
    figure(20)
    clf
    subplot(2,2,1)
    plot(paramDistrib(1,:))
    ylabel('s param')
    subplot(2,2,2)
    plot(paramDistrib(nparam-lp-2,:))
    ylabel('h param')
    subplot(2,2,3)
    plot(paramDistrib(nparam-lp-1,:))
    ylabel('vice')
    subplot(2,2,4)
    plot(paramDistrib(nparam-lp,:))
    ylabel('dfirn')
    
    figure(21)
    clf
    subplot(2,2,1)
    hist(paramDistrib(1,:))
    ylabel('s param')
    subplot(2,2,2)
    hist(paramDistrib(nparam-lp-2,:))
    ylabel('h param')
    subplot(2,2,3)
    hist(paramDistrib(nparam-lp-1,:))
    ylabel('vice')
    subplot(2,2,4)
    hist(paramDistrib(nparam-lp,:))
    ylabel('dfirn')
    
%     figure(22)
%     for i = burnin:length(paramDistrib(1,:))
%         plot(paramDistrib(2:nparam-lp-3,i)); hold on
%     end
    
%     figure(23)
%     n=length(paramDistrib(:,1))-lp-3;
%     plot(var(paramDistrib(2:n,:),0,1)); hold on
%     x = (1:niter+burnin-1)';
%     b = var(paramDistrib(2:n,:),0,1)';
%     fit = polyfit(x,b,1);
%     y=fit(1)*x+fit(2);
%     fit(1)
%     plot(x,y,'k');hold on
    

end

