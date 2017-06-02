function propParam = proposeParams2(nparam,paramRange,param,lp,H)
% DESCRIPTION:
% This function computes proposed parameter values 
% 
% INPUT:
% - *nparam* is the number of parameters to use
% - *paramRange* defines the min and max allowed parameter values
% - *param* is a initialized array of parameter (size nparam)
% 
% OUTPUT:
% - *propParam* is the new set of parameter values (size nparam)
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

%% Propose values for parameters as an offset from previous param value
    dp = 0.01; % proposal step size
    selecting=true;
    delta=(paramRange(:,2)-paramRange(:,1));
    %delta=(paramRange(1:nparam-lp,2)-paramRange(1:nparam-lp,1));
    propParam = nan(nparam,1);
    %disp(param(nparam-lp:nparam));
    while(selecting)
        % select Schwander params in proper range
        for i = 1:nparam-lp
            tic;
%             if i == 1 || i == 11 || i == 12 || i == 13
                selecting1 = true;
                while(selecting1)
                    propParam(i,1)=param(i,1)+dp*(0.5-rand(1,1)).*delta(i,1); %Schwander params

                   if abs(propParam(i,1)-paramRange(i,1)) <= delta(i,1) && ...
                            abs(propParam(i,1)-paramRange(i,2)) <= delta(i,1)
                        selecting1=false;
                    end 
                    t=toc;
                    if t > 15
                        disp('Took too long to find flow params')
                    end
                end
%             elseif i <= 10 && i >= 2
%                 selecting3 = true;
%                 while (selecting3)
%                     propParam(i,1) = normrnd(param(i,1), paramRange(i,2)); %sample from a normal distribution.
%                     if propParam(i,1) > 0
%                         selecting3 = false;
%                     end
%                 end
%             end
        end
        %select depths that are increasing, in proper range
        for i = 1:lp
%             propParam(nparam-lp+i,1) = param(nparam-lp+i)+dp*(0.5-rand(1,1))*delta(nparam-lp+i,1);
%              if i == 1              
%                 %choose shallowest depth like normal
%                 selecting3 = true;
%                 while(selecting3)
% %                     delta(nparam-lp+i) = paramn(nparam-lp+i+1) - 0; %delta is between surf (0) and layer 2 of last estimate
%                     %propParam(nparam-lp+i) = param(nparam-lp+i)+dp*(0.5-rand(1,1)).*delta(nparam-lp+i,1);
%                     %make sure this first horizon has a valid depth
%                     if (propParam(nparam-lp+i,1) <= paramRange(nparam-lp+i,2)) && ...
%                             (propParam(nparam-lp+i,1) >= paramRange(nparam-lp+i,1))
%                         selecting3=false;
%                     end
%                 end
%              else
                 %keep iterating until all depth values are increasing
                 selecting2 = true;
                 tic;
                 while(selecting2)
                    
                    % change possible range on lower Dr's to be deeper than those above
                    % leave upper bound the same
                    %delta(nparam-lp+i,1) = paramRange(nparam-lp+i,2) - propParam(nparam-lp+i-1);
                    %delta defined by previous accepted params
%                     if i == 1
%                         delta(nparam-lp+i) = param(nparam-lp+i+1) - 0; %delta is between surf (0) and layer 2 of last estimate
%                     elseif i < lp
%                         delta(nparam-lp+i,1) = param(nparam-lp+i+1) - param(nparam-lp+i-1);
%                     elseif i == lp
%                         delta(nparam-lp+i,1) = H - param(nparam-lp+i-1); %delta is between bottom and layer i - 1
%                     end
                    propParam(nparam-lp+i,1) = param(nparam-lp+i)+dp*(0.5-rand(1,1))*delta(nparam-lp+i,1);
                    %require each Dr to be deeper than the previous
                    %also require depths be shallower than ice thicknes
                    if i == 1
                        if (propParam(nparam-lp+i,1) <= paramRange(nparam-lp+i,2)) && ...
                                (propParam(nparam-lp+i,1) >= paramRange(nparam-lp+i,1))
                            selecting2 = false;
                            propParam(nparam-lp+1:nparam,1);
%                         else
%                             propParam(nparam-lp+1:nparam,1);
                        end
                    else
                        if (propParam(nparam-lp+i) > propParam(nparam-lp+i-1)) && ...
                                (propParam(nparam-lp+i,1) <= paramRange(nparam-lp+i,2)) && ...
                                (propParam(nparam-lp+i,1) >= paramRange(nparam-lp+i,1))
                            selecting2 = false;
                            propParam(nparam-lp+1:nparam,1);
%                         else
%                             propParam(nparam-lp+1:nparam,1);
                        end
                    end
                 end
                 t=toc;
                 if t > 15
                      disp('Took too long to find depths')
                 end 
%               end
        end
        if ~selecting1 && ~selecting2
            selecting = false;
        else
            disp('Danger if you''re here without a solution yet')
        end
    end
%     propParam
end

