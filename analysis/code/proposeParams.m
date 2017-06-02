function propParam = proposeParams(nparam,paramRange,param)
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
    dp = 1; % proposal step size
    selecting=true;
    delta=(paramRange(:,2)-paramRange(:,1));
    
    while(selecting)
        propParam=param+dp*(0.5-rand(nparam,1)).*delta;
        propParam(7) = paramRange(7,1); % this is a placeholder
        % if valid values (within range) for every parameter, quit the loop
        if (sum(abs(propParam-paramRange(:,1)) <= delta)==nparam) && (sum(abs(propParam-paramRange(:,2)) <= delta)==nparam)
            selecting=false;
        end 
    end
%% Propose value for age uncertainty separately
    if nparam == 7 
        %propParam(7) = normrnd(paramRange(7,1),paramRange(7,2)); % sample from normal distribution
        %propParam(7) = normrnd(paramRange(7,1),paramRange(7,2)/2);
        propParam(7) = gamrnd(1,paramRange(7,1));
        %propParam(7) = 0.03;
    end
end

