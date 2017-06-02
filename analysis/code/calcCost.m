function [cost] = calcLogLikelihood(H,D,z,obsAge,age)

%=========================================================================
% CALCULATE COST FOR EITHER MODEL
%=========================================================================
    
    % D is the vector of observed ages.
    DD=H-D; % make sure 0 depth is same for volcanic and modeled data

    % Index convenient for comparing obs ages and modeled ages. 
    % Model evaluated at 1m intervals, so round volcanic depths to compare
    ind=intersect(z,round(DD));

    sigAge = 0.03*obsAge;
    varAge=sigAge.^2;
    cost = -(sum((age(ind)-obsAge).^2)/(2*varAge));
end

