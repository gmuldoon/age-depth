function dePik = sampleDepth(numsteps,nsamp,pikDepthUnc,lp)
% DESCRIPTION:
% Sample a vector of horizon depths for each Metropolis iteration
%
% INPUT:
% - *numsteps* is number of Metropolis samples after burnin
% - *nsamp* is number of samples of vIce
% - *pikDepthUnc* is distribution of layer depths for each sample with unc in vIce
% - *lp* is the number of horizons with TWTT
%
% OUTPUT:
% - *dePik* is a sampled vector of horizon depths
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

%% Sample a vector of horizon depths for each Metropolis iteration
% Will be used to evaluate the age-depths of horizons for each Metropolis iteration
    dePik = nan(lp,numsteps);    
    for n=1:numsteps
        samp_depth = randi([1 nsamp]);          % random sample to choose a depth profile for the piks
        dePik(:,n) = pikDepthUnc(:,samp_depth); % choose the corresponding depth profile
    end

end

