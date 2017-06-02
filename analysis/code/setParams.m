function [paramRange,nparam] = setParams(accumFlag,H)
% DESCRIPTION:
% This function sets the number and meaning of ice flow parameters computed 
% using the Metropolis algorithm
%
% INPUT:
% - *core* is ice core to derive chronology for (Byrd or WD)
% - *accumFlag* specifies the accumulation function to use (defined below)
% - *H* is length of the core (m)
%
% OUTPUT:
% - *paramRange* defines the min and max allowed parameter values
% - *nparam* is the number of parameters to use
%
% NOTE:
% - *s* is the ratio of surface to basal velocity
% - For a flow divide, expect *h*~0.7 (i.e. higher) -- see Morse et al 2012
%   This indicates the stress will go linearly to zero somewhere in the bottom half of the ice sheet.
%   Cuffey and Patterson suggest 0.3 (from the bottom of the ice sheet).
%   A conversation with Lucas suggested no more than 0.4 (from bottom of ice sheet).
%   Schwander et al. cite Patterson 1994 says .33 to 0.5 (from bottom of ice sheet.
% - WAIS Divide project website indicates current accum. rate is 22 cm/yr
% - Byrd modern accumulation rate has been measured at 11 cm/yr
% - Using Morse profile as a rule of thumb, this might have gone down as much
%   as 40% during LGM
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

%% Set number of parameters according to accumulation function in play
% Schwander et al. ice flow model is used unless noted
    
    if accumFlag     == 0
        nparam = H+2;               % accum profile, Nye h
    elseif accumFlag == 1
        nparam = 2;                 % s, constant accum (prescribed Nye h)
    elseif accumFlag == 2  
        nparam = 2;                 % s, Morse accum
    elseif accumFlag == 3 
        nparam = 3;                 % s, constant accum, h
    elseif accumFlag == 4           % Morland ice flow model
        nparam = 4;                 % 4 accum. zones (prescribed h, s)
    elseif accumFlag == 5
        nparam = 5;                 % s, 4 accumulation zones                                         
    elseif accumFlag == 6
        nparam = 6;                 % s, 4 accumulation zones, h
    elseif accumFlag == 7 
        nparam = 7;                 % s, 4 accum. zones, h, obs age uncertainty
    elseif accumFlag == 100 
        nparam = 2+50/5;            % s, h, accum zone every 5k years
    else
        disp('You chose an invalid accumulation flag.')
        disp('Next time try 1, 2 or 3. But not so much with 3.')
    end

%% Bounded uniform priors for all uncertain parameters (except uncertainty on age)
    if accumFlag == 0
        paramRange(1,1:2)=[0.0 1.7];            % s
        paramRange(2,1:2)=[0.2 0.5];            % Nye h
        paramRange(3:H+2,1:2)=repmat([0.001 0.25],[H 1]); % independent accum at each meter depth       
    elseif accumFlag == 1 
        paramRange(1,1:2)=[0.0 1.7];            % s
        paramRange(2,1:2)=[0.05 0.25];          % constant accum
    elseif accumFlag == 2                 
        paramRange(1,1:2)=[0.0 1.7];            % s
        paramRange(2,1:2)=[0.2 0.5];            % h
    elseif accumFlag == 3                       
        paramRange(1,1:2)=[0.0 1.7];            % s
        paramRange(2,1:2)=[0.05 0.25];          % constant accum
        paramRange(3,1:2)=[0.2 0.5];            % h
    elseif accumFlag == 4
        paramRange(1,1:2)=[0.01 0.25];          % acccum 1294m< depth< 2191m
        paramRange(2,1:2)=[0.01 0.25];          %accum 1024m < depth < 1294 m
        paramRange(3,1:2)=[0.05 0.25];          %accum 150m<depth<1024m
        paramRange(4,1:2)=[0.05 0.25];          %accumulation depth < 150m
    elseif accumFlag == 5 
        paramRange(1,1:2)=[0.0 1.7];            % s                                       
        paramRange(2,1:2)=[0.03 0.1];           % acccum 1294m< depth< 2191m
        paramRange(3,1:2)=[0.03 0.15];          % accum 1024m < depth < 1294 m
        paramRange(4,1:2)=[0.05 0.15];          % accum 150m<depth<1024m
        paramRange(5,1:2)=[0.05 0.15];          % accumulation depth < 150m      
    elseif accumFlag == 6 
        paramRange(1,1:2)=[0.0 1.7];            % s
        paramRange(2,1:2)=[0.01 0.25];          % acccum 1294m< depth< 2191m
        paramRange(3,1:2)=[0.01 0.25];          % accum 1024m < depth < 1294 m
        paramRange(4,1:2)=[0.05 0.25];          % accum 150m<depth<1024m
        paramRange(5,1:2)=[0.05 0.25];          % accumulation depth < 150m
        paramRange(6,1:2)=[0.1 0.9];            % Nye h
    elseif accumFlag == 7 
        paramRange(1,1:2)=[0.0 1.7];            % s
        paramRange(2,1:2)=[0.01 0.1];           % acccum 1294m< depth< 2191m
        paramRange(3,1:2)=[0.03 0.1];           % accum 1024m < depth < 1294 m
        paramRange(4,1:2)=[0.07 0.15];          % accum 150m<depth<1024m
        paramRange(5,1:2)=[0.1 0.15];          % accumulation depth < 150m
        paramRange(6,1:2)=[0.1 0.5];            % Nye h 
        paramRange(7,1:2)=[0.03 0.03];          % obs age uncertainty; this is the mu and sigma for normal dist       
    elseif accumFlag == 100
        paramRange(1,1:2) = [0.0 1.7];          % s
        paramRange(2,1:2) = [0.1 0.5];          % h
        paramRange(3:nparam,1:2) = repmat([0.001 0.25],[10 1]); % accum. zone every 5kyr
    end
end

