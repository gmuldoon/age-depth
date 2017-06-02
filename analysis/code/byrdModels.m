function age = byrdModels(param,z,H,accumFlag)
% DESCRIPTION:
%  This function contains two simple ice flow models for the Byrd ice core.
% 
% INPUT:
% - *param* is ensemble of values for each parameter (size: nparam)
% - *z* is height above the bed (m)
% - *H* is length of the core (m)
% - *accumFlag* specifies the accumulation function to use, as defined in setParams.m
%
% OUTPUT:
% - *age* is the profile of simulated ages (size: H)
%
% NOTE:
% - *h* is Nye parameter
% - *hh* is _depth_ of the shear layer described by Schwander et al
% - *h*  is fraction _above_ the bed
% - See setParams.m for more info on what each accumFlag is
%
% Author: Gail Muldoon
% Email: gailrmuldoon@gmail.com
% Last revision: 15 Feb 2017

%% Initialize some useful things
    lz=length(z);
    age = nan(lz,1);
    
%% CALCULATE AGE BASED ON SCHWANDER ET AL 2001 ice flow model
    if accumFlag ~= 4
        acc=nan(lz,1);  % layer thickness for the Schwander Model
        HH=length(z)-1; % loop index
        
        % Define depth to the shear layer, hh
        if accumFlag == 0  
            h = param(2,1); 
            hh=(1-h)*H; 
        elseif accumFlag == 1 || accumFlag == 5
            hh=H-1200; 
        elseif accumFlag == 2
            h = param(2,1);
            hh=(1-h)*H; 
        elseif accumFlag == 3   
            h = param(3,1);               
            hh=(1-h)*H;                     
        elseif accumFlag == 6 || accumFlag == 7 
            h = param(6,1);                 
            hh=(1-h)*H;                                      
        elseif accumFlag == 100
            h = param(2,1);
            hh=(1-h)*H;
        end
        
        % Define other Schwander ice flow model parameters
        q=param(1,1);   % Schwander's velocity ratio  (calling it q like in the paper)
        k = 2 / (2*H - hh*(1-q)); % parameter in Schwander's model (see paper)
        age(lz,1) = 0;            % age at surface is 0
        acc(H) = 0.11;            % modern accum at byrd is 11 cm/yr
        accum0=acc(H);
        
        % Define accumulations depending on the accumFlag
        if accumFlag == 5 || accumFlag == 6 || accumFlag == 7
            a = H - 150;    % these are H - x because I'm reversing the origin
            b = H - 1024;   % the values here correspond to thresholds in Hammer et al 1994's
            c = H - 1200;   % piecewise linear accumulation function

            % Assign layer thickness based on parameters selected by Metropolis
            acc(z>a)        = param(5,1);         % shallowest part
            acc(z<=a & z>b) = param(4,1);
            acc(z<=b & z>c) = param(3,1);
            acc(z<=c)       = param(2,1);         % deepest part
           
        elseif accumFlag == 3 || accumFlag == 1
            acc(1:H) = param(2,1); % a single value for the whole column
            
        elseif accumFlag == 0      % variable accum for every meter of the column
            acc(1:H) = param(3:end,1); 
        end
       
        if accumFlag == 2 % Morse accumulation functional form (requires quadratic formula to solve)
                pt1=[10000,1];
                pt2=[17000,0.4];
                pt3=[80000,1];                         

            for i=HH:-1:1 % start just below the surface and integrate toward the bed
                if (i >= hh) && (i <= HH)                           % upper part of the ice (high values of i,z)
                    strainrate_z = 1 - k*(H-z(i));
                elseif (i > 0) && (i < hh)                          % lower part of the ice (low values of z)
                    strainrate_z = k*z(i)*(q+(1-q)/(2*hh)*z(i));    % no need for a dz term because dz=1 at each step
                end 
                if age(i+1,1) <= 10000 || age(i+1,1) > 80000 % during interglacial times
                    age(i,1) = age(i+1,1) + 1/(strainrate_z*accum0); %this is the typical D-J relation (from Schwander et al 2001)
                else %computing age given morse et al. function of accum
                    if age(i+1,1) > 10000 && age(i+1,1) <= 18000
                        %need the quadratic formula here, so defining the params a, b, c for that
                        slope = (pt1(2)-pt2(2))/(pt1(1)-pt2(1));
                        intercept =  pt1(2)-slope*pt1(1);
                    elseif age(i+1,1) > 18000 && age(i+1,1) <= 80000
                        slope = (pt3(2)-pt2(2))/(pt3(1)-pt2(1));
                        intercept = pt3(2)-slope*pt3(1);
                    end
                    a=slope;
                    b=intercept-age(i+1,1)*slope;
                    c=-1/(strainrate_z*accum0)-intercept*age(i+1,1);

                    age_tmp = (-b + sqrt(b^2-4*a*c))/(2*a);  % quadratic formula
                    age_tmp2 = (-b - sqrt(b^2-4*a*c))/(2*a); % other half of the quadratic formula

                    if isreal(age_tmp)  && age_tmp >=  age(i+1,1)% && ~isreal(age_tmp2)
                         age(i,1) = age_tmp;
                    elseif isreal(age_tmp2) && age_tmp2 >= age(i+1,1) %&& ~isreal(age_tmp)
                        age(i,1) = age_tmp2;
                    else
                        disp('You''ve got problems.')
                        disp(age(i+1,1))
                        disp(age_tmp)
                        disp(age_tmp2)
                    end
                end
            end
        else % No need for the quadratic formula to use these other accum functions
            for j=HH:-1:1
                if accumFlag == 100 % do this here because it's easier in the loop
                    if age(j+1,1) <= 5000
                        acc(j) = param(3);
                    elseif age(j+1,1) <= 10000
                        acc(j) = param(4);
                    elseif age(j+1,1) <= 15000
                        acc(j) = param(5);
                    elseif age(j+1,1) <= 20000
                        acc(j) = param(6);
                    elseif age(j+1,1) <= 25000
                        acc(j) = param(7);
                    elseif age(j+1,1) <= 30000
                        acc(j) = param(8);
                    elseif age(j+1,1) <= 35000
                        acc(j) = param(9);
                    elseif age(j+1,1) <= 40000
                        acc(j) = param(10);
                    elseif age(j+1,1) <= 45000
                        acc(j) = param(11);
                    elseif age(j+1,1) <= 50000
                        acc(j) = param(12);
                    end
                end
                % this is the strain part of the equation
                if j >=hh && j <= H
                    sz = 1 - k*(H-z(j));
                elseif j < hh && j > 0
                    sz = k*z(j)*(q+(1-q)/(2*hh)*z(j));    
                end
                % this is the age equation according to Schwander et al 2001
                age(j,1) = age(j+1,1) + (1/(sz*acc(j)))*(z(j+1)-z(j));
            end
        end
    end
%% CALCULATE AGE BASED ON MORLAND ET AL 2009 ice flow model
    if accumFlag == 4 
        hc=H;   % core length (meters) from Montagnat and Duval 2000
        h0=hc;  % This is the surface, not a Nye thing.
        b0= 0;  % 0.007065; % basal melt rate, inferred from temperature profile
                % and water in the borehole at the time of drilling

        a=150;  % thresholds from Hammer for layer thickness regimes;
        b=1024; % the origin for the Morland model is now at the base of the 
        c=1294; % ice sheet, so these are reversed from Schwander

        zMorland = H - z; % establishing coord system where z = 0 is at the base
        q=nan(lz,1);      % q is accumulation in the Morland model

        if accumFlag == 4
            q(zMorland  < a)                = param(1,1); % deepest part
            q(zMorland >= a & zMorland < b) = param(2,1);
            q(zMorland >= b & zMorland < c) = param(3,1);
            q(zMorland >= c)                = param(4,1); % shallowest part
        elseif accumFlag == 1    % constant accum
            q(1:H) = param(1,1); % a single value for the whole column
        end
        s0=(q-b0)/h0; % see the paper about this one
        %sd=param(1,1);     % velocity ratio
        sd = 0.722; % from morland paper. 
        s=s0*sd;         % again, see the paper.
        r=b0./q;       % parameter defined for the model

        % Evaluation of the age equation from Morland et al. 2009
        for i=lz:-1:1
           if i == lz
               %starting value for the age
               age(i) = - (1/s(i))*(log(1-(zMorland(i)/h0)*(1-r(i))));
           else
               %integrate age; second term is d/dz*dz
               age(i) = age(i+1) + (1/s(i))*((1-r(i))/h0)*(1/(1-(zMorland(i)/h0)*(1-r(i))))*(zMorland(i)-zMorland(i+1));
           end
        end
    end 
end

