function age = byrdModels2(param,z,H,accumFlag,lp,varargin)
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
    
    %if there's no depth to evaluate to, go to the bottom
    if nargin < 6
        z_eval = 1;
    else
        z_eval = round(varargin{1});
    end
    
    HH=H-1; % loop index
    lz=length(z)-z_eval+1;
    age = nan(lz,1);
    model ='schwander';
        
    
%% CALCULATE AGE BASED ON SCHWANDER ET AL 2001 ice flow model
    if strcmp(model,'schwander')
        if accumFlag == 7 
            % Define thickness of the shear layer, hh 
            h = param(length(param(:,1))-lp-2,1); %percent of ice thickness                
            hh = h*H; %shear layer thickness
        end

        % Define other Schwander ice flow model parameters
        q = param(1,1);             % Schwander's velocity ratio  (call q like in the paper)
        k = 2 / (2*H - hh*(1-q)); % parameter in Schwander's model (see paper)

        if accumFlag == 7 
            % Define accumulation layers in depth
            a = 200;
            b = 400;
            c = 600;
            d = 800;
            e = 1000;
            f = 1200;
            g = 1400;
            o = 1600;
            t = 1800;
                      
            % Set accumulation layer bounds in z
            az  = H - a;   % these are H - x because I'm reversing the origin.    
            bz  = H - b;   % values correspond to thresholds in Hammer 94
            cz  = H - c;   % piecewise linear layther thickness function
            dz  = H - d;
            ez  = H - e;
            fz  = H - f;
            gz  = H - g;
            oz  = H - o;
            tz  = H - t;
                      
            % Assign accumulation based on parameters selected by Metropolis
            acc(z > az)           = param(11,1); %shallowest part
            acc(z <= az & z > bz) = param(10,1);
            acc(z <= bz & z > cz) = param(9,1);
            acc(z <= cz & z > dz) = param(8,1);
            acc(z <= dz & z > ez) = param(7,1);
            acc(z <= ez & z > fz) = param(6,1);
            acc(z <= fz & z > gz) = param(5,1);
            acc(z <= gz & z > oz) = param(4,1);
            acc(z <= oz & z > tz) = param(3,1);
            acc(z <= tz)          = param(2,1);  %deepest part

            %acc(H) = 0.14;            % surf accum at byrd
            acc(H) = param(11,1); %set surface to be same as shallowest
        end

        age(H,1) = 0;            % age at surface is 0
        for j=HH:-1:z_eval

            % this is the strain part of the equation
            if j >= hh && j <= HH  % for the upper part of the ice col
                sz = 1 - k*(H-z(j));
            elseif j < hh && j > 0 % for the shear layer part 
                sz = k*z(j)*(q+(1-q)/(2*hh)*z(j));
            else
                disp('Houston we have a problem')
                hh
                HH
                j
            end
            % this is the age equation according to Schwander et al 2001
            age(j,1) = age(j+1,1) + (1/(sz*acc(j)))*(z(j+1)-z(j));
        end
    end    
%%  CALCULATE AGE BASED ON MORLAND ET AL 2009 ice flow model
    if strcmp(model,'morland')
        hc = H;   % core length (meters) from Montagnat and Duval 2000
        h0 = hc;  % This is the surface, not a Nye thing.
        b0 = 0;  % 0.007065; % basal melt rate, inferred from temperature profile
                % and water in the borehole at the time of drilling

        a = 150;  % thresholds from Hammer for layer thickness regimes;
        b = 1024; % the origin for the Morland model is now at the base of the 
        c = 1294; % ice sheet, so these are reversed from Schwander

        az = H - a;   % these are H - x because I'm reversing the origin.    
        bz = H - b;   % values correspond to thresholds in Hammer 94
        cz = H - c;   % piecewise linear layther thickness function

        %zMorland = H - z; % establishing coord system where z = 0 is at the base
    %     if strcmp(accum,'hammer')
            q=nan(lz,1);      % q is accumulation in the Morland model
            q(z > az)           = param(5,1); % shallowest part
            q(z <= az & z > bz) = param(4,1);
            q(z <= bz & z > cz) = param(3,1);
            q(z <= cz)          = param(2,1); % deepest part
            q(H)                = 0.14;       %according to Morland et al
    %     end
        s0=(q-b0)/h0; % see the paper about this one
        %sd=param(1,1);     % velocity ratio
        sd = 0.722;   % from morland paper. 
        s=s0*sd;      % again, see the paper.
        r=b0./q;       % parameter defined for the model

        % Evaluation of the age equation from Morland et al. 2009
        age(H,1) = 0;
        for i=HH:-1:z_eval
           %integrate age; second term is d/dz*dz
           age(i) = age(i+1) + (1/s(i))*((1-r(i))/h0)*(1/(1-(z(i)/h0)*(1-r(i))))*(z(i)-z(i+1));
        end
    end

%% CALCULATE AGE BASED ON MORLAND ET AL 2009 ice flow model
    if strcmp(model,'morland')
%         accumFlag = 4;
        hc=H;   % core length (meters) from Montagnat and Duval 2000
        h0=hc;  % This is the surface, not a Nye thing.
        b0= 0;  % 0.007065; % basal melt rate, inferred from temperature profile
                % and water in the borehole at the time of drilling

%         a=150;  % thresholds from Hammer for layer thickness regimes;
%         b=1024; % the origin for the Morland model is now at the base of the 
%         c=1294; % ice sheet, so these are reversed from Schwander

        a = 200;
        b = 400;
        c = 600;
        d = 800;
        e = 1000;
        f = 1200;
        g = 1400;
        o = 1600;
        t = 1800;
        
        az  = H - a;   % these are H - x because I'm reversing the origin.    
        bz  = H - b;   % values correspond to thresholds in Hammer 94
        cz  = H - c;   % piecewise linear layther thickness function
        dz  = H - d;
        ez  = H - e;
        fz  = H - f;
        gz  = H - g;
        oz  = H - o;
        tz  = H - t;

        zMorland = H - z; % establishing coord system where z = 0 is at the base
        q=nan(lz,1);      % q is accumulation in the Morland model

%         if accumFlag == 4
            q(zMorland  < az)                 = param(11,1); % shallowest part
            q(zMorland >= az & zMorland < bz) = param(10,1);
            q(zMorland >= bz & zMorland < cz) = param(9,1);
            q(zMorland >= cz & zMorland < dz) = param(8,1);
            q(zMorland >= dz & zMorland < ez) = param(7,1);
            q(zMorland >= ez & zMorland < fz) = param(6,1);
            q(zMorland >= fz & zMorland < gz) = param(5,1);
            q(zMorland >= gz & zMorland < oz) = param(4,1);
            q(zMorland >= oz & zMorland < tz) = param(3,1);
            q(zMorland >= tz)                 = param(2,1); % deepest part
            
%         elseif accumFlag == 1    % constant accum
%             q(1:H) = param(1,1); % a single value for the whole column
%         end
        s0=(q-b0)/h0; % see the paper about this one
        sd=param(1,1);     % velocity ratio
        %sd = 0.722;   % from morland paper. 
        s=s0*sd;      % again, see the paper.
        r=b0./q;      % parameter defined for the model

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
    
%%
if strcmp(model,'morse')

    % Define depth to the shear layer, hh 
    h = param(6,1);                 
    hh=(1-h)*H;  
    % Define other Schwander ice flow model parameters
    q=param(1,1);   % Schwander's velocity ratio  (call q like in the paper)
    k = 2 / (2*H - hh*(1-q)); % parameter in Schwander's model (see paper)
    age(H,1) = 0;

    pt1=[10000,1];
    pt2=[17000,0.4];
    pt3=[80000,1];
    accum0=0.11;

    for i=HH:-1:z_eval % start just below the surface and integrate toward the bed
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
end
if nargin == 6
    age = age(z_eval);
    %disp(nargin)
%     else
%         disp(nargin)
%         disp('nope')
end
    
    
end

