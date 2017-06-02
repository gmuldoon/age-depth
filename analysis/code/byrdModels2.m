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
        %disp(length(param(:,1))-lp-2)
        if accumFlag == 7 
            % Define depth to the shear layer, hh 
            h = param(length(param(:,1))-lp-2,1);                 
            hh=(1-h)*H;  
        end

        % Define other Schwander ice flow model parameters
        q=param(1,1);   % Schwander's velocity ratio  (call q like in the paper)
        k = 2 / (2*H - hh*(1-q)); % parameter in Schwander's model (see paper)

        if accumFlag == 7 
            % Define accumulation layers depending on the accumFlag
            a  = 100; % m depth
            b  = 200; % m depth
            c  = 300;% m depth
            d  = 400;
            e  = 500;
            f  = 600;
            g  = 700;
            o  = 800;
            t  = 900;
            u  = 950;
            v  = 1000;
            w  = 1050;
            x  = 1100;
            y  = 1150;
            aa = 1200;
            bb = 1250;
            cc = 1300;
            dd = 1350;
            ee = 1400;
            ff = 1450;
            gg = 1500;
            ll = 1550;
            mm = 1575;
            oo = 1600;
            pp = 1625;
            qq = 1650;
            rr = 1675;
            ss = 1700;
            tt = 1725;
            uu = 1750;
            vv = 1775;
            ww = 1800;
            xx = 1825;
            yy = 1850;
            aaa = 1875;
            bbb = 1900;
            ccc = 1950;
            ddd = 2000;
            eee = 2050;
            fff = 2100;
            ggg = 2150;
                      

            az  = H - a;   % these are H - x because I'm reversing the origin.    
            bz  = H - b;   % values correspond to thresholds in Hammer 94
            cz  = H - c;   % piecewise linear layther thickness function
            dz  = H - d;
            ez  = H - e;
            fz  = H - f;
            gz  = H - g;
            oz  = H - o;
            tz  = H - t;
            uz  = H - u;
            vz  = H - v;
            wz  = H - w;
            xz  = H - x;
            yz  = H - y;
            aaz = H - aa;
            bbz = H - bb;
            ccz = H - cc;
            ddz = H - dd;
            eez = H - ee;
            ffz = H - ff;
            ggz = H - gg;
            llz = H - ll;
            mmz = H - mm;
            ooz = H - oo;
            ppz = H - pp;
            qqz = H - qq;
            rrz = H - rr;
            ssz = H - ss;
            ttz = H - tt;
            uuz = H - uu;
            vvz = H - vv;
            wwz = H - ww;
            xxz = H - xx;
            yyz = H - yy;
            aaaz= H - aaa;
            bbbz= H - bbb;
            cccz= H - ccc;
            dddz= H - ddd;
            eeez= H - eee;
            fffz= H - fff;
            gggz=H - ggg;
          
            
            % Assign accumulation based on parameters selected by Metropolis
            acc(z >  az)            = param(40,1);    % shallowest part
            acc(z <= az & z > bz)   = param(39,1);
            acc(z <= bz & z > cz)   = param(38,1);
            acc(z <= cz & z > dz)   = param(37,1);         
            acc(z <= dz & z > ez)   = param(36,1);
            acc(z <= ez & z > fz)   = param(35,1);
            acc(z <= fz & z > gz)   = param(34,1);
            acc(z <= gz & z > oz)   = param(33,1);
            acc(z <= oz & z > tz)   = param(32,1);
            acc(z <= tz & z > uz)   = param(31,1);
            acc(z <= uz & z > vz)   = param(30,1);         
            acc(z <= vz & z > wz)   = param(29,1);
            acc(z <= wz & z > xz)   = param(28,1);
            acc(z <= xz & z > yz)   = param(27,1);
            acc(z <= yz & z > aaz)  = param(26,1);
            acc(z <= aaz & z > bbz) = param(25,1);
            acc(z <= bbz & z > ccz) = param(24,1);
            acc(z <= ccz & z > ddz) = param(23,1);
            acc(z <= ddz & z > eez) = param(22,1);
            acc(z <= eez & z > ffz) = param(21,1);     % deepest part
            acc(z <= ffz & z > ggz) = param(20,1);
            acc(z <= ggz & z > llz) = param(19,1);
            acc(z <= llz & z > mmz) = param(18,1);
            acc(z <= mmz & z > ooz) = param(17,1);
            acc(z <= ooz & z > ppz) = param(16,1);
            acc(z <= ppz & z > qqz) = param(15,1);
            acc(z <= qqz & z > rrz) = param(14,1);
            acc(z <= rrz & z > ssz) = param(13,1);
            acc(z <= ssz & z > ttz) = param(12,1);
            acc(z <= ttz & z > uuz) = param(11,1);
            acc(z <= uuz & z > vvz) = param(10,1);
            acc(z <= vvz & z > wwz) = param(9,1);
            acc(z <= wwz & z > xxz) = param(8,1);
            acc(z <= xxz & z > yyz) = param(7,1);
            acc(z <= yyz & z > aaaz) = param(6,1);
            acc(z <= aaaz & z > bbbz) = param(5,1);
            acc(z <= bbbz & z > cccz) = param(4,1);
            acc(z <= cccz & z > dddz) = param(3,1);
            acc(z <= dddz)            = param(2,1);
            
    %         acc(z > az)           = 0.14;         % shallowest part
    %         acc(z <= az & z > bz) = 0.14;
    %         acc(z <= bz & z > cz) = 0.13;
    %         acc(z <= cz)          = 0.11;         % deepest part

            acc(H) = 0.14;            % surf accum at byrd
        end
        %acc(:) = 0.14;

        age(H,1) = 0;            % age at surface is 0
        for j=HH:-1:z_eval
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

