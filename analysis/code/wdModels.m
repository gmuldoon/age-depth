function age = wdModels(param,z,H,accumFlag)
% This function uses the Schwander et al. 2001 and Parrenin et al 2007
% method to invert for flow parameters based on matching to the observed
% ice core chronology. 

%Note that z=0 at bedrock and increases upward (z=H is at the ice surface).
persistent WDdep WDz accum
%% Initialize parameters and set some constants
%read in the params and set up the age vector.
q=param(1);
if accumFlag == 1 
    accum=param(3);                                                        % sample a constant accumulation rate
elseif accumFlag == 2
    accum0 = 0.22;                                                         % modern surface accumulation at WD in m/yr
end
%h=0; % I think this should theoretically be 0 if the bed is sliding...WD
%ice core papers show melting at the bottom but assume h=0.2? Read more...
h=param(2);
%h=0.28;
%h = 0.7;                                                                   % set value of D-J kink for now. will invert for it later.
HH=H-1;                                                                    % useful index for just below the surface (already know age at surface)

age=nan(H,1);                                                              % age profile for the full ice column will be contained here
age(H,1) = 0;                                                              % age at the surface is constrained to be 0

%% Calculate modeled age
%points from Morse et al. 2002 which determine accumulation function
pt1=[10000,1];
pt2=[17000,0.4];
pt3=[80000,1];

hh=(1-h)*H; %h is fraction above the bed, but hh from Schwander et al is depth to the transition
k=2/(2*H - hh*(1-q));                                                   % parameter in the strain rate calculation, defined in Schwander et al. 2001

for i=HH:-1:1 % start just below the surface and integrate toward the bed
    if (i >= hh) && (i <= HH)                                              % this is the upper part of the ice (high values of i,z)
        strainrate_z = 1 - k*(H-z(i));
    elseif (i > 0) && (i < hh)                                             % this is the lower part of the ice (low values of z)
        strainrate_z = k*z(i)*(q+(1-q)/(2*hh)*z(i));                             %note there is no need for a dz term because dz=1 at each step
    end    
    if accumFlag == 1                                                           % when inverting for accum
        age(i,1) = age(i+1,1) + 1/(strainrate_z*accum);                    % integrating the age starting from just below ice surface H (corresponding to index HH)
    elseif accumFlag == 2                                                       % when assuming we know the functional form of accum over time
        if age(i+1,1) <= 10000 || age(i+1,1) > 80000 %ie during interglacial times
            age(i,1) = age(i+1,1) + 1/(strainrate_z*accum0); %this is the typical DJ relation, from Schwander et al 2001, eg
        else %computing age given morse et al. function of accum
            if age(i+1,1) > 10000 && age(i+1,1) <= 18000
                %need the quadratic formula here, so defining the params a, b,c for that
                slope = (pt1(2)-pt2(2))/(pt1(1)-pt2(1));
                intercept =  pt1(2)-slope*pt1(1);
            elseif age(i+1,1) > 18000 && age(i+1,1) <= 80000
                slope = (pt3(2)-pt2(2))/(pt3(1)-pt2(1));
                intercept = pt3(2)-slope*pt3(1);
            end
            a=slope;
            b=intercept-age(i+1,1)*slope;
            c=-1/(strainrate_z*accum0)-intercept*age(i+1,1);
            
            age_tmp = (-b + sqrt(b^2-4*a*c))/(2*a); %this is the quadratic formula
            age_tmp2 = (-b - sqrt(b^2-4*a*c))/(2*a); %this is the other half of the quadratic formula
     
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

if accumFlag == 3
        lthk=nan(H,1); % layer thickness for the Schwander Model

        a = H - 0;   % these are H - x because origin to be at the bed, H is at surface, but numbers correspond to depth below surface
        b = H - 1150;  % the values here correspond to thresholds in Hammer's
        c = H - hh;  % piecewise linear accumulation function

        % Assign layer thickness based on parameters selected by Metropolis
        for i=HH:-1:1
            if i > a
                lthk(i) = param(5,1);        % this is the shallowest part (high index values means closer to surface, H)
            elseif (i <= a) && (i > b)
                lthk(i) = param(4,1);
            elseif (i <= b) && (i > c)
                lthk(i) = param(3,1);
            elseif (i <= c)
                lthk(i) = param(2,1);        % this is the deepest part
            end   
        end

        % Calculate age based on Schwander model. Two loops are for above and
        % below the shear zone at h.
        for j=HH:-1:1
            if j >=hh && j <= HH   % upper part of the ice sheet, high values of z, j  
                strainrate_z = 1 - k*(H-z(j));
                age(j,1) = age(j+1,1) + (1/(strainrate_z*lthk(j)))*(z(j+1)-z(j));
            elseif j < hh && j > 0
                strainrate_z = k*z(j)*(q+(1-q)/(2*hh)*z(j));
                age(j,1) = age(j+1,1) + (1/(strainrate_z*lthk(j)))*(z(j+1)-z(j));
            end
        end

        %Change depth coord so depth = 0 at the surface
        d = H - z;
elseif accumFlag == 4
    if isempty(WDdep) || isempty(WDz) || isempty(accum); %empty instead of exists bc setting as persistent var makes it exist, but it's still empty
        disp('loading WD accum data')
        fid = fopen('../data/WD2014_Accumulation_9-11-2014_noheader.csv');
        WD=cell2mat(textscan(fid,'%f %f %f','Delimiter',','));
        fclose(fid);
   
        accum = WD(i,3); % accum in m ice equiv per year (from inverse firn densification model)
        WDdep = WD(:,1); % depth in m below the surface, ice age in years before 1950,
        WDz = flipud(H - WDdep); %changing to a coord sys with z=0 at the bed, want WDz to be large at high indices   
    end
    
    %calculate the age to each WD accum age
    WDage = zeros(length(WDdep),1);
    for i=length(WDdep)-1:-1:1 % -2 because the last data point is just slightly more than H, so ignore it
        if WDz(i) >= hh && WDz(i) <= HH  
            strainrate_z = 1 - k*(H-WDz(i));
        elseif WDz(i) < hh && WDz(i) > 0
            strainrate_z = k*WDz(i)*(q+(1-q)/(2*hh)*WDz(i));  
        else
            disp('')
            disp('Danger, Will Robinson! Not at a valid depth in WD accum data')
        end
        WDage(i,1) = WDage(i+1,1) + (1/(strainrate_z*accum));%*(WDz(i+1)-WDz(i));
    end
    
    age = WDage; %for the sake of the cost calculation
    
end

%% Interpolate the age-depth profile at each 1m in the ice core
if accumFlag == 4
    interp_age = interp1(WDz,age,z);
    % interp leads nan for last point, so subsitute a reasonable endpoint
    interp_age(end) = interp_age(end-1) + (interp_age(end-1)-interp_age(end-2));
    age=flipud(interp_age);
else
    age=flipud(age);
end

