function pikAge = pikAge(pikDepth,ageDistrib,z,numsteps,core,H)

    pikZ = H - pikDepth;
    
    for n=1:numsteps
        if strcmp(core,'byrd')
            % sample an index for age model
            samp_age = randi([1 numsteps]); 

            % evaluate model at each dePik depth to find corresponding age
            [~,ind] = intersect(z,round(pikZ(:,n)));
            pikAge(:,n) = ageDistrib(ind,samp_age);
            
            %deAge(:,n) = deAge(:,n); %put top of the ice column at lower indices to match dePik
            
        elseif strcmp(core,'wd')
            % evaluate model at each dePik depth to find corresponding age
            [~,ind] = intersect(z,round(pikZ(:,n)));
            %assume WD ages are the observed w/ gaussian errors from metrop 
            pikAge(:,n) = normrnd(obsAge1950(ind),obsAge1950Unc(ind));
        end
        %flip age bc here age is ordered by z index, which puts deeper ages
        %first in the volum
        pikAge(:,n) = flipud(pikAge(:,n));
    end

end

