function [D,H,obsAge1950,obsAge1950Unc] = interpWDobs(WD)
% Read in the lengthy chronology (with unc) of WD ice core and reduce it to
% age-depth evaluated at every 1 m

    [D,ind]=unique(WD(:,1)); %some of the depths are repeats, so take only the unique values
    obsAge1950=WD(ind,2);   
    obsAge1950Unc=WD(ind,3);

    H = round(max(D)); % ice thickness

    %% get the chronology and unc only at every 1m
    interpAge = interp1(D,obsAge1950,(1:H)');
    interpAgeUnc= interp1(D,obsAge1950Unc,(1:H)');
    
    %% interp gives nan for the last point, so subsitute reasonable endpoint
    interpAge(end) = interpAge(end-1) + (interpAge(end-1)-interpAge(end-2));
    interpAgeUnc(end) = interpAgeUnc(end-1) + (interpAgeUnc(end-1)-interpAgeUnc(end-2));
    
    obsAge1950=interpAge;
    obsAge1950Unc=interpAgeUnc;

end

