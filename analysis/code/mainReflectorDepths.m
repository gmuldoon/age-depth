
%Calculate reflector depth at various sites

plotting = 0;


%pik TWTT to use
% core_S_T
Byrd_BB01a_13190 = [844,958,1254,1755,2242];
% Byrd_BB01a_13125 = [844,958,1272,1725];
% Byrd_BB01a_12950 = [809,914,1219,1674];
% Byrd_BB01a_13480 = [914,1054,1394,1882];
% WD_ARCH1a_9050 = [1555,1709,2266,3108];
% 

reflector_LM_times{1,:} = Byrd_BB01a_13190;
% reflector_LM_times{2,:} = Byrd_BB01a_13125;
% reflector_LM_times{3,:} = Byrd_BB01a_12950;
% reflector_LM_times{4,:} = Byrd_BB01a_13480;
% reflector_LM_times{5,:} = WD_ARCH1a_9050;
% 
pik_depth = nan(length(reflector_LM_times), length(Byrd_BB01a_13190));

for n = 1:length(reflector_LM_times)
   pik = reflector_LM_times{n}/100.;
   
   if n < 5
       core = 'byrd';
       H = 2164;
   else
       core = 'wd';
       %%Get H
       %WD ice core chronology (from Buizert, personal email): depth, age, 2sig age unc
        filename='../data/WD_chronologyunc.csv';
        ncol=3;
        WD=readFloats(filename,ncol);
        %interpolate WD chronology to be at 1m intervals
        [~,H,~,~] = interpWDobs(WD);
   end
   [~,~,pik_depth(n,:),~,~] = radarDepth(pik,core,H,plotting);
end
    


