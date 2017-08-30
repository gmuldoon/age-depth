function [keAge,Em,ageBest,ageSamp,paramSamp] = calculate_keAge(cost,Param,z,H,accumFlag,lp)
%% Select param values

%Find parameters corresponding to lowest-cost iteration
[costSort,minind] = min(cost(:,2)); %minimum of age cost
paramBest = Param(:,minind);



% best guess std on parameters
paramStd(1,1)  = 0.05;
paramStd(2:11,1) = 0.03;
paramStd(12,1) = 0.1;
paramStd(13,1) = 1e5;
paramStd(14,1) = 10;
paramStd(15:18,1) = 10;

%% Find depths of volcanic data at which to evaluate model
filename='../data/byrd_volcage.txt';
ncol=2;
byrdObs=readFloats(filename,ncol);
D1968=byrdObs(:,1);
accum = 0.11;
D = D1968+(2013-1968)*accum;
[Dtrim,~] = unique(ceil(D)); %Only use those depths which are different when rounded
Zvolc = H - Dtrim;
  
%% Compute keAge
% ke = 2*<E>^2/var(E)

%compute perfect model of ages
ageBest = byrdModels2(paramBest,z,H,accumFlag,lp);
cjparam(:,1) = zeros(1,1);
cjparam(:,2) = ones(1,1);

%Resample parameters and compute modeled ages
nsamp = 1000;
ageSamp = zeros(H,nsamp);
csj_m = zeros(nsamp,1);
for i = 1:nsamp
    paramSamp = normrnd(paramBest,paramStd/10);
    ageSamp(:,i) = byrdModels2(paramSamp,z,H,accumFlag,lp);      
    cj_samp = normrnd(cjparam(:,1),cjparam(:,2));
    csj_m(i) = sum(cj_samp.^2/2);
end
varSamp = var(ageSamp(Zvolc,:),1,2);
sqrt(varSamp)
Em = sum((ageSamp(Zvolc,:) - ageBest(Zvolc)).^2./(2*varSamp));

%compute var(E)
ageVar = var(Em)
varcj = var(csj_m);

%compute <E>
ageMean = mean(Em)

%Calc ke
keAge = 2*ageMean^2/ageVar;
kecj = 2*(mean(csj_m))^2/varcj

figure(60)
clf
plot(Em)

figure(61)
clf
hist(Em,20)

Zvolc;
figure
clf
plot(std(ageSamp(Zvolc,:),1,2))

figure
clf
hist(ageSamp(269,:)')

%ageSamp(10,:)
end



 
 

