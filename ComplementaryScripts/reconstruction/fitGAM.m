%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAM = fitGAM(model)
%
%   Modified from GECKO, under MIT License:
%   https://github.com/SysBioChalmers/SLIMEr/blob/master/models/scaleAbundancesInModel.m
%
% 2019-01-22    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,GAM] = fitGAM(model)

%Load chemostat data:
fid         = fopen('../../ComplementaryData/data/validationData.csv');
fluxData    = textscan(fid,'%f32 %f32 %s %s','Delimiter',',','HeaderLines',1);
fluxData    = [fluxData{1} fluxData{2}];
fclose(fid);
fluxData(17:end,:) = [];

%GAMs to span:
disp('Estimating GAM:')
GAM = 70:15:160;

%1st iteration:
GAM = iteration(model,GAM,fluxData);

%2nd iteration:
GAM = iteration(model,GAM-10:1:GAM+10,fluxData);

%3rd iteration:
GAM = iteration(model,GAM-1:0.1:GAM+1,fluxData);

model = setGAM(model,GAM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GAM = iteration(model,GAM,fluxData)

fitting = ones(size(GAM))*1000;

for i = 1:length(GAM)
    %Simulate model and calculate fitting:
    mod_data   = abs(simulateChemostat(model,fluxData,GAM(i)));
    R          = (mod_data - fluxData)./fluxData;
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(['GAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[~,best] = min(fitting);

if best == 1 || best == length(GAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    GAM = GAM(best);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mod_data = simulateChemostat(model,fluxData,GAM)

model = setGAM(model,GAM);

%pos = getIndexes(model',{'r_4041','r_1714','r_1654'},'rxns');
pos = getIndexes(model',{'r_4041','r_1714'},'rxns');

%Simulate chemostats:
mod_data = zeros(size(fluxData));
for i = 1:7
    %Fix glucose and maximize biomass
    model = setParam(model,'lb',pos(2),-fluxData(i,2));
    model = setParam(model,'ub',pos(2),0);
    model = setParam(model,'obj',pos(1),1);
    sol   = solveLP(model);
    %Store relevant variables:
    mod_data(i,:) = sol.x(pos)';

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = setGAM(model,GAM)

xr_pos = getIndexes(model','r_4041','rxns');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,xr_pos) = sign(S_ix) * GAM;
    end
end
end