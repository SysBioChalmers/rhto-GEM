%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sol,model] = simulateGrowth(model,fluxData)
%
%   Modified from SLIMEr, under MIT License:
%   https://github.com/SysBioChalmers/SLIMEr/blob/master/simulations/simulateGrowth.m
%
% 2018-09-25    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol,model] = simulateRhtoGrowth(model,fluxData)

%Constrain all fluxes with exp. data
%The variation should never be more than the abundance: 
stdev = min([fluxData.stdevs,abs(fluxData.averages)/1.96],[],2);

%The variation should be at least a 5% (to avoid unfeasible simulations):
stdev = max([stdev,abs(fluxData.averages)/1.96*0.05],[],2);
LB    = fluxData.averages - 1.96*stdev;      %C.I. of 95%
UB    = fluxData.averages + 1.96*stdev;      %C.I. of 95%
model = setParam(model,'lb',fluxData.rxnIDs,LB);
model = setParam(model,'ub',fluxData.rxnIDs,UB);

% Maximize maintenance while reducing number of fluxes
ngamIdx = getIndexes(model,'r_4046','rxns');
model = setParam(model,'lb',ngamIdx,0);
model = setParam(model,'ub',ngamIdx,1000);
model = setParam(model,'obj',ngamIdx,1);
sol   = solveLP(model,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
