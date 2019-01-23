%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sol,model] = simulateGrowth(model,fluxData)
%
%   Simulates growth as part of the scaleAbundancesRhto script.
%   Modified from SLIMEr, under MIT License:
%   https://github.com/SysBioChalmers/SLIMEr/blob/master/simulations/simulateGrowth.m
%
% 2019-01-21    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol,model] = simulateRhtoGrowth(model,fluxData)

%Constrain all fluxes with exp. data
%The variation should never be more than the abundance: 
%stdev = min([fluxData.stdevs,abs(fluxData.averages)/1.96],[],2);

%The variation should be at least a 5% (to avoid unfeasible simulations):
stdev = abs(fluxData.averages)*0.05; %max([stdev,abs(fluxData.averages)/1.96*0.05],[],2);
LB    = fluxData.averages - stdev;      %C.I. of 95%
UB    = fluxData.averages + stdev;      %C.I. of 95%
model = setParam(model,'lb',fluxData.rxnIDs,LB);
model = setParam(model,'ub',fluxData.rxnIDs,UB);

% Remove GAM from biomass equation
mets  = getIndexes(model,{'ATP[c]','ADP[c]','phosphate[c]','H2O[c]','H+[c]'},'metscomps');
bmIdx = getIndexes(model,'r_4041','rxns');
model.S(mets,bmIdx) = 0;

% Maximize maintenance while reducing number of fluxes
ngamIdx = getIndexes(model,'r_4046','rxns');
model = setParam(model,'lb',ngamIdx,0);
model = setParam(model,'ub',ngamIdx,1000);
model = setParam(model,'obj',ngamIdx,1);
sol   = solveLP(model,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

