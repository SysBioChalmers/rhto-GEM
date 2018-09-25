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
for i = 1:length(fluxData.rxnIDs)
    %The variation should never be more than the abundance: 
    stdev = min([fluxData.stdevs(i),abs(fluxData.averages(i))/1.96]);
    %The variation should be at least a 5% (to avoid unfeasible simulations):
    stdev = max([stdev,abs(fluxData.averages(i))/1.96*0.05]);
    LB    = fluxData.averages(i) - 1.96*stdev;      %C.I. of 95%
    UB    = fluxData.averages(i) + 1.96*stdev;      %C.I. of 95%
    model = setParam(model,'lb',fluxData.rxnIDs(i),LB);
    model = setParam(model,'ub',fluxData.rxnIDs(i),UB);
end

%Simulation 1: Maximize maintenance
model = setParam(model,'lb','r_4046',0);
model = setParam(model,'ub','r_4046',1000);
model = setParam(model,'obj','r_4046',1);
sol   = solveLP(model);

model    = changeRxnBounds(model,model.rxns(posMaint),0,'l');
model    = changeRxnBounds(model,model.rxns(posMaint),+1000,'u');
model    = changeObjective(model,model.rxns(posMaint),+1);
sol      = optimizeCbModel(model);

%Simulation 2: Force maintenance, minimize sum(abs(fluxes))
obj   = sol.x(posMaint);
model = changeRxnBounds(model,model.rxns(posMaint),obj*0.999,'l');
model = changeRxnBounds(model,model.rxns(posMaint),obj*1.001,'u');
sol   = optimizeCbModel(model,'min','one');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
