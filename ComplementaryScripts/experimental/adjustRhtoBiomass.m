%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create
%
% Usage: model=
%
% Eduard Kerkhoven. Last update: 2018-07-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newModel, k] = adjustRhtoBiomass(model,data)

% Change lipid backbone composition
rxnIdx                  =   getIndexes(model,'r_4063','rxns');
metIdx                  =   getIndexes(model,data.lipidData.metIDs,'mets');
bbIdx                   =   getIndexes(model,'s_3746','mets');
model.S(:,rxnIdx)       =   0;
model.S(metIdx,rxnIdx)  =   -data.lipidData.abundance-data.lipidData.std;
model.S(bbIdx,rxnIdx)   =   1;

% Change lipid chain length distribution
rxnIdx                  =   getIndexes(model,'r_4065','rxns');
metIdx                  =   cell2mat(getIndexes(model,strcat('C',data.chainData.chain,' chain'),'metnames'));
bbIdx                   =   getIndexes(model,'s_3747','mets');
model.S(:,rxnIdx)       =   0;
% Normalize lipid chain data to lipid backbone level, to get reasonable
% estimate before final scaling
scaling = sum([data.lipidData.abundance;data.lipidData.std]);
scaling = scaling / sum([data.chainData.abundance;data.chainData.std]);
scaling = (data.chainData.abundance + data.chainData.std) * scaling;
model.S(metIdx,rxnIdx)  =   -scaling;
model.S(bbIdx,rxnIdx)   =   1;

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model = setParam(model,'ub',[chainExIdx,backbExIdx],1000);

sol=solveLP(model,1)

[newModel,k] = scaleAbundancesRhto(model,data,'tails');
end