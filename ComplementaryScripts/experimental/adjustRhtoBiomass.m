%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newModel, k] = adjustRhtoBiomass(model,data)
%
%   Adjusts the biomass equation to match the lipid backbone and chain
%   length distribution data that is provided. Chain length data is scaled
%   to agree with the lipid class measurements.
%
% 2019-01-21    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newModel, k] = adjustRhtoBiomass(model,data)

% Change lipid backbone composition
rxnIdx                  =   getIndexes(model,'r_4063','rxns');
metIdx                  =   getIndexes(model,data.lipidData.metIds,'mets');
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
scaling = scaling / sum([data.chainData.abundance;data.chainData.std],'omitnan');
scaling = (data.chainData.abundance + data.chainData.std) * scaling;
model.S(metIdx,rxnIdx)  =   -scaling;
model.S(bbIdx,rxnIdx)   =   1;

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model = setParam(model,'ub',[chainExIdx,backbExIdx],1000);

[newModel,k] = scaleLipidsRhto(model,'tails');
end