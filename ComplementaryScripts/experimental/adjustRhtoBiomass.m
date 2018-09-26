%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create
%
% Usage: model=
%
% Eduard Kerkhoven. Last update: 2018-07-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load model
model=importModel('../../Rhodosporidium_toruloides-GEM/ModelFiles/xml/rhto.xml');

data=readTiukovaData(1);

% Change lipid backbone composition
rxnIdx                  =   getIndexes(model,'r_4063','rxns');
metIdx                  =   getIndexes(model,data.lipidData.metIDs,'mets');
bbIdx                   =   getIndexes(model,'s_3746','mets');
model.S(:,rxnIdx)       =   0;
model.S(metIdx,rxnIdx)  =   -data.lipidData.abundance;
model.S(bbIdx,rxnIdx)   =   1;

% Change lipid chain length distribution
rxnIdx                  =   getIndexes(model,'r_4065','rxns');
metIdx                  =   getIndexes(model,data.chainData.metIDs','mets');
bbIdx                   =   getIndexes(model,'s_3747','mets');
model.S(:,rxnIdx)       =   0;
model.S(metIdx,rxnIdx)  =   -data.chainData.abundance;
model.S(bbIdx,rxnIdx)   =   1;

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model = setParam(model,'ub',[chainExIdx,backbExIdx],1000);

model = setParam(model,'obj','r_2111',1);
sol=solveLP(model,1);

[model,k] = scaleAbundancesRhto(model,data,'tails');