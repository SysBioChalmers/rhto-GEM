%% Copy pseudoreactions
load('../../scrap/model_r1.mat');
load('../../scrap/modelTemplate.mat');

rxns=modelSce.rxns(contains(modelSce.rxnNames,'pseudoreaction'));
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'Modeling reaction',1); % Add reactions and metabolites

model=addRxnsGenesMets(model,modelSce,'r_4046',false,...
    'Modeling reaction',1); % Add reactions and metabolites

rxns=modelSce.rxns(contains(modelSce.rxnNames,'SLIME'));
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'SLIME reaction',1); % Add reactions and metabolites

%% Add all exchange rxns
% These were not gene annotated, and therefore not added in draft.
% Might not require all exchange rxns, but easier to remove unconnected ones later.
rxns=getExchangeRxns(modelSce);
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'Modelling reaction',1); % Add reactions and metabolites

%% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx=find(cellfun(@isempty,modelSce.grRules)); % Which rxns have no genes
rxnIdx=find(getTransportRxns(modelSce));
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns=modelSce.rxns(rxnIdx); % Obtain reaction IDs
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'Modeling reaction required for intercellular transport, gene unknown',1); % Add reactions and metabolites

save('../../scrap/model_r2.mat','model');
cd('..'); newCommit(model);cd('reconstruction')