%% Improvement of reaction and metabolite annotation
% Use preliminary functionality from RAVEN, available from add_MetaNetX
% branch (https://github.com/SysBioChalmers/RAVEN/tree/feat/add_MetaNetX),
% which queries metabolites and reactions to MetaNetX and subsequently
% incorporates the relevant annotations.
load('../../scrap/model_r7.mat');

MNXref=buildMNXref('both');
MNXfields=mapToMNX(model,true,MNXref,false)

newModel = addMNXannot(modelCb,MNXfields,MNXref);
newModel = convertMiriams(newModel);
model.metMiriams = newModel.metMiriams;
model.rxnMiriams = newModel.rxnMiriams;

save('../../scrap/model_r8.mat','model');
cd('..'); newCommit(model); cd('reconstruction')
