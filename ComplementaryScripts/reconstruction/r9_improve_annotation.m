clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Improvement of reaction and metabolite annotation
% Use preliminary functionality from RAVEN, available from add_MetaNetX
% branch (https://github.com/SysBioChalmers/RAVEN/tree/feat/add_MetaNetX),
% which queries metabolites and reactions to MetaNetX and subsequently
% incorporates the relevant annotations.
load([root '/scrap/model_r8.mat']);
load([root '/scrap/modelTemplate.mat']); % yeast-GEM 8.2.0

[match, matchIdx]   = ismember(modelSce.rxns,model.rxns);
model.rxnMiriams(matchIdx(match)) = modelSce.rxnMiriams(match);

[match, matchIdx]   = ismember(modelSce.mets,model.mets);
model.metMiriams(matchIdx(match)) = modelSce.metMiriams(match);

modelCb = ravenCobraWrapper(model);
MNXref  = buildMNXref('both');
MNXfields = mapToMNX(modelCb,true,MNXref,false);

newModel = addMNXannot(modelCb,MNXfields,MNXref);

newModel = convertMiriams(newModel);
model.metMiriams = newModel.metMiriams;
model.rxnMiriams = newModel.rxnMiriams;

save([root '/scrap/model_r9.mat'],'model');

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
cd('..'); newCommit(model); cd('reconstruction')
