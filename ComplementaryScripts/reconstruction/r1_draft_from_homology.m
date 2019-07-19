clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% This script generates a first draft of a genome-scale model (here: R.
% toruloides) using the getModelFromHomology function in RAVEN. It uses
% S. cerevisiae as a template model, 

% BLAST against Sce genome. Rhto protein fasta obtained from JGI, modified
% to have genes in format RHTO_00001. Sce protein fasta
%  blastSce=getBlast('rhto',[data '/genome/rhto_np11.faa'],'sce',[data '/genome/sce_s288c.faa']);
%  blastYli=getBlast('rhto',[data '/genome/rhto_np11.faa'],'yli',[data '/genome/yli_clib122.faa']);
% mkdir([root '/scrap'])
% save([root '/scrap/blastStruct.mat'],'blast*');
load([root '/scrap/blastStruct.mat']);

% Load S. cerevisiae model (corrected by Benjamin), downloaded from
% https://github.com/SysBioChalmers/yeast-GEM/raw/master/ModelFiles/xml/yeastGEM.xml
% (version 8.2.0).
modelSce        = importModel([data '/reconstruction/yeastGEM_820.xml'],true);
modelSce.id     = 'sce';
% Remove compartment information from metabolite ids (redundant).
modelSce.mets   = regexprep(modelSce.mets,'\[[a-z]+\]$','');

% Change reversibility fields to match boundaries. Prevents problems with
% MENECO.
for i=1:length(modelSce.rxns)
	if modelSce.lb(i)==0 && not(modelSce.ub(i)==0)
        modelSce.rev(i)=0;
    elseif modelSce.lb(i)<=0 && modelSce.ub(i)>=0
        modelSce.rev(i)=1;
	end
end
% Confirm that the model is functional, set objective to growth.
modelSce        = setParam(modelSce,'obj','r_4041',1);
solveLP(modelSce)

% Yarrowia 
modelYli        = importModel([data '/reconstruction/iYali_411.xml'],true);
modelYli.id     = 'yli';
[modelYli.grRules, modelYli.rxnGeneMat]=standardizeGrRules(modelYli);
for i=1:length(modelYli.rxns)
	if modelYli.lb(i)==0 && not(modelYli.ub(i)==0)
        modelYli.rev(i)=0;
    elseif modelYli.lb(i)<=0 && modelYli.ub(i)>=0
        modelYli.rev(i)=1;
	end
end
solveLP(modelYli)

modelYli.rxns = regexprep(modelYli.rxns,'y00','r_');
modelYli = removeReactions(modelYli,contains(modelYli.rxns,'y10'),true,true,true);

save([root '/scrap/modelTemplate.mat'], 'model*');

%% Generate draft model, based on homology.
model   = getModelFromHomology(modelSce,blastSce,'rhto',{},1,false,10^-20,150,35);

%% Make model based on MetaPhOrs data
orthologs = readtable([data '\reconstruction\MetaPhOrs_convertedIDs.csv']);
% Only keep orthologs if phylomeDB tree had more than 2 members and
% evidence is 1.
moreThanTwo         = str2double(regexprep(orthologs.phylome,'.*/([0-9]+)','$1'))>1;
evidence            = str2double(regexprep(orthologs.phylome(:),'\/[0-9]+',''))>=0.5;
keep                = find(evidence & moreThanTwo);
orthologList(:,1)   = orthologs.SaccharomycesCerevisiae(keep);
orthologList(:,2)   = orthologs.RhodosporidiumToruloides(keep);
blastStructure      = makeFakeBlastStructure(orthologList,'sce','rhto');
modelOrth           = getModelFromHomology(modelSce,blastStructure,'rhto');
modelOrth.id        = 'ortho';

modelComb           = mergeModels({model,modelOrth});
model               = contractModel(modelComb);

%% Add some reactions as based on homology with Yarrowia lipolytica
modelYli    = getModelFromHomology(modelYli,blastYli,'rhto',{},1,false,10^-20,150,35);

% Discard reactions that were already in draft rhto-GEM
modelYli    = removeReactions(modelYli,contains(modelYli.rxns,model.rxns),true,true,true);

% Focus on reactions derived from yeast-GEM
tmp         = removeReactions(modelYli,cellfun(@isempty,regexp(modelYli.rxns,'r_\d{4}$')),true,true,true);

% How the Yarrowia model was constructed, there is a set of new metabolites
% that were introduced by simplifying lipid metabolism. Discard these
% metabolites and associated reactions.
tmp         = removeMets(tmp,contains(tmp.mets,'m'),false,true,true,true);
tmp         = removeReactions(tmp,~ismember(tmp.rxns,modelSce.rxns));
model       = addRxnsGenesMets(model,modelSce,tmp.rxns,tmp.grRules,'Identified from homology to Yarrowia lipolytica',2);

% Add Yarrowia specific reactions
tmp                 = removeReactions(modelYli,~contains(modelYli.rxns,'y'),true,true,true);
% Replace old identifiers with new format
oldIdx              = find(contains(tmp.mets,'m'));
tmp.mets(oldIdx)    = generateNewIds(model,'mets','m_',numel(oldIdx));
tmp.metNames        = regexprep(tmp.metNames,'alpha-D-ribose 1-phosphate','alpha-D-ribose 1-phosphate(2-)');

modelComb           = mergeModels({model,tmp});
model               = contractModel(modelComb);

%% Add R. toruloides GEM meta data
model.annotation.defaultLB    = -1000;
model.annotation.defaultUB    = +1000;
model.annotation.taxonomy     = 'taxonomy/1130832';
model.annotation.givenName    = 'Eduard';
model.annotation.familyName   = 'Kerkhoven';
model.annotation.email        = 'eduardk@chalmers.se';
model.annotation.organization = 'Chalmers University of Technology';
model.annotation.note         = 'Rhodotorula toruloides - strain NP11';
model.id                      = 'rhto';
model.description             = 'Genome-scale metabolic model of Rhodotorula toruloides';

save([root '/scrap/model_r1.mat'],'model');
disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model); cd('reconstruction')
