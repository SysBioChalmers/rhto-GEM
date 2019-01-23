%% This script generates a first draft of a genome-scale model (here: R.
% toruloides) using the getModelFromHomology function in RAVEN. It uses
% S. cerevisiae as a template model, 

% BLAST against Sce genome. Rhto protein fasta obtained from JGI, modified
% to have genes in format RHTO_00001. Sce protein fasta
%blastedRhto=getBlast('rhto','../../ComplementaryData/genome/rhto_np11.faa',...
%    'sce','../../ComplementaryData/genome/sce_s288c.faa');

%save('../../scrap/blastedRhto.mat','blastedRhto');
load('../../scrap/blastedRhto.mat');

% Load S. cerevisiae model (corrected by Benjamin), downloaded from
% https://github.com/SysBioChalmers/yeast-GEM/raw/master/ModelFiles/xml/yeastGEM.xml
% (version 8.2.0).
modelSce=importModel('../../ComplementaryData/reconstruction/yeastGEM_820.xml',true);
modelSce.id='sce';
modelSce.mets=regexprep(modelSce.mets,'\[[a-z]+\]$','');
[modelSce.grRules, modelSce.rxnGeneMat]=standardizeGrRules(modelSce);
% Put all genes in nucleus
modelSce.geneComps=zeros(length(modelSce.genes),1);
modelSce.geneComps(1:end)=11;
% r_1896 ss exchange reaction, but was named transport reaction. Later on,
% all transport reactions are selected by their name, so names should be
% correct.
modelSce.rxnNames(find(strcmp(modelSce.rxns,'r_1896')))={'L-homoserine exchange'};
% Remove shortNames
modelSce=rmfield(modelSce,'geneShortNames');

% Change reversibility fields to match boundaries. Prevents problems with
% MENECO.
for i=1:length(modelSce.rxns)
	if modelSce.lb(i)==0 && not(modelSce.ub(i)==0);
        modelSce.rev(i)=0;
    elseif modelSce.lb(i)<=0 && modelSce.ub(i)>=0;
        modelSce.rev(i)=1;
	end
end

% Confirm that the model is functional, set objective to growth.
modelSce=setParam(modelSce,'obj','r_4041',1);
sol=solveLP(modelSce)
printFluxes(modelSce,sol.x)
%save('../../scrap/modelSce.mat', 'modelSce');

%% Generate draft model, based on homology.
model=getModelFromHomology(modelSce,blastedRhto,'rhto',{},1,false,10^-20,150,35);

%% Make model based on MetaPhOrs data
orthologs=readtable('..\..\ComplementaryData\reconstruction\MetaPhOrs_convertedIDs.csv');
% Only keep orthologs if phylomeDB tree had more than 2 members and
% evidence is 1.
moreThanTwo=str2double(regexprep(orthologs.phylome,'.*/([0-9]+)','$1'))>1;
evidence=str2double(regexprep(orthologs.phylome(:),'\/[0-9]+',''))>=0.5;
keep=find(evidence & moreThanTwo);
orthologList(:,1)=orthologs.SaccharomycesCerevisiae(keep);
orthologList(:,2)=orthologs.RhodosporidiumToruloides(keep);
blastStructure=makeFakeBlastStructure(orthologList,'sce','rhto');
modelOrth=getModelFromHomology(modelSce,blastStructure,'rhto');
modelOrth.id='ortho';

modelComb=mergeModels({model,modelOrth});
model=contractModel(modelComb);

clear moreThanTwo evidence keep orthologList modelComb modelOrth orthologs

%% Add R. toruloides GEM meta data
model.annotation.defaultLB    = -1000;
model.annotation.defaultUB    = +1000;
model.annotation.taxonomy     = 'taxonomy/1130832';
model.annotation.givenName    = 'Eduard';
model.annotation.familyName   = 'Kerkhoven';
model.annotation.email        = 'eduardk@chalmers.se';
model.annotation.organization = 'Chalmers University of Technology';
model.annotation.note         = 'Rhodosporidium toruloides - strain NP11';
model.id                      = 'rhto';
model.description             = 'Genome-scale metabolic model of Rhodosporidium toruloides';

save('../../scrap/model_r1.mat','model');
cd('..'); newCommit(model); cd('reconstruction')