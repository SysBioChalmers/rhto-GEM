%% This script generates a first draft of a genome-scale model (here: R.
% toruloides) using the getModelFromHomology function in RAVEN. It uses
% S. cerevisiae as a template model, 

% BLAST against Sce genome. Rhto protein fasta obtained from JGI, modified
% to have genes in format RHTO_00001. Sce protein fasta 
%blastedRhto=getBlast('rhto','../../ComplementaryData/genome/rhto_np11.faa',...
%    'sce','../../ComplementaryData/genome/sce_s288c.faa');
%save('../../ComplementaryData/genome/reconstruction/blastedRhto.mat','blastedRhto');
load('../../ComplementaryData/genome/blastedRhto.mat');

% Load S. cerevisiae model (corrected by Benjamin), downloaded from
% https://github.com/SysBioChalmers/yeast-GEM/raw/master/ModelFiles/xml/yeastGEM.xml
% (version 8.2.0).
modelSce=importModel('../../ComplementaryData/reconstruction/yeastGEM_820.xml',true);
modelSce.id='sce';
modelSce.mets=regexprep(modelSce.mets,'\[[a-z]+\]$','');
[modelSce.grRules, modelSce.rxnGeneMat]=standardizeGrRules(modelSce);

% Confirm that the model is functional, set objective to growth.
modelSce=setParam(modelSce,'obj','r_4041',1);
sol=solveLP(modelSce)
printFluxes(modelSce,sol.x)

%% Do some cosmetic adjustments on the Sce model, before generating Yli model
%modelSce.metNames=regexprep(modelSce.metNames,'>',''); % Character is not allowed
%modelSce.rxnNames=regexprep(modelSce.rxnNames,'>',''); % Character is not allowed
%modelSce.metNames=regexprep(modelSce.metNames,', cytoplasmic','');  % Comma is problematic

% Set boundaries between -1000 and 1000.
%modelSce.lb(isinf(-modelSce.lb))=-1000;
%modelSce.ub(isinf(modelSce.ub))=1000;
modelSce=setParam(modelSce,'eq','r_4046',0);
% Put all genes in nucleus
modelSce.geneComps=zeros(length(modelSce.genes),1);
modelSce.geneComps(1:end)=11;

% r_1896 ss exchange reaction, but was named transport reaction. Later on,
% all transport reactions are selected by their name, so names should be
% correct.
modelSce.rxnNames(find(strcmp(modelSce.rxns,'r_1896')))={'L-homoserine exchange'};

% Change reversibility fields to match boundaries. Prevents problems with
% MENECO.
for i=1:length(modelSce.rxns)
	if modelSce.lb(i)==0 && not(modelSce.ub(i)==0);
        modelSce.rev(i)=0;
    elseif modelSce.lb(i)<=0 && modelSce.ub(i)>=0;
        modelSce.rev(i)=1;
	end
end

% Confirm that model is stil functional
solveLP(modelSce)
printFluxes(modelSce,ans.x)
% Remove shortNames
modelSce=rmfield(modelSce,'geneShortNames');

%% Generate draft model, based on homology.
modelRhto=getModelFromHomology(modelSce,blastedRhto,'rhto',{},1,false,10^-20,150,35);

%% Make model based on MetaPhOrs data
orthologs=readtable('..\..\ComplementaryData\reconstruction\MetaPhOrs_convertedIDs.csv');
% Only keep orthologs if phylomeDB tree had more than 2 members and
% evidence is 1.
moreThanTwo=str2double(regexprep(orthologs.phylome,'.*/([0-9]+)','$1'))>2;
evidence=str2double(regexprep(orthologs.phylome(:),'\/[0-9]+',''))==1;
keep=find(evidence & moreThanTwo);
orthologList(:,1)=orthologs.SaccharomycesCerevisiae(keep);
orthologList(:,2)=orthologs.RhodosporidiumToruloides(keep);
blastStructure=makeFakeBlastStructure(orthologList,'sce','rhto');
modelRhtoOrth=getModelFromHomology(modelSce,blastStructure,'rhto');
modelRhtoOrth.id='ortho';

modelComb=mergeModels({modelRhto,modelRhtoOrth});
modelRhto=contractModel(modelComb);
%% Copy pseudoreactions
rxns=modelSce.rxns(contains(modelSce.rxnNames,'pseudoreaction'));
modelRhto=addRxnsGenesMets(modelRhto,modelSce,rxns,false,...
    'Modeling reaction',1); % Add reactions and metabolites

rxns=modelSce.rxns(contains(modelSce.rxnNames,'SLIME'));
modelRhto=addRxnsGenesMets(modelRhto,modelSce,rxns,false,...
    'SLIME reaction',1); % Add reactions and metabolites

%% Add all exchange rxns (these were not gene annotated, and therefore not added in
% draft). Might not require all exchange rxns, but easier to remove
% unconnected ones later.
rxns=getExchangeRxns(modelSce);

modelRhto=addRxnsGenesMets(modelRhto,modelSce,rxns,false,...
    'Modeling reaction',1); % Add reactions and metabolites

%% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx=find(cellfun(@isempty,modelSce.grRules)); % Which rxns have no genes
rxnIdx=find(getTransportRxns(modelSce));
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns=modelSce.rxns(rxnIdx); % Obtain reaction IDs
modelRhto=addRxnsGenesMets(modelRhto,modelSce,rxns,false,...
    'Modeling reaction required for intercellular transport, gene unknown',1); % Add reactions and metabolites

%% Add ATP synthase and ox phosp, problematic to gap-fill
rxns={'r_0226','r_0438','r_0439'};
modelRhto=addRxnsGenesMets(modelRhto,modelSce,rxns,false,...
    'Energy metabolism',1); % Add reactions and metabolites

%     
% %% Add biomass and lipid pseudoreaction
% % Currently, biomass equations and lipid pseudoreactions are taken from S.
% % cerevisiae. However, it is possible that the R. toruloides biomass and
% % lipids have different components (that are not part of S. cerevisiae, or
% % vice versa). This should ideally be improved first.
% 
% rxnIdx=regexp(modelSce.rxnNames,'biomass|pseudoreaction','all'); % Which reactions contain biomass or pseudoreaction in rxnNames
% rxnIdx=find(~cellfun('isempty',rxnIdx)); % Get reaction indexes
% rxns=modelSce.rxns(rxnIdx); % Get reaction IDs
% 
% modelRhto=addRxnsGenesMets(modelRhto,modelSce,rxns,false,...
%     'Modeling reaction','1'); % Add reactions and metabolites

modelRhto=setParam(modelRhto,'eq','r_1714',-1);
modelRhto=setParam(modelRhto,'obj','r_2111',1);

%% Add old reactions
modelRhtoOld=importModel('../../../rhto_old.xml');
withGenes=~cellfun(@isempty,modelRhtoOld.grRules);
rxns=modelRhtoOld.rxns(withGenes);
oldRxns=rxns(find(~ismember(rxns,modelRhto.rxns)));
% 
% oldRxns_r=oldRxns(strncmp(oldRxns,'r_',2));
% oldRxns_r=regexprep(oldRxns_r,'_EX.*','');
% modelRhto=addRxnsGenesMets(modelRhto,modelSce,oldRxns_r,false);

oldRxns_mo=oldRxns(strncmp(oldRxns,'mo_',3));
oldRxns_t=oldRxns(strncmp(oldRxns,'t_',2));
oldRxns_t(42:end)='';
% oldRxns_kegg=oldRxns(strncmp(oldRxns,'R',1));

modelRhtoOld.metNames(getIndexes(modelRhtoOld,'s_2877','mets'))={'palmitoleoyl-CoA(4-)'};
modelRhtoOld.metNames(getIndexes(modelRhtoOld,...
    {'s_0260';'s_0415';'s_0416';'s_1188'},'mets'))=...
    {'3-phosphonato-D-glycerate(3-)',...
    'alpha-D-ribose 1-phosphate(2-)';'alpha-D-ribose 1-phosphate(2-)',...
    'N-[(R)-4-phosphonopantothenoyl]-L-cysteine [cytoplasm]'};
modelRhto=addRxnsGenesMets(modelRhto,modelRhtoOld,oldRxns_mo,true);
modelRhto=addRxnsGenesMets(modelRhto,modelRhtoOld,oldRxns_t,true);

exportModel(modelRhto,'../../scrap/rhtoDraft.xml',true);
modelRhtoBck=modelRhto;
modelRhto=modelRhtoBck;
%% Run MENECO
% MENECO requires the target compounds to already be part of the draft
% model. This should be fine here, as we added the whole biomass equation
% above.

% Find targets: any substrate for the pseudoreactions
rxnIdx=find(contains(modelRhto.rxnNames,'pseudoreaction'));
targets=find(any(modelRhto.S(:,rxnIdx)<0,2));
%targets=modelRhto.mets(targets);
%[modelRhto.mets(targets), modelRhto.metNames(targets)]
targetSBML=strcat('<species id="M_',modelRhto.mets(targets),...
    '" name="',modelRhto.metNames(targets),'"/>');

% Identified by MENECO (see meneco_20161220.out for output file).
% A minimum of 62 reactions are required. There are different combinations
% of 62 reactions. Not to favour one reaction over the other, as we don't
% have any prove at the moment which one is more likely to be present, we
% will add the union of reactions.
%with all old reactions
menecoRxns={'r_0737','r_2119','r_0347','r_0678','r_0510','r_0851','r_0938','r_0736','r_0350','r_2030','r_2118','r_1027','r_1083','r_0904','r_0537','r_0346','r_0236','r_1667','r_1051','r_0883','r_0067','r_0109','r_0913','r_0237','r_0066','r_0738','r_2117','r_1063','r_0735','r_0942','r_1887','r_0195','r_0453','r_1665','r_0344'};

menecoRxns=getAllRxnsFromGenes(modelSce,menecoRxns);
modelRhto=addRxnsGenesMets(modelRhto,modelSce,menecoRxns,true,...
    'Identified by MENECO to produce biomass components',1); % Add reactions and metabolites

% Test meneco results
modelRhto_tmp=addExchangeRxns(modelRhto,'in',{'s_0394','s_0397','s_0458','s_0796','s_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452','s_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589','s_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602','s_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
targets=modelRhto.mets(targets);
out=cell(length(targets),2);
for i=1:length(targets)
    tmpmodel=addExchangeRxns(modelRhto_tmp,'out',targets(i));
    tmpmodel=setParam(tmpmodel,'obj',numel(tmpmodel.rxns),1);
    out(i,1)=tmpmodel.metNames(getIndexes(tmpmodel,targets(i),'mets'));
    soltmp=solveLP(tmpmodel,1);
    out(i,2)=num2cell(soltmp.f);
end
out

%exportModel(modelRhto,'rhto_afterFirstMeneco_20180728.xml');

tmpSce=addExchangeRxns(modelSce,'out','s_0529');
tmpSce=setParam(tmpSce,'obj',length(tmpSce.rxns),1);
sol=solveLP(tmpSce,1)


sceFlux=tmpSce.rxns(~sol.x==0); % List of reactions that carry flux to make CoA in Sce.
sceFlux(numel(sceFlux))=[]; % Remove last entry, is exchange reaction
forLipids=sceFlux(~ismember(sceFlux,modelRhto.rxns))

forLipids=getAllRxnsFromGenes(modelSce,forLipids)

% check what reactiosn are required to make lipids in modelSce
modelRhto=addRxnsGenesMets(modelRhto,modelSce,forLipids,true,...
    'Identified to produce lipids',1); % Add reactions and metabolites
sol=solveLP(modelRhto)
cd('..')
newCommit(modelRhto)

%% Rename all non-Yeast metIDs
nonS_xxxMet=~contains(modelRhto.mets,'s_');
newMetIds=cellfun(@(x) sprintf('st_%04s',x), ...
    string(1:length(find(nonS_xxxMet))), 'uni', false);
modelRhto.mets(nonS_xxxMet)=newMetIds;

nonS_xxxMet=~contains(modelRhto.mets,'s_');
newMetIds=cellfun(@(x) sprintf('st_%04s',x), ...
    string(1:length(find(nonS_xxxMet))), 'uni', false);
modelRhto.mets(nonS_xxxMet)=newMetIds;


%% Add more SLIMEr reactions
metsToAdd.mets={'MAGbb'};
metsToAdd.metNames={'monoglyceride backbone'};
metsToAdd.metFormulas={'C3H6O2'};
metsToAdd.compartments={'c'};
modelRhto=addMets(modelRhto,metsToAdd);

newCommit(modelRhto)

rxnsToAdd.rxns={'MAGbb1','MAGbb2','MAGbb3','MAGbb4'};
rxnsToAdd.rxnNames={'monoglyceride (16:0) [cytosol] SLIME rxn',...
    'monoglyceride (16:1) [cytosol] SLIME rxn',...
    'monoglyceride (18:0) [cytosol] SLIME rxn',...
    'monoglyceride (18:1) [cytosol] SLIME rxn'};
rxnsToAdd.equations={'1-monoglyceride (16:0)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C16:0 chain[c]',...
    '1-monoglyceride (16:1)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C16:1 chain[c]',...
    '1-monoglyceride (18:0)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C18:0 chain[c]',...
    '1-monoglyceride (18:1)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C18:1 chain[c]'};
modelRhto=addRxns(modelRhto,rxnsToAdd,3,'',false);
    
metsToAdd.mets={'st_0600','st_0601','st_0602'};
metsToAdd.metNames={'C14:0 chain','C18:2 chain','C18:3 chain'};
metsToAdd.metFormulas={'C14H28O2','C18H32O2','C18H30O2'};
metsToAdd.compartments={'c';'c';'c'};
modelRhto=addMets(modelRhto,metsToAdd);


%% Reactions to find gene associations for.
% getModelFromHomology left some OLD_sce genes that it could not find
% orthologs for. Additionally, the reactions that were added by MENECO and
% manually added for coenzyme A are still annotated with the Sce gene (not
% prefixed by OLD_sce_!) Try to find the responsible R. toruloides genes.

% All Sce genes have a Y in the name, while Rhto genes do not.
rxnIdx=strfind(modelRhto.grRules,'Y');
rxnIdx=~cellfun('isempty',rxnIdx); % Get reaction indexes

out=cell(length(find(rxnIdx)),3);
out(:,1)=modelRhto.rxns(rxnIdx);
out(:,2)=modelRhto.rxnNames(rxnIdx);
out(:,3)=modelRhto.grRules(rxnIdx);
out
save('modelRhto_20161220.mat','modelRhto');
load('modelRhto_20161220.mat')
%% Model from KEGG
modelKEGG=getKEGGModelForOrganism('uma','reRhoto1_prot.fasta','D:\KEGGdump\euk100_kegg80','D:\KEGGdump\rhto',false,false,false)
save('modelKEGG_20161220.mat','modelKEGG');