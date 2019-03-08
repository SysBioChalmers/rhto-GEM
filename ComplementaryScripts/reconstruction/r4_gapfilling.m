%% Run MENECO
% MENECO requires the target compounds to already be part of the draft
% model. This should be fine here, as we added the whole biomass equation
% above.
load('../../scrap/model_r3.mat','model');
load('../../scrap/modelSce.mat');
exportModel(modelSce,'../meneco/sceRepair.xml')

%% Export model for meneco
exportModel(model,'../meneco/preMeneco.xml')

% Find targets: any substrate for the pseudoreactions, us the following
% text to reconstruct menecoTargets.sbml.
rxnIdx=find(contains(model.rxnNames,'pseudoreaction'));
targets=find(any(model.S(:,rxnIdx)<0,2));
[model.mets(targets), model.metNames(targets)]
targetSBML=strcat('<species id="M_',model.mets(targets),...
    '" name="',model.metNames(targets),'"/>');

% Identified by MENECO (see meneco.txt for output file).
% A minimum of 13 reactions are required, with different combinations of 19
% reactions. Not to favour one reaction over the other, as we don't
% have any prove at the moment which one is more likely to be present, we
% will add the union of reactions.
fid     = fopen ('../../ComplementaryData/meneco/menecoRxns.txt');
menecoRxns = textscan(fid,'%s'); fclose(fid);
menecoRxns = menecoRxns{1};

% If these reactions are present, that means that their respective enzymes
% are present. Any other reaction annotated to the same enzymes should also
% be added.
menecoRxns=getAllRxnsFromGenes(modelSce,menecoRxns);
model=addRxnsGenesMets(model,modelSce,menecoRxns,true,...
    'Identified by MENECO to produce biomass components',1); % Add reactions and metabolites

% Test meneco results, using the same exchange reactions as in menecoSeeds.sbml.
model_tmp=addExchangeRxns(model,'in',{'s_0394','s_0397','s_0458','s_0796',...
's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
out=cell(length(targets),2);
for i=1:length(targets)
    tmpmodel=addExchangeRxns(model_tmp,'out',targets(i));
    tmpmodel=setParam(tmpmodel,'obj',numel(tmpmodel.rxns),1);
    out(i,1)=tmpmodel.metNames(getIndexes(tmpmodel,targets(i),'mets'));
    soltmp=solveLP(tmpmodel,1);
    out(i,2)=num2cell(soltmp.f);
end
out
model=setParam(model,'obj','r_2111',1);
sol=solveLP(model,1)

%% Manual curation identified complex IV missing
model=addRxnsGenesMets(model,modelSce,'r_0438',true,...
    'Energy metabolism',1);

%% Growth on xylose
% Through manual curation identified that r_1094 (xylulokinase) is missing
% for supporting growth on xylose. Corresponding not known in R. toruloides
model=addRxnsGenesMets(model,modelSce,'r_1094',false,...
    'Identified to assimilate xylose',1); % Add reactions and metabolites

save('../../scrap/model_r4.mat','model');
cd('..'); newCommit(model); cd('reconstruction')