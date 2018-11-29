%% Run MENECO
% MENECO requires the target compounds to already be part of the draft
% model. This should be fine here, as we added the whole biomass equation
% above.

% Find targets: any substrate for the pseudoreactions
rxnIdx=find(contains(model.rxnNames,'pseudoreaction'));
targets=find(any(model.S(:,rxnIdx)<0,2));
%targets=model.mets(targets);
%[model.mets(targets), model.metNames(targets)]
targetSBML=strcat('<species id="M_',model.mets(targets),...
    '" name="',model.metNames(targets),'"/>');

% Identified by MENECO (see meneco.out for output file).
% A minimum of 2 reactions are required, with different combinations of 3
% reactions. Not to favour one reaction over the other, as we don't
% have any prove at the moment which one is more likely to be present, we
% will add the union of reactions.
menecoRxns={'r_1665','r_2030','r_1667','r_1887'};

menecoRxns=getAllRxnsFromGenes(modelSce,menecoRxns);
model=addRxnsGenesMets(model,modelSce,menecoRxns,true,...
    'Identified by MENECO to produce biomass components',1); % Add reactions and metabolites

% Test meneco results
model_tmp=addExchangeRxns(model,'in',{'s_0394','s_0397','s_0458','s_0796','s_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452','s_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589','s_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602','s_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
%targets=model.mets(targets);
out=cell(length(targets),2);
for i=1:length(targets)
    tmpmodel=addExchangeRxns(model_tmp,'out',targets(i));
    tmpmodel=setParam(tmpmodel,'obj',numel(tmpmodel.rxns),1);
    out(i,1)=tmpmodel.metNames(getIndexes(tmpmodel,targets(i),'mets'));
    soltmp=solveLP(tmpmodel,1);
    out(i,2)=num2cell(soltmp.f);
end
out

%% Further gap-filling
% Almost all biomass components can be made, only lipids are problematic.
% This seemingly has to do with coenzyme A synthesis. As S. cerevisiae is
% known to produce coenzyme A, we'll see what reactions are required for
% this, and use this to curate the Rhto model.
tmpSce=addExchangeRxns(modelSce,'out','s_0529');
tmpSce=setParam(tmpSce,'obj',length(tmpSce.rxns),1);
sol=solveLP(tmpSce,1)

model=setParam(model,'eq','r_1714',-1);
model=setParam(model,'obj','r_2111',1);

sceFlux=tmpSce.rxns(~sol.x==0); % List of reactions that carry flux to make CoA in Sce.
sceFlux(numel(sceFlux))=[]; % Remove last entry, is exchange reaction
forLipids=sceFlux(~ismember(sceFlux,model.rxns));
forLipids=getAllRxnsFromGenes(modelSce,forLipids);

% check what reactiosn are required to make lipids in modelSce
model=addRxnsGenesMets(model,modelSce,forLipids,true,...
    'Identified to produce lipids',1); % Add reactions and metabolites
sol=solveLP(model)

%% Growth on xylose
% Through manual curation identified that r_1094 (xylulokinase) is missing
% for supporting growth on xylose.
model=addRxnsGenesMets(model,modelSce,'r_1094',true,...
    'Identified to assimilate xylose',1); % Add reactions and metabolites

save('../../scrap/model_r4.mat','model');
cd('..'); newCommit(model); cd('reconstruction')