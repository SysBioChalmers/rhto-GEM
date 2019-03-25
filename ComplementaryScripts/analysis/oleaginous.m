%% This script evaluates importance of reactions for oleaginous phenotype
% (or more precisely: TAG production).
model = importModel('../../ModelFiles/xml/rhto.xml');
backup=model;
model=backup;
% Add exchange reactions for triglyceride (18:0/18:1/18:0-TAG).
idx = getIndexes(model, {'triglyceride (1-18:0, 2-18:1, 3-18:0)[erm]', ...
    'triglyceride (1-18:0, 2-18:1, 3-18:0)[lp]'}, 'metscomps');
% Add exchange reactions for products
rxnsToAdd.rxns          = 'exch_TAG';
rxnsToAdd.mets          = model.mets(idx);
rxnsToAdd.stoichCoeffs  = {[-1, -1]}; 
rxnsToAdd.lb            = 0;
model = addRxns(model,rxnsToAdd);

model = setParam(model,'obj','exch_TAG',1);
model = setParam(model,'lb','r_4046',0);
sol=solveLP(model,1)
printFluxes(model,sol.x)

% Make sure that COBRA Toolbox version >3 is installed
initCobraToolbox()

[Glc.grRatio, Glc.grRateKO, Glc.grRateWT, Glc.hasEffect, Glc.delRxn, Glc.fluxSolution] = singleRxnDeletion(model,'FBA');

modelXyl = setParam(model,'lb',{'r_1714','r_1718'},[0,-1]);
[Xyl.grRatio, Xyl.grRateKO, Xyl.grRateWT, Xyl.hasEffect, Xyl.delRxn, Xyl.fluxSolution] = singleRxnDeletion(modelXyl,'FBA');

idx = find(Glc.grRatio < 0.90 | Xyl.grRatio < 0.90);

out = [num2cell(Glc.grRatio(idx)*100) num2cell(Xyl.grRatio(idx)*100) ...
    model.rxns(idx) model.rxnNames(idx) constructEquations(model,idx)];

fid = fopen('ComplementaryData/results/oleaginous.tsv','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["glucose" "xylose" "rxns" "rxnName" "eqn"]);
for j=1:length(idx)
    fprintf(fid,'%d\t%d\t%s\t%s\t%s\n',out{j,:});
end
fclose(fid)