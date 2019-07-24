clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% This script evaluates importance of reactions for oleaginous phenotype
% (or more precisely: TAG production).
model = importModel([root '/ModelFiles/xml/rhto.xml']);

%  Add exchange reactions for triglyceride (18:0/18:1/18:0-TAG).
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
initCobraToolbox(false)

model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [-1, 0, 0]);
Glc.grRatio = singleRxnDeletion(model,'FBA');

modelXyl    = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, -1, 0]);
Xyl.grRatio = singleRxnDeletion(modelXyl,'FBA');

modelGly    = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, 0, -1]);
Gly.grRatio = singleRxnDeletion(modelGly,'FBA');

idx = find(Glc.grRatio < 0.90 | Xyl.grRatio < 0.90 | Gly.grRatio < 0.90);

out = [num2cell(Glc.grRatio(idx)*100) num2cell(Xyl.grRatio(idx)*100) ...
    num2cell(Gly.grRatio(idx)*100) model.rxns(idx) model.rxnNames(idx) ...
    constructEquations(model,idx)];

fid = fopen([data '/results/oleaginous.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',["glucose" "xylose" "glycerol" "rxns" "rxnName" "eqn"]);
for j=1:length(idx)
    fprintf(fid,'%d\t%d\t%d\t%s\t%s\t%s\n',out{j,:});
end
fclose(fid);