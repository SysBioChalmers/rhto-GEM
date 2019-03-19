%% Gene essentiality prediction based on T-DNA random insertion mutant
% library, from Coradetti et al (2018) eLife. doi:10.7554/eLife.32110
model       = importModel('../../ModelFiles/xml/rhto.xml');

% Make sure that COBRA Toolbox version >3 is installed
cd('../../../cobratoolbox')
initCobraToolbox()
modelCobra = ravenCobraWrapper(model); 

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(modelCobra, 'FBA');

%% Load experimental data
fid         = fopen('../rhtoGEM/ComplementaryData/validation/essential_Coradetti2018.csv');
essData     = textscan(fid,'%s %s','Delimiter',',','HeaderLines',1);
genes       = essData{1};
essential   = essData{2};
fclose(fid); clear essData

%% Some stats
gOverlap = length(find(ismember(genes,model.genes)));
gModel   = length(model.genes);
gExp     = length(genes);
disp(['Number of genes in model:  ' num2str(gModel)])
disp(['Number of genes with data: ' num2str(gExp)])
disp(['Model genes with data:     ' num2str(gOverlap) ' (' num2str((gOverlap/gModel)*100) '%)'])

%%
exp_ess = model.genes(le(grRatio,0.33) | isnan(grRatio));
exp_non = model.genes(grRatio>0.33);
dat_ess = genes(ismember(essential,'Yes'));
dat_non = genes(ismember(essential,'No'));

TP = length(intersect(exp_ess,dat_ess));
TN = length(intersect(exp_non,dat_non));
FP = length(intersect(exp_ess,dat_non));
FN = length(intersect(exp_non,dat_ess));

disp(sprintf(['TP: ' num2str(TP) '\nTN: ' num2str(TN) '\nFP: ' num2str(FP)...
    '\nFN: ' num2str(FN) '\nAccuracy: ' num2str(TN/(TN+FP))]));
