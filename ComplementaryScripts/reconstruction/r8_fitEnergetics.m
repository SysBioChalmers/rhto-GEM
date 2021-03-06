clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Fit to continuous glucose-limited chemostat cultivation, using data
% from Shen et al, J Biotech, 2013, doi:10.1016/j.jbiotec.2013.08.010,
% table 1.
load([root '/scrap/model_r7.mat']);

% Plot growth rate vs glucose uptake rate from Shen et al., 
fid         = fopen([data '/validation/bioreactor_growth.csv']);
fluxData    = textscan(fid,'%f32 %f32 %q %q','Delimiter',',','HeaderLines',1);
fluxData    = [fluxData{1} fluxData{2}];
fclose(fid);

b   = polyfit(fluxData(1:16,1),fluxData(1:16,2),1);
% Offset at zero growth rate, this is the glucose requirement of NGAM
b(2);

% Determine how much ATP can be produced from glucose in the model
modelTmp = setParam(model,'obj','r_4046',1);
modelTmp = setParam(modelTmp,'ub','r_4046',1000);
modelTmp = setParam(modelTmp,'lb','r_1714',-1);
sol      = solveLP(modelTmp,1)
% Confirm that nothing else is excreted
printFluxes(modelTmp,sol.x)

% Now, set the NGAM lower boundary to match NGAM
NGAM    = -sol.f * b(2);
disp(['NGAM is set to: ' num2str(NGAM)])
model   = setParam(model,'lb','r_4046',NGAM);

% Fit GAM
cd ../experimental
[model,GAM] = fitGAM(model);
disp(['GAM is set to: ' num2str(GAM)])

save([root '/scrap/model_r8.mat'],'model');
cd ../reconstruction
disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model); cd('reconstruction')
