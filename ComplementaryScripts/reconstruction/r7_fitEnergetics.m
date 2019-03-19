%% Fit to continuous glucose-limited chemostat cultivation, using data
% from Shen et al, J Biotech, 2013, doi:10.1016/j.jbiotec.2013.08.010,
% table 1.
load('../../scrap/model_r6.mat');

% Plot growth rate vs glucose uptake rate from Shen et al., 
fid         = fopen('../../ComplementaryData/validation/bioreactor_growth.csv');
fluxData    = textscan(fid,'%f32 %f32 %s %s','Delimiter',',','HeaderLines',1);
fluxData    = [fluxData{1} fluxData{2}];
fclose(fid);

b=polyfit(fluxData(1:16,1),fluxData(1:16,2),1);
% Offset at zero growth rate, this is the glucose requirement of NGAM
b(2);


% Determine how much ATP can be produced from glucose in the model
% Remove NADH-dependent succinate dehydrogenase
model = removeReactions(model,'r_4264',true,true,true);
% Set ICDH, THFS and all fatty-acid-CoA ligases as irreversible
model = setParam(model,'lb',{'r_0659','r_0446'},0);
model = setParam(model,'lb',contains(model.rxnNames,'fatty-acid--CoA ligase'),0);
model = setParam(model,'ub','r_4046',1000);
modelTmp = setParam(model,'obj','r_4046',1);
modelTmp = setParam(modelTmp,'ub','r_4046',1000);
modelTmp = setParam(modelTmp,'lb','r_1714',-1);
sol=solveLP(modelTmp,1)
% Confirm that nothing else is excreted
printFluxes(modelTmp,sol.x)

% Now, set the NGAM lower boundary to match NGAM
NGAM = -sol.f * b(2);
disp(['NGAM is set to: ' num2str(NGAM)])
model = setParam(model,'lb','r_4046',NGAM);

% Fit GAM
cd ../experimental
[model,GAM]=fitGAM(model);
disp(['GAM is set to: ' num2str(GAM)])

save('../../scrap/model_r7.mat','model');
cd('..'); newCommit(model); cd('reconstruction')
