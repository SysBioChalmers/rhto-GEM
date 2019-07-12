clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Generate condition-specific models by modifying biomass composition and
% constraining with measured fluxes
model   = importModel([root '\ModelFiles\xml\rhto.xml']);

%% Readjust some reactions
% Ensure that exchange reactions for carbon sources are reversible
model = setParam(model,'rev',{'r_1634','r_1718','r_1808'},1);
model = setParam(model,'eq',{'r_1126','r_1127','r_1128'},0); % Fix citrate antiporter
model = setParam(model,'eq','y300057',0); %L-Glutamate 5-semialdehyde:NAD+ oxidoreductase
model = setParam(model,'eq','t_0082',0); %phosphate acyltransferase
model = setParam(model,'lb','r_4046',0); % NGAM set to zero, maximized later on

% Add protein excretion reation.
model           = addExchangeRxns(model,'out','s_3717');
lipidPos        = find(contains(model.rxnNames,{'lipid backbone pseudoreaction','lipid chain pseudoreaction'}));
[scaleLp,~]     = find(model.S(:,lipidPos)<0);
ngamIdx = getIndexes(model,'r_4046','rxns');
[~,allExIdx] = getExchangeRxns(model,'both');
allExIdx = [allExIdx; getIndexes(model,'r_4046','rxns')]; % Include NGAM

%% Load experimental data
cd([scripts '/experimental/'])
expDat      = readLahtveeData();

%% Make condition specific models, set rates and adjust biomass equation
clear out
for i=1:length(expDat.sample)
    models{i} = model;
    %% Set exchange fluxes
    models{i} = setGAM(models{i},0);
    % Set carbon source uptake, with 5% variance (95 - 105%).
    models{i} = setParam(models{i},'lb',{'r_1634','r_1808','r_1718','r_1714'},[-expDat.rates(i,2:4),0]);
    %% Scale lipids in biomass to fit measured total lipid content
    fL = 1.01;
    while abs(fL-1) > 0
        fL0 = fL;
        models{i}.S(scaleLp,lipidPos) = full(models{i}.S(scaleLp,lipidPos))*fL;
        [~,~,~,~,~,L]   = sumBioMass(models{i});
        fL              = expDat.lipids(i)/L;
    end
    %% Scale protein & other biomass components
    [models{i},GAMpol]  = changeOtherComp(models{i},expDat,i);
    %models{i}           = setGAM(models{i},GAMpol);
    %% Protein excretion
    if ~(expDat.rates(i,9) == 0)
        [~,P]  = sumBioMass(models{i}); % Weight of protein pseudometabolite
        pIdx   = getIndexes(models{i},'EXC_OUT_s_3717','rxns');
        models{i}.S(find(models{1}.S(:,pIdx)),pIdx) = -1/P; % Make flux represent 1 g/gDCW/h
        models{i} = setParam(models{i},'lb','EXC_OUT_s_3717',expDat.rates(i,9));
    end
    %% Xylitol excretion
    models{i} = setParam(models{i},'lb','r_2104',expDat.rates(i,7)); % Xylitol    
    models{i} = setParam(models{i},'ub','r_1672',expDat.rates(i,5)); % CO2
    models{i} = setParam(models{i},'lb','r_1992',expDat.rates(i,6)); % O2
    models{i} = setParam(models{i},'ub','r_1654',-expDat.rates(i,8)); % Ammonium    
    % Growth rate
    sol(1,i)=solveLP(models{i},1);
    gasIdx = getIndexes(models{i},{'r_1672','r_1992'},'rxns');
    atBound = find(sol(1,i).x(gasIdx) == transpose(expDat.rates(i,5:6)));
    if any(atBound)
       if atBound == 1
           models{i} = setParam(models{i},'ub','r_1672',1000);
       elseif atBound == 2
           models{i} = setParam(models{i},'lb','r_1992',-1000);
       end
    end
    sol(1,i)=solveLP(models{i},1);
    out(:,i)=sol(1,i).x;
    disp(['Growth rate in sample ' num2str(i) ': ' num2str(-sol(1,i).f)]);
end

%% Compare simulated growth with measured growth. Set the lower value as
% lower bound and then maximize NGAM
for i=1:length(models)
   sim = sol(1,i).x(getIndexes(models{i},'r_2111','rxns'));
   exp = expDat.rates(i,1);
   models{i} = setParam(models{i},'lb','r_2111',min([sim,exp]));
   models{i} = setParam(models{i},'obj','r_4046',1);
   sol(1,i) = solveLP(models{i},1);
   fba(:,i)=sol(1,i).x;
end

%% Random sampling
% Fix measured exchange fluxes around 5% of the value from FBA
exIdx = getIndexes(model,{'r_1634','r_1808','r_1718','r_1714',...
    'EXC_OUT_s_3717','r_2104','r_1672','r_1992','r_1654','r_2111','r_4046'},'rxns');

nsamples = 5000;%
for i=1:numel(models)
    % Fix measured exchange fluxes + NGAM around 5% of the value from FBA.
    fluxes = sol(1,i).x(exIdx);
    models{i} = setParam(models{i},'var',exIdx,fluxes,5);
%     if i==1
%         [~, goodRxns1] = randomSampling(models{i},1,true,true,true);
%     elseif i==5
%         [~, goodRxns5] = randomSampling(models{i},1,true,true,true);
%     elseif i==9
%         [~, goodRxns9] = randomSampling(models{i},1,true,true,true);
%     end
    if i<5
        rs{i} = randomSampling(models{i},nsamples,true,true,true,goodRxns1);
    elseif i<9
        rs{i} = randomSampling(models{i},nsamples,true,true,true,goodRxns5);
    else
        rs{i} = randomSampling(models{i},nsamples,true,true,true,goodRxns9);
    end
end

for i=1:length(rs)
    fluxMean(:,i) = full(mean(rs{i},2));
    fluxSD(:,i) = full(std(rs{i},0,2));
end

save([root '/scrap/randsampl.mat'], 'fba', 'rs', 'fluxMean', 'models', 'goodRxns*','fluxSD')


%% Write files
%FBA results, only exchange fluxes. For internal fluxes, refer to mean from
%random sampling.
out = out(allExIdx,:);
rmIdx = find(sum(cell2mat(out(:,3:end)),2) == 0);
out(rmIdx,:) =[];
fid = fopen([data '/Lahtvee2019/exp_pFBA_exchFluxes.tsv'],'w');
fprintf(fid,['%s' repmat('\t%s',1,13) '\n'],["reaction" "reactionName" string(expDat.sample')]);
for j=1:size(out,1)
    fprintf(fid,['%s\t%s' repmat('\t%f',1,12) '\n'],out{j,:});
end
fclose(fid);

%% Write file with all mean fluxes
clear out;for i=1:numel(models); out(:,i) = fluxMean(:,i); end

out = [models{1}.rxns models{1}.rxnNames num2cell(out)];
fid = fopen([data '/Lahtvee2019/exp_randSampl_allFluxes.tsv'],'w');
fprintf(fid,['%s' repmat('\t%s',1,13) '\n'],["rxnID" "reactionName" string(expDat.sample')]);
for j=1:length(out)
    fprintf(fid,['%s\t%s' repmat('\t%d',1,12) '\n'],out{j,:});
end
fclose(fid);

% Write exchange fluxes only
out = out(allExIdx,:);
rmIdx = find(sum(cell2mat(out(:,3:end)),2) == 0);
out(rmIdx,:) =[];
fid = fopen([data '/Lahtvee2019/exp_randSampl_exchFluxes.tsv'],'w');
fprintf(fid,['%s' repmat('\t%s',1,13) '\n'],["rxnID" "reactionName" string(expDat.sample')]);
for j=1:length(out)
    fprintf(fid,['%s\t%s' repmat('\t%d',1,12) '\n'],out{j,:});
end
fclose(fid);

%% NADH and NADPH reactions
for i={'NADPH','NADH'}
    [fluxes, rxnIdx] = getMetProduction(model,i,fluxMean,true);
    clear out
    out.rxns    = model.rxns(rxnIdx);
    out.rxnNames= model.rxnNames(rxnIdx);
    out.rxnEqns = constructEquations(model,rxnIdx);
    out.fluxes  = num2cell(fluxes);
    out = [out.rxns out.rxnNames out.rxnEqns out.fluxes];
    fid = fopen([data '/Lahtvee2019/exp_randSampl_' i{1} '_productionFluxes.tsv'],'w');
    fprintf(fid,['%s' repmat('\t%s',1,14) '\n'],["rxnID" "rxnName" ...
        "rxnEqn" string(expDat.sample')]);
    for j=1:length(out)
        fprintf(fid,['%s\t%s\t%s' repmat('\t%f',1,12) '\n'],out{j,:});
    end
    fclose(fid);
end

%% Export models
for i=1:numel(models)
    fn = ['rhto_' expDat.sample{i} '.xml'];
    exportModel(models{i},[data '/Lahtvee2019/' fn]);
end