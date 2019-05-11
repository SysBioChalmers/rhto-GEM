root  = regexprep(pwd(),'(.*)\\[^\\]*\\.*','$1');
scripts = [root '/ComplementaryScripts'];
data    = [root '/ComplementaryData'];
cd([scripts '/experimental/']);
%load([root '/scrap/model_r8.mat']);
model = importModel([root '\ModelFiles\xml\rhto.xml']);
expDat  = readLahtveeData();
%% Generate glucose model from Tiukova data
tiukovaData = readTiukovaData(model,2);
model = adjustRhtoBiomass(model,tiukovaData);
%% Readjust some reactions
% Ensure that exchange reactions for carbon sources are reversible
model.rev(getIndexes(model,{'r_1634','r_1718','r_1808'},'rxns')) = 1;
% Add protein excretion reation.
model = addExchangeRxns(model,'out','s_3717');
lipidPos        = find(strcmp(model.rxnNames,'lipid backbone pseudoreaction'));
scaleLp         = getIndexes(model,{'s_3733','s_1524','s_0666','s_0694'},'mets');

for i=1:length(expDat.sample)
    %% Scale lipids
    models{i} = model;
    %% Set exchange fluxes
    models{i} = setGAM(models{i},0);
    models{i} = setParam(models{i},'var',{'r_1634','r_1808','r_1718','r_1714'},...
        [-expDat.rates(i,2:4),0],5);
    fL = 1.01;
    while abs(fL-1) > 0
        fL0 = fL;
        models{i}.S(scaleLp,lipidPos) = full(models{i}.S(scaleLp,lipidPos))*fL;
        [~,~,~,~,~,L]   = sumBioMass(models{i});
        fL              = expDat.lipids(i)/L;
        disp(['Scaling lipid data: fL = ' num2str(fL0) ' -> difference = ' num2str(fL)])
    end
    %% Scale protein & other biomass components
    [models{i},GAMpol] = changeOtherComp(models{i},expDat,i);
    chainExIdx  = getIndexes(models{i},'r_4064','rxns');
    backbExIdx  = getIndexes(models{i},'r_4062','rxns');
    models{i} = setParam(models{i},'ub',[chainExIdx,backbExIdx],1000);
    %sol=solveLP(models{i});
    [models{i},k] = scaleLipidsRhto(models{i},'tails');
    models{i} = setParam(models{i},'eq',[chainExIdx,backbExIdx],0);
    models{i} = setParam(models{i},'lb','r_4046',0); % NGAM
    %% Protein excretion
    if ~(expDat.rates(i,9) == 0)
        [~,P]  = sumBioMass(models{i}); % Weight of protein pseudometabolite
        pIdx   = getIndexes(models{i},'EXC_OUT_s_3717','rxns');
        models{i}.S(find(models{1}.S(:,pIdx)),pIdx) = -1/P; % Make flux represent 1 g/gDCW/h
        models{i} = setParam(models{i},'var','EXC_OUT_s_3717',expDat.rates(i,9),5);
    end
    % Growth rate
    sol(1,i)=solveLP(models{i},1);
    disp(['Growth rate in sample ' num2str(i) ': ' num2str(-sol(1,i).f)]);
end
models2=models;
models=models2;
%% Set further constraints.
for i=1:numel(models)
    disp(['NGAM of condition ' expDat.sample{i} ' after fixing:']);
    % Set NGAM as objective function and track this while setting
    % constraints
    models{i} = setParam(models{i},'var','r_2111',expDat.rates(i,1),5); % Growth
    models{i} = setParam(models{i},'obj','r_4046',1);
    sols(i)=solveLP(models{i});
    one=num2str(-sols(i).f,4);
    models{i} = setParam(models{i},'var','r_2104',expDat.rates(i,7),5); % Xylitol
    sols(i)=solveLP(models{i},1);
    two=num2str(-sols(i).f,4);
    models{i} = setParam(models{i},'ub','r_1672',expDat.rates(i,5)); % CO2
    sols(i)=solveLP(models{i});
    three=num2str(-sols(i).f,4);
    models{i} = setParam(models{i},'lb','r_1992',expDat.rates(i,6)); % O2
    sols(i)=solveLP(models{i});
    four=num2str(-sols(i).f,4);
    models{i} = setParam(models{i},'ub','r_1654',-expDat.rates(i,8)); % Ammonium
    sols(i)=solveLP(models{i},1);
    five=num2str(-sols(i).f,4);    
    disp('Growth  Xylitol CTR     OTR     Ammonium');
    disp([one '   ' two '   ' three '   ' four '   ' five]);
end

%% Export models
for i=1:numel(models)
    fn = ['rhto_' expDat.sample{i} '.xml'];
    exportModel(models{i},[data '/Lahtvee2019/' fn]);
end

%% Summarize results
[~,idx] = getExchangeRxns(model,'both');
idx = [idx; getIndexes(model,'r_4046','rxns')]; % Include NGAM
clear out;for i=1:numel(models); out(:,i) = sols(i).x(idx); end

rmIdx = sum(out,2) == 0;
idx = idx(~rmIdx);
out = out(~rmIdx,:);

out = [model.rxnNames(idx) num2cell(out)];
fid = fopen([data '/Lahtvee2019/exp_exchangeRxns.tsv'],'w');
fprintf(fid,['%s' repmat('\t%s',1,12) '\n'],["reaction" string(expDat.sample')]);
for j=1:length(out)
    fprintf(fid,['%s' repmat('\t%d',1,12) '\n'],out{j,:});
end
fclose(fid);