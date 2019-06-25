%% Curate lipid metabolism, by modifying reactions with metabolites that
% have an acyl-chain.
load('../../scrap/model_r4.mat');
load('../../scrap/modelTemplate.mat'); % yeast-GEM 8.2.0

%% Remove all SLIME reactions, will be manually added if needed
model = removeReactions(model,contains(model.rxnNames,'SLIME rxn'));

%% Curate reactions specified in lipidTemplates.txt
cd([scripts '/reconstruction'])
% Load text file that contains template reactions, in which compartments
% they are located, which grRules should be applied (specified for each
% compartment), and which sets of acyl-chains are used.
fid         = fopen ([data '/reconstruction/lipidTemplates.txt']);
loadedData  = textscan(fid,[repmat('%q ',[1,21]) '%q'],'delimiter','\t');
fclose(fid);
templNames = loadedData{1};   templEqns = loadedData{2};   comps = loadedData{3};
grRules    = loadedData{4};   chains = {};
for k = 1:length(loadedData)-4
    chains(:,k) = loadedData{k+4};
end

% Remove reactions that match a lipid template reaction (ignoring acyl-chains)
toRemove = regexprep(templNames,'CHAIN.*','');
% Specifically make sure that all MLCL reactions are removed
toRemove = [toRemove; 'MLCL'];
toRemove = startsWith(model.rxnNames,toRemove);
model = removeReactions(model,toRemove);

% Loop through all template reactions, and make new reactions for each set
% of acyl-chains, and for each specified compartment.
for i = 1:length(templNames)
    clear rxnsToAdd
    grRule  = strtrim(split(grRules{i},','));
    comp    = strtrim(split(comps{i},','));
    chain   = chains(i,:);
    % Replace compartments and acyl-chains in template reactions
    [newEqns, newNames] = makeLipidRxns(templEqns(i),templNames(i),chain,comp);
    % Repeat grRule for all reactions generated above. If no grRule is
    % specified, make empty matrix.
    if length(grRule)==1
        rxnsToAdd.grRules = cell(repmat(grRule,length(newEqns),1));
    else
        idx = 1:1:length(newEqns);
        idx = reshape(idx,[length(idx)/length(comp),length(comp)]);
        for j = 1:length(comp)
            rxnsToAdd.grRules(idx(:,j),1) = grRule(j);
        end
    end
    % Query reaction name to see if it already exists in the model. If so,
    % then don't add it.
    [Lia, Locb] = ismember(newNames, model.rxnNames);
    rxnsToAdd.grRules(Lia)= [];
    rxnsToAdd.equations = newEqns(~Lia);
    rxnsToAdd.rxnNames  = newNames(~Lia);
    rxnsToAdd.rxns      = generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
    % Query reaction name to see if it also existed in yeast-GEM. If so,
    % then use the reaction ID from yeast-GEM.
    [Lia,Locb]=ismember(rxnsToAdd.rxnNames, modelSce.rxnNames);
    if any(Locb)
        rxnsToAdd.rxns(Lia) = modelSce.rxns(Locb(Lia));
    end
    model = addRxns(model,rxnsToAdd,3,'','m_',true);
end

%% New lipid transport - similar as above, convert template reactions based
% on acyl-chains.
fid     = fopen([data '/reconstruction/lipidTransport.txt']);
loadedData    = textscan(fid,[repmat('%q ',[1,19]) '%q'],'delimiter','\t');
fclose(fid);
templNames = loadedData{1};   templEqns = loadedData{2}; chains = {};
for k = 1:length(loadedData)-2
    chains(:,k) = loadedData{k+2};
end

% Remove existing lipid transport reactions
rxnNames = regexprep(templNames,'CHAIN.*transport.*','');
rxnNames = startsWith(model.rxnNames,rxnNames);
model = removeReactions(model,rxnNames);

for i = 1:length(templNames)
    clear rxnsToAdd
    chain   = chains(i,:);
    [newEqns, newNames] = makeLipidRxns(templEqns(i),templNames(i),chain);
    [Lia, Locb] = ismember(newNames, model.rxnNames);
    rxnsToAdd.equations = newEqns(~Lia);
    rxnsToAdd.rxnNames  = newNames(~Lia);
    rxnsToAdd.rxns      = generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
    [Lia,Locb]=ismember(rxnsToAdd.rxnNames, modelSce.rxnNames);
    if any(Locb)
        rxnsToAdd.rxns(Lia) = modelSce.rxns(Locb(Lia));
    end
    model = addRxns(model,rxnsToAdd,3,'','m_',true);
end

%% Add SLIME reactions
cd ../experimental
% First gather which acyl-chains and lipid classes are present in the
% biomass measurements.
expData = readTiukovaData(model,1);
cd ../reconstruction

% Load SLIME template reactions.
fid     = fopen([data '/reconstruction/SLIMERtemplates.txt']);
loadedData    = textscan(fid,[repmat('%q ',[1,21]) '%q'],'delimiter','\t');
fclose(fid);
lipids = loadedData{2};	abbrev = loadedData{1};   comps = loadedData{3};
rxnNames = loadedData{4}; chains = {};
for k = 1:length(loadedData)-4
    chains(:,k) = loadedData{k+4};
end

for i = 1:length(lipids)
    clear rxnsToAdd
    % Find MW of backbone
    bbIdx = find(ismember(expData.lipidData.metAbbrev,abbrev(i)));
    bbMW = expData.lipidData.MW(bbIdx);
    
    % Go through list of chains per backbone
    chainList = chains(i,:);
    chainList = chainList(~cellfun('isempty',chainList));
    for j = 1:length(chainList)
        chain   = split(chainList(j),',');
        lipid   = lipids(i);
        for k = 1:length(chain)
            lipid = regexprep(lipid,['CHAIN' num2str(k)],chain{k});
        end
        [chainCount,~,n]    = unique(chain);
        [n,~]       = histc(n,unique(n));
        [~,coeff]   = ismember(chainCount,expData.chainData.chain);
        products    = [expData.lipidData.metIds(bbIdx); ...
                      expData.chainData.metIds(coeff)];
        totalMW     = sum([expData.lipidData.MW(bbIdx); ...
                      expData.chainData.MW(coeff).*n]);
        coeff       = [totalMW; expData.chainData.MW(coeff).*n]/1000;
        substrate   = model.mets(getIndexes(model,[lipid{1} ...
                      '[' comps{i} ']'], 'metscomps'));
        
        rxnsToAdd.mets          = [substrate; products];
        rxnsToAdd.stoichCoeffs  = [-1; coeff];
        rxnsToAdd.rxnNames      = {[lipid{1} ' ' rxnNames{i}]};
        rxnsToAdd.rxnMiriams    = struct('name',{{'sbo'}},'value',{{'SBO:0000395'}});
        rxnsToAdd.lb            = 0;
        rxnsToAdd.rxnConfidenceScores = 1;
        rxnsToAdd.rxns          = generateNewIds(model,'rxns','t_',1);
        [Lia,Locb]=ismember(rxnsToAdd.rxnNames, modelSce.rxnNames);
        if any(Locb)
            rxnsToAdd.rxns(Lia) = modelSce.rxns(Locb(Lia));
        end        
        model = addRxns(model,rxnsToAdd);
    end
end
% Add transport of free fatty acids between cytoplasm and ER membrane
model           = addTransport(model,'c','erm',expData.chainData.FA,true,false,'t_');
model           = addTransport(model,'c','ce',expData.chainData.FA,true,false,'t_');

% Add SLIME reactions for free fatty acids
clear rxnsToAdd
for i=1:length(expData.chainData.chain)
    mets            = model.mets(getIndexes(model,[expData.chainData.FA{i} '[c]'],'metscomps'));
    mets            = [mets, model.mets(getIndexes(model,'fatty acid backbone[c]','metscomps'))];
    mets            = [mets, model.mets(getIndexes(model,['C' expData.chainData.chain{i} ' chain[c]'],'metscomps'))];
    coeffs          = expData.chainData.MW(i);
    coeffs          = [-1, (coeffs+1.008)/1000, coeffs/1000];
    rxnsToAdd.mets  = mets;
    rxnsToAdd.stoichCoeffs = coeffs;
    rxnsToAdd.rxnNames     = strcat(expData.chainData.FA(i), ' SLIME rxn');
    rxnsToAdd.rxnMiriams   = struct('name',{{'sbo'}},'value',{{'SBO:0000395'}});
    rxnsToAdd.lb           = 0;
    rxnsToAdd.rxnConfidenceScores = 1;
    rxnsToAdd.rxns         = generateNewIds(model,'rxns','t_',1);
    model = addRxns(model,rxnsToAdd);
end

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model = setParam(model,'ub',[chainExIdx,backbExIdx],1000);
sol=solveLP(model,1)

cd ../experimental
[model,k] = adjustRhtoBiomass(model,expData);

save('../../scrap/model_r5.mat','model');
cd('..'); newCommit(model); cd('reconstruction')