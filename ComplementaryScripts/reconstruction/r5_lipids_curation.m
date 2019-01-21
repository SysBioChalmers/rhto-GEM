%% Remove reactions involving C16:1, this chain is only present at very low
% levels and is not considered to reduce the complexity of lipid metabolism
load('../../scrap/model_r4.mat');
% Remove reactions where one of the reactants contains 16:1 or similar in
% % % the metNames
% toRemove = find(contains(model.metNames,{'16:1','palmitoleate','palmitoleoyl'}));
% [row,col] = find(model.S(toRemove,:));
% col(ismember(model.rxns(col),'r_4065')) = []; % Keep lipid pseudoreaction
% model = removeReactions(model,col,true,true,true);

% % Remove all SLIME reactions, will be manually added if needed
% model = removeReactions(model,contains(model.rxnNames,'SLIME rxn'));

%% Add 18:3 synthesis and transport (degradation is missing for now)
metsToAdd.metNames={'linoleate';'linolenate';'linolenoyl-CoA'};
metsToAdd.compartments={'c';'c';'c'};
metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
%metsToAdd.metFormulas
model=addMets(model,metsToAdd); clear metsToAdd;

model=addTransport(model,'c',{'erm','lp','mm'},'linolenoyl-CoA',true,false);
model=addTransport(model,'c',{'erm','lp','mm'},'linoleoyl-CoA',true,false);

model=addTransport(model,'c',{'erm','lp','mm'},'linolenate',true,false);
model=addTransport(model,'c',{'erm','lp','mm'},'linoleate',true,false);

clear rxnsToAdd
rxnsToAdd.equations={'H+[erm] + oxygen[erm] + NADH[erm] + linoleoyl-CoA[erm] => 2 H2O[erm] + NAD[erm] + linolenoyl-CoA[erm]';...
    'coenzyme A[erm] + ATP[erm] + linolenate[erm] <=> AMP[erm] + diphosphate[erm] + linolenoyl-CoA[erm]';...
    'coenzyme A[lp] + ATP[lp] + linolenate[lp] <=> diphosphate[lp] + AMP[lp] + linolenoyl-CoA[lp]';...
    'ATP[p] + coenzyme A[p] + linolenate[p] <=> AMP[p] + diphosphate[p] + linolenoyl-CoA[p]';...
    'H2O[p] + linolenoyl-CoA[p] => coenzyme A[p] + H+[p] + linolenate[p]';...
    'ATP[c] + H2O[c] + linolenoyl-CoA[c] => ADP[c] + H+[c] + phosphate[c] + linolenoyl-CoA[p]'};
rxnsToAdd.rxnNames={'linoleoyl-CoA desaturase (n-C18:2CoA - n-C18:3CoA), ER membrane';...
    'fatty-acid--CoA ligase (octadecatrienoate), ER membrane';...
    'fatty-acid--CoA ligase (octadecatrienoate), lipid particle';...
    'fatty-acid--CoA ligase (octadecatrienoate), peroxisome';...
    'peroxisomal acyl-CoA thioesterase (18:3)';...
    'fatty acyl-CoA transport via ABC system (C18:3)'};
rxnsToAdd.rxns=generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));;
model=addRxns(model,rxnsToAdd,3,'',true,false); clear rxnsToAdd

%% Curate reactions specified in lipidTemplates.txt
fid     = fopen ('lipidTemplates.txt');
data    = textscan(fid,'%s %s %s %s %s %s %s %s','delimiter','\t');
fclose(fid);
templNames = data{1};   templEqns  = data{2};   chain1     = data{3};
chain2     = data{4};   chain3     = data{5};   chain4     = data{6};
comps      = data{7};   grRules    = data{8};
for i = 1:length(templNames)
    clear rxnsToAdd
    
    grRule  = strtrim(split(grRules{i},','));
    comp    = strtrim(split(comps{i},','));
    ch1     = strtrim(split(chain1{i},','));
    ch2     = strtrim(split(chain2{i},','));
    ch3     = strtrim(split(chain3{i},','));
    ch4     = strtrim(split(chain4{i},','));
    [newEqns, newNames] = makeRxns(templEqns(i),templNames(i),true,...
        comp,ch1,ch2,ch3,ch4);
    if length(grRule)==1
        rxnsToAdd.grRules = cell(repmat(grRule,length(newEqns),1));
    else
        idx = 1:1:length(newEqns);
        idx = reshape(idx,[length(idx)/length(comp),length(comp)]);
        for j = 1:length(comp)
            rxnsToAdd.grRules(idx(:,j),1) = grRule(j);
        end
    end
    [Lia, Locb] = ismember(newNames, model.rxnNames); % Find if reaction already exists
    rxnsToAdd.grRules(Lia)= [];
    rxnsToAdd.equations = newEqns(~Lia);
    rxnsToAdd.rxnNames  = newNames(~Lia);
    rxnsToAdd.rxns      = generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
    [Lia,Locb]=ismember(rxnsToAdd.rxnNames, modelSce.rxnNames);
    if any(Locb)
        rxnsToAdd.rxns(Lia) = modelSce.rxns(Locb(Lia));
    end
    model = addRxns(model,rxnsToAdd,3,'','st_',true);
end

%% Add lipid transport reactions
fid     = fopen ('lipidTransport.txt');
data    = textscan(fid,'%s %s %s %s %s','delimiter','\t');
fclose(fid);
templNames = data{1};   templEqns  = data{2};   chain1     = data{3};
chain2     = data{4};   chain3     = data{5};
for i = 1:length(templNames)
    clear rxnsToAdd
    ch1     = strtrim(split(chain1{i},','));
    ch2     = strtrim(split(chain2{i},','));
    ch3     = strtrim(split(chain3{i},','));
    [newEqns, newNames] = makeRxns(templEqns(i),templNames(i),true,...
        '',ch1,ch2,ch3);
    [Lia, Locb] = ismember(newNames, model.rxnNames); % Find if reaction already exists
    rxnsToAdd.equations = newEqns(~Lia);
    rxnsToAdd.rxnNames  = newNames(~Lia);
    rxnsToAdd.rxns      = generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
    [Lia,Locb]=ismember(rxnsToAdd.rxnNames, modelSce.rxnNames);
    if any(Locb)
        rxnsToAdd.rxns(Lia) = modelSce.rxns(Locb(Lia));
    end
    model = addRxns(model,rxnsToAdd,3,'','st_',true);
end

% Some manual reactions

%% Add SLIME metabolites
metsToAdd.metNames={'1-monoglyceride backbone','C18:2 chain','C18:3 chain'};
metsToAdd.metFormulas={'C3H6O2','C18H32O2','C18H30O2'};
metsToAdd.compartments={'lp';'c';'c'};
metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
model=addMets(model,metsToAdd); clear metsToAdd;

%% Add SLIME reactions
cd ../experimental
data = readTiukovaData(1);
cd ../reconstruction
clear rxnsToAdd
for i=1:length(data.lipidData.metNames);
    % Find compartment of component in lipid species part of lipid
    % pseudoreaction
    if ~strcmp('ergosterol',data.lipidData.metNames{i});
        compId = model.metComps(getIndexes(model,data.lipidData.metIDs(i),'mets'));
        rxnComp = model.comps(compId);
        % Find lipid metabolites in that compartment
        mets = find(startsWith(model.metNames,regexprep(data.lipidData.metNames(i),' backbone','')));
        mets(~(model.metComps(mets) == compId)) = [];
        mets(contains(model.metNames(mets),'backbone')) = [];
        % Loop through each metabolite
        for j = 1:length(mets)
            chains      = regexp(model.metNames(mets(j)),'(\d\d:\d)','match');
            [chains,~,n]= unique(chains{1});
            [n,~]       = histc(n,unique(n));
            [~,coeff]   = ismember(chains,data.chainData.chain);
            products    = [data.lipidData.metIDs(i); data.chainData.metIDs(coeff)];
            totalMW     = sum([data.lipidData.MW(i);data.chainData.MW(coeff).*n]);
            coeff       = [totalMW; data.chainData.MW(coeff).*n]/1000;
            rxnsToAdd.mets          = [model.mets(mets(j)); products];
            rxnsToAdd.stoichCoeffs  = [-1; coeff];
            rxnsToAdd.rxnNames      = strcat(model.metNames(mets(j)), ' SLIME rxn');
            rxnsToAdd.rxnMiriams    = struct('name',{{'sbo'}},'value',{{'SBO:0000395'}});
            rxnsToAdd.lb            = 0;
            rxnsToAdd.rxnConfidenceScores = 1;
            rxnsToAdd.rxns          = generateNewIds(model,'rxns','t_',1);
            model = addRxns(model,rxnsToAdd);
        end
    end
end


clear rxnsToAdd
for i=5:length(data.chainData.chain)
    % Make sure transport reactions exist
    mets            = getIndexes(model,data.chainData.FA{i},'metnames');
    metComps        = model.comps(model.metComps(mets));
    metComps(ismember(metComps,{'c','p','e'})) = [];
    model           = addTransport(model,'c',metComps,data.chainData.FA{i},true,false);
    
    % Add SLIME reaction from cytosolic species
    mets            = model.mets(getIndexes(model,[data.chainData.FA{i} '[c]'],'metscomps'));
    mets            = [mets, model.mets(getIndexes(model,'fatty acid backbone[c]','metscomps'))];
    mets            = [mets, model.mets(getIndexes(model,['C' data.chainData.chain{i} ' chain[c]'],'metscomps'))];
    coeffs          = data.chainData.MW(i);
    coeffs          = [-1, (coeffs+1.008)/1000, coeffs/1000];
    rxnsToAdd.mets  = mets;
    rxnsToAdd.stoichCoeffs = coeffs;
    rxnsToAdd.rxnNames     = strcat(data.chainData.FA(i), ' SLIME rxn');
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

save('../../scrap/model_r5.mat','model');
cd('..'); newCommit(model); cd('reconstruction')
