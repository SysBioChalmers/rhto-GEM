%% Remove reactions involving C16:1, this chain is only present at very low
% levels and is not considered to reduce the complexity of lipid metabolism
load('../../scrap/model_r4.mat');
% Remove reactions where one of the reactants contains 16:1 or similar in
% the metNames
toRemove = find(contains(model.metNames,{'16:1','palmitoleate','palmitoleoyl'}));
[row,col] = find(model.S(toRemove,:));
col(ismember(model.rxns(col),'r_4065')) = []; % Keep lipid pseudoreaction
model = removeReactions(model,col,true,true,true);

%% Add 18:3 synthesis and transport (degradation is missing for now)
metsToAdd.metNames={'linolenate';'linolenate';'linolenate';'linolenoyl-CoA';'linolenoyl-CoA';'linolenoyl-CoA';'linolenoyl-CoA';'linoleoyl-CoA';'linoleoyl-CoA'};
metsToAdd.compartments={'erm';'lp';'p';'c';'erm';'lp';'p';'lp';'mm'};
metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
%metsToAdd.metFormulas
model=addMets(model,metsToAdd); clear metsToAdd;

clear rxnsToAdd
rxnsToAdd.equations={'H+[erm] + oxygen[erm] + NADH[erm] + linoleoyl-CoA[erm] => 2 H2O[erm] + NAD[erm] + linolenoyl-CoA[erm]';...
    'coenzyme A[erm] + ATP[erm] + linolenate[erm] <=> AMP[erm] + diphosphate[erm] + linolenoyl-CoA[erm]';...
    'coenzyme A[lp] + ATP[lp] + linolenate[lp] <=> diphosphate[lp] + AMP[lp] + linolenoyl-CoA[lp]';...
    'ATP[p] + coenzyme A[p] + linolenate[p] <=> AMP[p] + diphosphate[p] + linolenoyl-CoA[p]';...
    'H2O[p] + linolenoyl-CoA[p] => coenzyme A[p] + H+[p] + linolenate[p]';...
	'ATP[c] + H2O[c] + linolenoyl-CoA[c] => ADP[c] + H+[c] + phosphate[c] + linolenoyl-CoA[p]';...
    'linoleoyl-CoA[c] <=> linoleoyl-CoA[erm]';...
    'linoleoyl-CoA[c] <=> linoleoyl-CoA[lp]';...
    'linoleoyl-CoA[c] <=> linoleoyl-CoA[mm]'};
rxnsToAdd.rxnNames={'linoleoyl-CoA desaturase (n-C18:2CoA - n-C18:3CoA), ER membrane';...
    'fatty-acid--CoA ligase (octadecatrienoate), ER membrane';...
    'fatty-acid--CoA ligase (octadecatrienoate), lipid particle';...
    'fatty-acid--CoA ligase (octadecatrienoate), peroxisome';...
    'peroxisomal acyl-CoA thioesterase (18:3)';...
    'fatty acyl-CoA transport via ABC system (C18:3)';...
    'linolenoyl-CoA transport, cytoplasm-ER membrane';...
    'linolenoyl-CoA transport, cytoplasm-lipid particle';...
    'linolenoyl-CoA transport, cytoplasm-mitochondrial membrane'};
rxnsToAdd.rxns=generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));;
model=addRxns(model,rxnsToAdd,3,'',false,false); clear rxnsToAdd
% 
% %% Add more SLIMEr reactions
% metsToAdd.metNames={'monoglyceride backbone','C14:0 chain','C18:2 chain','C18:3 chain'};
% metsToAdd.metFormulas={'C3H6O2','C14H28O2','C18H32O2','C18H30O2'};
% metsToAdd.compartments={'c';'c';'c';'c'};
% metsToAdd.mets=generateNewIds(model,'mets','st_',length(metsToAdd.metNames));
% model=addMets(model,metsToAdd); clear metsToAdd;
% 
% rxnsToAdd.rxnNames={'monoglyceride (16:0) [cytosol] SLIME rxn',...
%     'monoglyceride (16:1) [cytosol] SLIME rxn',...
%     'monoglyceride (18:0) [cytosol] SLIME rxn',...
%     'monoglyceride (18:1) [cytosol] SLIME rxn'};
% rxnsToAdd.equations={'1-monoglyceride (16:0)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C16:0 chain[c]',...
%     '1-monoglyceride (16:1)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C16:1 chain[c]',...
%     '1-monoglyceride (18:0)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C18:0 chain[c]',...
%     '1-monoglyceride (18:1)[c] => 0.8574 monoglyceride backbone[c] + 0.25441 C18:1 chain[c]'};
% rxnsToAdd.lb=[0,0,0,0];
% rxnsToAdd.ub=[1000,1000,1000,1000];
% rxnsToAdd.rxns=generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
% model=addRxns(model,rxnsToAdd,3,'',false);
% 
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

%% Reactions to find gene associations for.
% % getModelFromHomology left some OLD_sce genes that it could not find
% % orthologs for. Additionally, the reactions that were added by MENECO and
% % manually added for coenzyme A are still annotated with the Sce gene (not
% % prefixed by OLD_sce_!) Try to find the responsible R. toruloides genes.
% 
% % All Sce genes have a Y in the name, while Rhto genes do not.
% rxnIdx=strfind(model.grRules,'Y');
% rxnIdx=~cellfun('isempty',rxnIdx); % Get reaction indexes
% 
% out=cell(length(find(rxnIdx)),3);
% out(:,1)=model.rxns(rxnIdx);
% out(:,2)=model.rxnNames(rxnIdx);
% out(:,3)=model.grRules(rxnIdx);
% out
% save('model_20161220.mat','model');
% load('model_20161220.mat')
% %% Model from KEGG
% modelKEGG=getKEGGModelForOrganism('uma','reRhoto1_prot.fasta','D:\KEGGdump\euk100_kegg80','D:\KEGGdump\rhto',false,false,false)
% save('modelKEGG_20161220.mat','modelKEGG');

save('../../scrap/model_r5.mat','model');
cd('..'); newCommit(model); cd('reconstruction')
