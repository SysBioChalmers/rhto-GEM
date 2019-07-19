clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Rhodotorula toruloides specific reactions
load([root,'/scrap/model_r2.mat'],'model');

% Manually lipid pseudometabolites, need specific metabolic ID for later scripts
metsToAdd.metNames      = {'1-monoglyceride backbone','C18:2 chain','C18:3 chain'};
metsToAdd.compartments  = {'lp','c','c'};
metsToAdd.metFormulas   = {'C3H6O2','',''};
metsToAdd.mets          = {'m_0100','m_0102','m_0103'};
metsToAdd.metMiriams    = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000649'}})},1,3);
model                   = addMets(model,metsToAdd); clear metsToAdd;

% Colated a list of reactions from caretonoid metabolism, mitochondrial
% beta-oxidation and lipid metabolism related to C18:2

fid         = fopen([data '/reconstruction/rhtoSpecificMets.txt']);
loadedData  = textscan(fid,'%q %q %q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

clear metsToAdd
metsToAdd.metNames      = loadedData{1};
metsToAdd.compartments  = loadedData{2};
metsToAdd.metFormulas   = loadedData{3};
metsToAdd.mets          = generateNewIds(model,'mets','m_',length(metsToAdd.metNames));
model                   = addMets(model,metsToAdd); clear metsToAdd;

fid         = fopen([data '/reconstruction/rhtoSpecificRxns.txt']);
loadedData  = textscan(fid,'%q %q %q %q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

clear rxnsToAdd
rxnsToAdd.equations     = regexprep(loadedData{1},'***','');
rxnsToAdd.rxnNames      = regexprep(loadedData{2},'***','');
rxnsToAdd.grRules       = regexprep(loadedData{3},'***','');
rxnsToAdd.subSystems    = regexprep(loadedData{4},'***','');
for i=1:numel(rxnsToAdd.subSystems)
    rxnsToAdd.subSystems{i}=rxnsToAdd.subSystems(i);
end
rxnsToAdd.eccodes       = regexprep(loadedData{5},'***','');
rxnsToAdd.rxns          = generateNewIds(model,'rxns','t_',length(rxnsToAdd.rxnNames));
model                   = addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd

save([root '/scrap/model_r3.mat'],'model');

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model); cd('reconstruction')
