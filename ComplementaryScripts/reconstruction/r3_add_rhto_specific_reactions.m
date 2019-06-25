if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end

%% Rhodosporidium specific reactions
load([root,'/scrap/model_r2.mat'],'model');

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
%cd('..'); newCommit(model); cd('reconstruction')
