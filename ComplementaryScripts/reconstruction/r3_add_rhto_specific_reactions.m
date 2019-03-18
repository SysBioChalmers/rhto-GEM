%% Rhodosporidium specific reactions
load('../../scrap/model_r2.mat','model');

% Colated a list of reactions from caretonoid metabolism, mitochondrial
% beta-oxidation and lipid metabolism related to C18:2

fid     = fopen ('../../ComplementaryData/reconstruction/rhtoSpecificMets.txt');
data    = textscan(fid,'%s %s %s %s','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

metsToAdd.metNames = data{1};
metsToAdd.compartments = data{2};
metsToAdd.metFormulas = data{3};
metsToAdd.mets = generateNewIds(model,'mets','m_',length(metsToAdd.metNames));
model=addMets(model,metsToAdd); clear metsToAdd;

fid     = fopen ('../../ComplementaryData/reconstruction/rhtoSpecificRxns.txt');
data    = textscan(fid,'%s %s %s %s %s','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

rxnsToAdd.equations = regexprep(data{1},'***','');
rxnsToAdd.rxnNames = regexprep(data{2},'***','');
rxnsToAdd.grRules = regexprep(data{3},'***','');
rxnsToAdd.subSystems = regexprep(data{4},'***','');
for i=1:numel(rxnsToAdd.subSystems)
    rxnsToAdd.subSystems{i}=rxnsToAdd.subSystems(i);
end
rxnsToAdd.eccodes = regexprep(data{5},'***','');
rxnsToAdd.rxns = generateNewIds(model,'rxns','t_',length(rxnsToAdd.rxnNames));
model=addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd

%model=deleteUnusedGenes(model);

save('../../scrap/model_r3.mat','model');
cd('..'); newCommit(model); cd('reconstruction')