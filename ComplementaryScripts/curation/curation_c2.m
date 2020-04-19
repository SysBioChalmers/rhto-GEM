clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
% This script takes rhto-GEM v1.2.1 and performs additional curation
model = importModel([root '/ModelFiles/xml/rhto.xml']);

%% Manual curate grRules to reduce invalid isoenzymes
fid         = fopen([data '/curation/curateGrRules.txt']);
loadedData  = textscan(fid,'%q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

rxns = loadedData{1};
grRules = loadedData{2};
model = changeGrRules(model,rxns,grRules);

%% Remove duplicate reactions
% Duplicate FA transport reactions
model = removeReactions(model,{'t_0871','t_0872','t_0873','t_0876','t_0877','t_0878'});

% Acetyl-CoA acyltransferase in mitochondria was duplicated.
model = setParam(model,'rev','t_0052',0);
model = setParam(model,'lb','t_0052',0);
% And incorrect ec number was assigned to mitochondrial 3-ketoacyl-CoA
% thiolases.
model.eccodes(getIndexes(model,{'t_0052','t_0053','t_0054','t_0055', ...
    't_0056','t_0057','t_0058','t_0059','t_0060','t_0061','t_0062', ...
    't_0063'},'rxns')) = {'2.3.1.16'};

save([root '/scrap/model_c2.mat'])
cd ..; newCommit(model); cd curation
