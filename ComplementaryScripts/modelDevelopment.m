%% Load first model, and make uniform upper and lower boundaries
load('ModelFiles/mat/rhto.mat')
% Uniform upper and lower bounds
model.lb(model.lb == -99999) = -1000;
model.ub(model.ub == 99999) = 1000;
run ComplementaryScripts\saveModel.m


%% Update metabolite annotations
% From KEGG model
load('ModelFiles/mat/rhto.mat')
keggModel = getModelFromKEGG()

[tf, idx] = ismember(model.metNames, keggModel.metNames);
idx = idx(idx ~= 0);
model.metFormulas(tf) = keggModel.metFormulas(idx);
model.metMiriams = cell(length(model.mets),1);
model.metMiriams(tf) = transpose(keggModel.metMiriams(idx));
model.inchis = repmat({''},length(model.mets),1);
model.inchis(tf) = keggModel.inchis(idx);
model.inchis(cellfun(@isempty,model.inchis)) = {''};

% From Yeast model (load YeastMetabolicNetwork-GEM from SysBioChalmers)
sce = importModel('../YeastMetabolicNetwork-GEM/ModelFiles/xml/yeastGEM.xml')
[tf, idx] = ismember(model.metNames, sce.metNames);
idx = idx(idx ~= 0);
model.metFormulas(tf) = sce.metFormulas(idx);
model.metMiriams = cell(length(model.mets),1);
model.metMiriams(tf) = transpose(sce.metMiriams(idx));
model.metCharges = zeros(length(model.mets),1);
model.metCharges(tf) = sce.metCharges(idx);
model.metCharges(isnan(model.metCharges)) = 0;
run ComplementaryScripts\saveModel.m

%% Update rxn annotations
% From KEGG
load('ModelFiles/mat/rhto.mat')
[tf, idx] = ismember(model.rxns, keggModel.rxns);
idx = idx(idx ~= 0);
model.eccodes = repmat({''},length(model.rxns),1);
model.eccodes(tf) = keggModel.eccodes(idx);
model.rxnMiriams = cell(length(model.rxns),1);
model.rxnMiriams(tf) = keggModel.rxnMiriams(idx);
model.subSystems = cell(length(model.rxns),1);
model.subSystems(tf) = keggModel.subSystems(idx);
model.subSystems(cellfun(@isempty,model.subSystems)) = {''};

% From yeast
[tf, idx] = ismember(model.rxns, sce.rxns);
idx = idx(idx ~= 0);
model.rxnMiriams(tf) = sce.rxnMiriams(idx);
model.eccodes(tf) = sce.eccodes(idx);
%model.subSystems = cell(length(model.rxns),1);
%model.subSystems(tf) = sce.subSystems(idx);
model.eccodes(cellfun(@isempty,model.eccodes)) = {''};

run ComplementaryScripts\saveModel.m

%% Curation of duplicate reactions
load('ModelFiles/mat/rhto.mat')
% Duplicate reaction, add gene association to prefered reaction
model = removeReactions(model, {'R01600_EXP_2', 'R01786_EXP_2'},...
    true, true, true);
model = changeGeneAssoc(model, 'r_0534', '(RHTO_06072) or (RHTO_06870)');

model = removeReactions(model, 'R02737',...
    true, true, true);
model = changeGeneAssoc(model, 'r_0195', '(RHTO_01529) or (RHTO_06117)');

model = removeReactions(model, 'R00959_EXP_2',...
    true, true, true);
model = changeGeneAssoc(model, 'r_0888', '(RHTO_07820) or (RHTO_01766)');

model = removeReactions(model, {'R03656_EXP_2','R03648_EXP_2','R08218_EXP_2','R05577','R03650','R03656','R03905','R03659','R03654','R05578_EXP_2','R03664_EXP_2','R03658','R03662_EXP_2','R04212','R03664','R08218','R03940'},...
    true, true, true);
model = changeGeneAssoc(model, 'r_0665', '(RHTO_04405) or (RHTO_03176)');
model = changeGeneAssoc(model, 'r_0212', '(RHTO_04575) or (RHTO_02304)');
model = changeGeneAssoc(model, 'r_0995', '(RHTO_06886) or (RHTO_01083)');
model = changeGeneAssoc(model, 'r_0479', '(RHTO_03621) or (RHTO_08004)');
model = changeGeneAssoc(model, 'r_0729', '(RHTO_03139) or (RHTO_06168)');
model = changeGeneAssoc(model, 'r_1057', '(RHTO_08139) or (RHTO_03153)');
model = changeGeneAssoc(model, 'r_0995', '(RHTO_06886) or (RHTO_01083)');

model = removeReactions(model, 'R09246', true, true, true);
model = changeGeneAssoc(model, {'r_0281','r_0282','r_0283','r_0284',...
    'r_0285','r_0286','r_0287','r_0288','r_0289','r_0290','r_0291',...
    'r_0292','r_0293','r_0294','r_0295','r_0296','r_0297','r_0298',...
    'r_0299'}, '(RHTO_05050) or (RHTO_00464)');

run ComplementaryScripts\saveModel.m

%% Preferentially use yeast consensus network versions of reactions
load('ModelFiles/mat/rhto.mat')
% Reaction already exists as version from yeast consensus network.
% Preferably use that reaction instead.
model = addRxnsGenesMets(model, sce, 'r_0969', '(RHTO_07507)', 'Identified using KEGG HMM search', 2);
model = removeReactions(model, 'R01051',...
    true, true, true);
model = addRxnsGenesMets(model, sce, 'r_0544', '(RHTO_00071)', 'Identified using KEGG HMM search', 2);
model = removeReactions(model, 'R00650',...
    true, true, true);
model = addRxnsGenesMets(model, sce, 'r_0694', '(RHTO_07381)', 'Identified using KEGG HMM search', 2);
model = removeReactions(model, 'R00678',...
    true, true, true);

run ComplementaryScripts\saveModel.m

%% Strange localization suggested by KEGG, remove reactions
load('ModelFiles/mat/rhto.mat')
model = removeReactions(model, {'R01056', 'R01600', 'R01786'},...
    true, true, true);

% Duplicate / bad reactions (unspecific, outside of scope (=metabolism)).
model = removeReactions(model, {'R02396', 'R05071', 'R05068', 'R01262', ...
    'R02566','R02565','R02565_EXP_2','R02566_EXP_2', 'R01078', ...
    'R01005','R08368','R02657','R00100','R00100_EXP_2','R00021', ...
    'R03098','R03098_EXP_2','R02425','R07770','R04300','R02080', ...
    'R03105','R00896','R04996','R04996_EXP_2','R03980','R05731'},...
    true, true, true);

% Remove unconnected, isolated reactions (most are not gene annotated)
model = removeReactions(model, {'r_1581', 'r_1098', 'r_1792', ...
    'r_1164', 'r_1875', 'r_1878', 'r_1876', 'r_1877', 'r_2000', ...
    'r_1243', 'R02984', 'R03200', 'r_1554', 'r_1834', 'R01602', ...
    'r_1984', 'r_1776', 'r_2184', 'r_1228', 'r_1952'},...
    true, true, true);

run ComplementaryScripts\saveModel.m
