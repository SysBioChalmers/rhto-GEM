clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
% This script takes rhto-GEM v1.2.0 and performs curation based on e.g.
% literature data
model = importModel([root '/ModelFiles/xml/rhto.xml']);

% Set default UB and LB
model.annotation.defaultUB = 1000;
model.annotation.defaultLB = -1000;
model.lb(isinf(model.lb)) = -1000;
model.ub(isinf(model.ub)) = 1000;

%% Curate glucose and xylose transporters
% From proteomics measurements of cultivations on glucose and xylose, infer
% which genes code for transporters (https://doi.org/10.1186/s13068-019-1478-8).
model = changeGrRules(model,'r_1166','RHTO_06080 or RHTO_01923 or RHTO_04266 or RHTO_00984 or RHTO_06016 or RHTO_03685 or RHTO_05157',true);
model = changeGrRules(model,'r_1719','RHTO_01630 or RHTO_03448 or RHTO_07444 or RHTO_06801 or RHTO_06080 or RHTO_01923',true);
rxnIDs = getIndexes(model,{'r_1166','r_1719'},'rxns');
model.rxnNotes(rxnIDs) = {'Differential proteomics expression'};
model.rxnReferences(rxnIDs) = {'doi:10.1186/s13068-019-1478-8'};

%% Curate gene association based on previous reconstruction
model = changeGrRules(model, 'r_0300', 'RHTO_07345', true);
model = changeGrRules(model, 'r_0301', 'RHTO_06406', true);
model = changeGrRules(model, 'r_0545', 'RHTO_06717', true);
model = changeGrRules(model, 'r_0658', 'RHTO_01289 or RHTO_01290', true);
model = changeGrRules(model, 'r_1094', 'RHTO_04556', true);
model = changeGrRules(model, 'r_0205', 'RHTO_00641 or RHTO_07387 ', true);
model = changeGrRules(model, 'r_0454', 'RHTO_05714', true);
model = changeGrRules(model, 'r_0455', 'RHTO_01560', true);
rxnIDs = getIndexes(model,{'r_0300','r_0301','r_0545','r_0658','r_1094','r_0205','r_0454','r_0455'},'rxns');
model.rxnNotes(rxnIDs) = {'Manual curation'};
model.rxnReferences(rxnIDs) = {'doi:10.1186/s12934-015-0217-5'};

% Two additional reactions were identified
metsToAdd.metNames      = {'L-xylulose'};
metsToAdd.compartments  = {'c'};
metsToAdd.metFormulas   = {'C5H10O5'};
metsToAdd.mets          = generateNewIds(model,'mets','m_',length(metsToAdd.metNames));
model                   = addMets(model,metsToAdd); clear metsToAdd;

rxnsToAdd.equations     = {'L-arabinitol[c] + NAD[c] <=> L-xylulose[c] + NADH[c] + H+[c]', 'L-xylulose[c] + NADPH[c] + H+[c] <=> xylitol[c] + NADP(+)[c]'};
rxnsToAdd.rxnNames      = {'arabinitol 4-dehydrogenase', 'L-xylulose reductase'};
rxnsToAdd.grRules       = {'RHTO_01629','RHTO_00373'};
rxnsToAdd.subSystems    = {{'Pentose and glucuronate interconversions'},{'Pentose and glucuronate interconversions'}};
rxnsToAdd.eccodes       = {'1.1.1.12','1.1.1.10'};
rxnsToAdd.rxns          = generateNewIds(model,'rxns','t_',length(rxnsToAdd.rxnNames));
rxnsToAdd.rxnNotes      = {'Manual curation', 'Manual curation'};
rxnsToAdd.rxnReferences = {'doi:10.1186/s12934-015-0217-5', 'doi:10.1186/s12934-015-0217-5'};
model                   = addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd

%% Various corrections
% Duplicate reactions
model = removeReactions(model,{'r_4235', 'y200008', 'r_1000', 'r_4262'}, true, true);

% Rephrased grRules, removing unclear gene associations
model = changeGrRules(model, 'r_0831', '(RHTO_02312 and RHTO_07893 and RHTO_04225) or (RHTO_02312 and RHTO_07893 and RHTO_07860)', true);
model = changeGrRules(model, 'r_0832', '(RHTO_02312 and RHTO_07893 and RHTO_04225) or (RHTO_02312 and RHTO_07893 and RHTO_07860)', true);

%%
save([root '/scrap/model_c1.mat'])
cd ..; newCommit(model); cd curation
