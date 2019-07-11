clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Curate lipid metabolism, by modifying reactions with metabolites that
% have an acyl-chain.
load([root '/scrap/model_r4.mat']);
load([root '/scrap/modelTemplate.mat']); % yeast-GEM 8.2.0

%% Curate reactions specified in lipidTemplates.txt
cd([scripts '/reconstruction'])
% Load text file that contains template reactions, in which compartments
% they are located, which grRules should be applied (specified for each
% compartment), and which sets of acyl-chains are used.

% Acyl chains with abundance >5%: 16:0, 18:0, 18:1, 18:2, 18.3. Assume sn1
% position to be satured (16:0 and 18:0). Assume sn2 position to be
% unsatured (18:1, 18:2, 18:3). To reduce complexity of cardiolipins, represent
% each set of acyl chains once (while CL(16:0, 18:1, 16:0, 18:2) is included,
% note sn2 and sn4 positions, CL(16:0, 18:2, 16:0, 18:1) is not in the model).
% During remodeling, sn1 position is first replaced by 18:1, donated by
% PC(16:0, 18:1), followed by replace sn3 position using same PC.

fid         = fopen([data '/reconstruction/lipidTemplates.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', '\t');
fclose(fid);

% Reorganize the content so that it can be used by the addLipidReactions
% function.
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.grRules     = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');
% Remove reactions that match a lipid template reaction (ignoring acyl-chains)
toRemove    = regexprep(template.rxns,'CHAIN.*','');
toRemove    = find(startsWith(model.rxnNames,toRemove));
model       = removeReactions(model,toRemove);

% Now use the templates to add the relevant reactions to the model. If a
% reaction already existed in the S. cerevisiae template model, then it
% will use the same reaction identifier.
model = addLipidReactions(template,model,modelSce);

fid         = fopen([data '/reconstruction/lipidTransport.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', '\t');
fclose(fid);

clear template
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.chains = {};
for k = 1:numCols-2; template.chains(:,k) = loadedData{k+3}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');

model = addLipidReactions(template,model,modelSce);

%SLIMER SLIMER SLIMER
clear template
% First remove any SLIME reactions that might exist in the draft model.
model = removeReactions(model,contains(model.rxnNames,'SLIME rxn'));

% Load SLIME template reactions and parse it through the addSLIMEreactions
% function to amend the model.
fid             = fopen([data '/reconstruction/SLIMERtemplates.txt']);
firstLine       = fgets(fid);
numCols         = numel(strfind(firstLine,char(9))); % number of \t
loadedData      = textscan(fid,['%q %q %f' repmat(' %q',[1,numCols-2])],'delimiter','\t');
fclose(fid);
template.metName    = loadedData{1};	template.bbID   = loadedData{2};
template.bbMW       = loadedData{3};    template.comps  = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');

model=addSLIMEreactions(template,model,modelSce);

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model       = setParam(model,'ub',[chainExIdx,backbExIdx],1000);
solveLP(model,1)

model           = addTransport(model,'c','erm',expData.chainData.FA,true,false,'t_');
model           = addTransport(model,'c','ce',expData.chainData.FA,true,false,'t_');

cd ../experimental
expData = readTiukovaData(model,1);
[model,k]   = adjustRhtoBiomass(model,expData);
sol=solveLP(model,1)

save([root '/scrap/model_r5.mat'],'model');
%cd('..'); newCommit(model); cd('reconstruction')
