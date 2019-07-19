clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Manual curation
load([root '/scrap/model_r5.mat']);

% Set exchange reactions to alternative carbon sources to reversible
model = setParam(model,'rev',{'r_1634','r_1718','r_1808'},1);

% Remove NADH-dependent succinate dehydrogenase, uses CoQ instead.
model = removeReactions(model,'r_4264',true,true,true);
% Remove L-Glutamate 5-semialdehyde:NAD+ oxidoreductase, already
% represented with r_0012 and r_1887
model = removeReactions(model,'y300057',true,true,true);
% Remove citrate antiporters, curated in oleaginous yeast: 10.1016/j.ymben.2019.05.002
model = removeReactions(model,{'r_1126','r_1127','r_1128'},true,true,true);

% Set ICDH, THFS and all fatty-acid-CoA ligases as irreversible
model = setParam(model,'lb',{'r_0659','r_0446'},0);
model = setParam(model,'lb',contains(model.rxnNames,'fatty-acid--CoA ligase'),0);
model = setParam(model,'ub','r_4046',1000);

% % Mitochondrial tRNA synthetases: GEM doesn't model mitochondrial translation
% model = removeReactions(model,{'r_0210','r_0213','r_0540','r_0666',...
%     'r_0712','r_0730','r_0853','r_1043','r_1058','r_1067','r_1090',...
%     'r_0480','r_4155'},true,true,true);

%% Remove 'sce' from subsystems
model.subSystems = cellfun(@(x) regexprep(x,'sce[0-9]+ +',''),model.subSystems, 'UniformOutput', 0);

%% Remove unused metabolites
model = removeMets(model,all(model.S == 0,2),false,true,true,true);

%% Remove unconnected non-gene associated reactions
subGraphs = getAllSubGraphs(model);

% Find which reactions have no gene associated
rxnToRemove     = [];
for i=2:size(subGraphs,2)
    metIdx      = subGraphs(:,i);
    rxnIdx      = model.S(metIdx,:);
    [~,col,~]   = find(rxnIdx);
    col         = unique(col);
    grRules     = model.grRules(col);
    if isempty(grRules{1})
         rxnToRemove = [rxnToRemove; col];
    end
end

model = removeReactions(model,rxnToRemove,true,true,true);

save([root '/scrap/model_r6.mat'],'model');

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model); cd('reconstruction')
