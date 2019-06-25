if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Manual curation
load([root '/scrap/model_r5.mat']);

% Set xylose exchange to reversible, to allow xylose uptake
model = setParam(model,'rev','r_1718',1);

% Scope of model is central metabolism. Protein modification reactions are
% not part of this
model = removeReactions(model,{'r_0518','r_0521','r_0523'});
model = removeReactions(model,{'r_4240','r_4243','r_4244','r_4271',...
    'r_4272','r_4320','r_4165','r_4166','r_4314'});
model = removeReactions(model,{'r_4252','r_4323','r_4324'});
model = removeReactions(model,{'r_0281','r_0282','r_0283','r_0284',...
    'r_0285','r_0286','r_0287','r_0288','r_0289','r_0290','r_0291',...
    'r_0292','r_0293','r_0294','r_0295','r_0296','r_0297','r_0298',...
    'r_0299'});

% Remove NADH-dependent succinate dehydrogenase
model = removeReactions(model,'r_4264',true,true,true);
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
%cd('..'); newCommit(model); cd('reconstruction')