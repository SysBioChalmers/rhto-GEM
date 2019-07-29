clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% This script predicts metabolic engineering targets for increased
% production of triglycerides and linolenic acid.
model = importModel('../../ModelFiles/xml/rhto.xml');

% Add exchange reactions for linolenate and triglyceride. For the
% triglyceride, we specifically choose 18:0/18:1/18:0-TAG (so called SOS)
% as target. In addition to functioning as representative of triglycerides
% for biofuel production, this SOS-TAG is a major part of cocoa-butter.
idx = getIndexes(model, {'linolenate[c]','triglyceride (1-18:0, 2-18:1, 3-18:0)[erm]'}, 'metscomps');

% Add exchange reactions for products
model = addExchangeRxns(model, 'out', idx);
% Keep track of ids of exchange reactions
rxn1 = model.rxns(end-1);
rxn2 = model.rxns(end);

%% Perform FSEOF for linolenic acid and TAG on glucose
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [-1, 0, 0]);

targets{1} = FSEOF(model, 'r_2111', rxn1, 10, 0.9);
targets{2} = FSEOF(model, 'r_2111', rxn2, 10, 0.9);

%% Perform FSEOF for linolenic acid and TAG on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, -1, 0]);

targets{3} = FSEOF(model, 'r_2111', rxn1, 10, 0.9);
targets{4} = FSEOF(model, 'r_2111', rxn2, 10, 0.9);

%% Perform FSEOF for linolenic acid and TAG on glycerol
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, 0, -1]);

targets{5} = FSEOF(model, 'r_2111', rxn1, 10, 0.9);
targets{6} = FSEOF(model, 'r_2111', rxn2, 10, 0.9);

%% Summarize results in table
geneAssoc = ~cellfun('isempty',model.grRules);
for i=1:size(targets,2)
    target(:,i)=targets{i}.logical;
end
for i=1:size(targets,2)
    slope(:,i)=targets{i}.slope;
    slope(~target(:,i),i)=nan;
end

target  = find(sum(target,2) & geneAssoc);
[~,I]=sort(sum(slope(target,:),2,'omitnan'),'descend');
out     = [num2cell(slope(target(I),:)), model.rxnNames(target(I)), model.grRules(target(I))];

fid = fopen([data '/results/fseof_TAG_linolenicAcid.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',["glu_LA" "glu_TAG" "xyl_LA" "xyl_TAG" ...
    "gly_LA" "gly_TAG" "rxnName" "grRule"]);
for j=1:length(I)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n',out{j,:});
end
fclose(fid);

%% R. toruloides produces a number of carotenoid pigments including
% torularhodin, torulene, b-carotene and y-carotene.
% https://www.ncbi.nlm.nih.gov/pubmed/17898860
clear target targets out slope
idx = getIndexes(model, 'torularhodin[c]', 'metscomps');

model = addExchangeRxns(model, 'out', idx);
% Keep track of ids of exchange reactions
rxn = model.rxns(end);

%% Perform FSEOF for carotenoids on glucose
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [-1, 0, 0]);
targets{1}  = FSEOF(model, 'r_2111', rxn, 10, 0.9);

%% Perform FSEOF for carotenoids on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, -1, 0]);
targets{2}  = FSEOF(model, 'r_2111', rxn, 10, 0.9);

%% Perform FSEOF for carotenoids on glycerol
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, 0, -1]);
targets{3}  = FSEOF(model, 'r_2111', rxn, 10, 0.9);

%% Summarize results in table
geneAssoc = ~cellfun('isempty',model.grRules);
for i=1:size(targets,2)
    target(:,i)=targets{i}.logical;
end
for i=1:size(targets,2)
    slope(:,i)=targets{i}.slope;
    slope(~target(:,i),i)=nan;
end

target  = find(sum(target,2) & geneAssoc);
[~,I]=sort(sum(slope(target,:),2,'omitnan'),'descend');
out     = [num2cell(slope(target(I),:)), model.rxnNames(target(I)), model.grRules(target(I))];

fid = fopen([data '/results/fseof_torularhodin.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["glucose" "xylose" "glycerol" "rxnName" "grRule"]);
for j=1:length(I)
    fprintf(fid,'%d\t%d\t%d\t%s\t%s\n',out{j,:});
end
fclose(fid);