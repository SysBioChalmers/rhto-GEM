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
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [-1, 0]);

targets = FSEOF(model, 'r_2111', rxn1, 10, 0.9, 'fseof_linolenic acid_glc.tab');
targets = FSEOF(model, 'r_2111', rxn2, 10, 0.9, 'fseof_TAG_glc.tab');

%% Perform FSEOF for linolenic acid and TAG on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [0, -1]);

targets = FSEOF(model, 'r_2111', rxn1, 10, 0.9, 'fseof_linolenic acid_xyl.tab');
targets = FSEOF(model, 'r_2111', rxn2, 10, 0.9, 'fseof_TAG_xyl.tab');

%% R. toruloides produces a number of carotenoid pigments including
% torularhodin, torulene, b-carotene and y-carotene.
% https://www.ncbi.nlm.nih.gov/pubmed/17898860

idx = getIndexes(model, {'torularhodin[c]', 'torulene[c]', ...
    'beta-carotene[c]', 'gamma-carotene'}, 'metscomps');

model = addExchangeRxns(model, 'out', idx);
% Keep track of ids of exchange reactions
rxn1 = model.rxns(end-3);
rxn2 = model.rxns(end-2);
rxn3 = model.rxns(end-1);
rxn4 = model.rxns(end);

%% Perform FSEOF for linolenic acid and TAG on glucose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [-1, 0]);

targets = FSEOF(model, 'r_2111', rxn1, 10, 0.9, 'fseof_torularhodin_glc.tab');
targets = FSEOF(model, 'r_2111', rxn2, 10, 0.9, 'fseof_torulene_glc.tab');
targets = FSEOF(model, 'r_2111', rxn3, 10, 0.9, 'fseof_beta-carotene_glc.tab');
targets = FSEOF(model, 'r_2111', rxn4, 10, 0.9, 'fseof_gamma-carotene_glc.tab');

%% Perform FSEOF for linolenic acid and TAG on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [0, -1]);

targets = FSEOF(model, 'r_2111', rxn1, 10, 0.9, 'fseof_torularhodin_xyl.tab');
targets = FSEOF(model, 'r_2111', rxn2, 10, 0.9, 'fseof_torulene_xyl.tab');
targets = FSEOF(model, 'r_2111', rxn3, 10, 0.9, 'fseof_beta-carotene_xyl.tab');
targets = FSEOF(model, 'r_2111', rxn4, 10, 0.9, 'fseof_gamma-carotene_xyl.tab');
