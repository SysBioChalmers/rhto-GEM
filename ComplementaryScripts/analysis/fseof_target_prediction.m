%% This script predicts metabolic engineering targets for increased
% production of triglycerides and linolenic acid.
model = importModel('../ModelFiles/xml/rhto.xml');

% Add exchange reactions for linolenic acid and triglyceride. For the
% triglyceride, we specifically choose 18:0/18:1/18:0-TAG (so called SOS)
% as target. In addition to functioning as representative of triglycerides
% for biofuel production, this SOS-TAG is a major part of cocoa-butter.
model = addExchangeRxns(model,'out',{'st_0107','s_3032'});

%% Perform FSEOF for linolenic acid and TAG on glucose
targets = FSEOF(model, 'r_2111', 'EXC_OUT_st_0107',10,0.9,...
    'fseof_linolenic acid_glc.tab');
targets = FSEOF(model, 'r_2111', 'EXC_OUT_s_3032',10,0.9,...
    'fseof_TAG_glc.tab');

%% Perform FSEOF for linolenic acid and TAG on xylose
model       = setParam(model,'eq','r_1714',0);
model       = setParam(model,'eq','r_1718',-1);
targets = FSEOF(model, 'r_2111', 'EXC_OUT_s_3032',10,0.9,...
    'fseof_TAG_xyl.tab');
targets = FSEOF(model, 'r_2111', 'EXC_OUT_st_0107',10,0.9,...
    'fseof_linolenic acid_xyl.tab');
