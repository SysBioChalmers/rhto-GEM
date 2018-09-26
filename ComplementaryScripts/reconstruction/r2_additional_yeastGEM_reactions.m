%% Copy pseudoreactions
rxns=modelSce.rxns(contains(modelSce.rxnNames,'pseudoreaction'));
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'Modeling reaction',1); % Add reactions and metabolites

model=addRxnsGenesMets(model,modelSce,'r_4046',false,...
    'Modeling reaction',1); % Add reactions and metabolites

rxns=modelSce.rxns(contains(modelSce.rxnNames,'SLIME'));
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'SLIME reaction',1); % Add reactions and metabolites

%% Add all exchange rxns
% Tthese were not gene annotated, and therefore not added in draft.
% Might not require all exchange rxns, but easier to remove unconnected ones later.
rxns=getExchangeRxns(modelSce);
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'Modelling reaction',1); % Add reactions and metabolites

%% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx=find(cellfun(@isempty,modelSce.grRules)); % Which rxns have no genes
rxnIdx=find(getTransportRxns(modelSce));
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns=modelSce.rxns(rxnIdx); % Obtain reaction IDs
model=addRxnsGenesMets(model,modelSce,rxns,false,...
    'Modeling reaction required for intercellular transport, gene unknown',1); % Add reactions and metabolites

%% Add ATP synthase and ox phosp, problematic to gap-fill
rxns={'r_0226','r_0438','r_0439'};
model=addRxnsGenesMets(model,modelSce,rxns,true,...
    'Energy metabolism',1); % Add reactions and metabolites

model=setParam(model,'eq','r_1714',-1);
model=setParam(model,'obj','r_2111',1);

%% Some offline curation identified additional reactions
rxns={'r_0001';'r_0002';'r_0004';'r_0013';'r_0019';'r_0021';'r_0022';'r_0041';'r_0063';'r_0064';'r_0066';'r_0067';'r_0109';'r_0145';'r_0146';'r_0195';'r_0236';'r_0237';'r_0320';'r_0340';'r_0342';'r_0344';'r_0346';'r_0347';'r_0350';'r_0357';'r_0362';'r_0369';'r_0370';'r_0437';'r_0453';'r_0482';'r_0499';'r_0501';'r_0504';'r_0506';'r_0507';'r_0508';'r_0509';'r_0510';'r_0532';'r_0537';'r_0551';'r_0656';'r_0657';'r_0678';'r_0692';'r_0736';'r_0760';'r_0767';'r_0844';'r_0851';'r_0883';'r_0904';'r_0913';'r_0917';'r_0920';'r_0937';'r_0938';'r_0942';'r_0943';'r_0963';'r_0970';'r_0973';'r_0985';'r_0993';'r_1004';'r_1005';'r_1021';'r_1026';'r_1027';'r_1030';'r_1032';'r_1033';'r_1051';'r_1063';'r_1083';'r_1102';'r_1108';'r_1109';'r_1119';'r_1132';'r_1133';'r_1136';'r_1146';'r_1147';'r_1161';'r_1173';'r_1176';'r_1183';'r_1186';'r_1190';'r_1191';'r_1196';'r_1199';'r_1207';'r_1217';'r_1225';'r_1250';'r_1253';'r_1254';'r_1267';'r_1270';'r_1278';'r_2034';'r_2117';'r_2118';'r_2119';'r_2171';'r_2172';'r_2173';'r_2178';'r_2179';'r_2180';'r_2181';'r_2195';'r_2196';'r_2197';'r_2198';'r_2199';'r_2201';'r_2202';'r_2203';'r_2204';'r_2205';'r_2214';'r_2215';'r_2217';'r_2218';'r_2219';'r_2220';'r_2221';'r_2222';'r_2223';'r_2224';'r_2225';'r_2226';'r_2227';'r_2228';'r_2232';'r_2233';'r_2332';'r_2334';'r_3312';'r_3313';'r_3314';'r_3315';'r_3324';'r_3325';'r_3326';'r_3327';'r_3328';'r_3329';'r_3330';'r_3331';'r_3428';'r_3429';'r_3430';'r_3431';'r_3432';'r_3433';'r_3434';'r_3435';'r_3436';'r_3437';'r_3438';'r_3439';'r_3440';'r_3441';'r_3442';'r_3443';'r_3444';'r_3445';'r_3446';'r_3447';'r_3448';'r_3449';'r_3450';'r_3451';'r_3452';'r_3453';'r_3454';'r_3455';'r_3456';'r_3457';'r_3458';'r_3459';'r_3460';'r_3461';'r_3462';'r_3463';'r_3464';'r_3465';'r_3466';'r_3467';'r_3468';'r_3469';'r_3470';'r_3471';'r_3472';'r_3473';'r_3474';'r_3475';'r_3476';'r_3477';'r_3478';'r_3479';'r_3480';'r_3481';'r_3482';'r_3483';'r_3484';'r_3485';'r_3486';'r_3487';'r_3488';'r_3489';'r_3490';'r_3491';'r_3492';'r_3493';'r_3494';'r_3495';'r_3496';'r_3497';'r_3498';'r_3499';'r_3500';'r_3501';'r_3502';'r_3503';'r_3504';'r_3505';'r_3506';'r_3507'};
grRules={'(RHTO_06352 and RHTO_05208) or (RHTO_02645 and RHTO_05208)';'RHTO_02645 and RHTO_05208';'RHTO_05208 and RHTO_02605';'YEL038W and RHTO_05673';'RHTO_01031';'RHTO_06687 and RHTO_03730 and RHTO_05146 and RHTO_00339 and RHTO_06918 and RHTO_00813 and RHTO_01107';'RHTO_06687 and RHTO_03730 and RHTO_05146 and RHTO_00339 and RHTO_06918 and RHTO_00813 and RHTO_01107';'RHTO_07838';'RHTO_06605';'RHTO_00098';'RHTO_01811';'RHTO_04877';'RHTO_02161 and RHTO_02004';'RHTO_05506';'RHTO_02358';'RHTO_01529 or RHTO_06117';'RHTO_07171';'RHTO_07171';'RHTO_05579';'YPL087W';'YBR183W or YPL087W';'RHTO_02268';'RHTO_03440';'RHTO_00942';'RHTO_00942';'RHTO_02997';'RHTO_02306 or RHTO_07144';'RHTO_02861';'RHTO_01931 or RHTO_02076';'RHTO_05208 and RHTO_06193';'RHTO_06990';'RHTO_04201 and RHTO_03104';'RHTO_01888';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'(RHTO_04065 and RHTO_05749) or (RHTO_04065 and RHTO_05749) or (RHTO_05749 and RHTO_04065) or (RHTO_05749 and RHTO_04065)';'RHTO_06687 and RHTO_03730 and RHTO_05146 and RHTO_00339 and RHTO_06918 and RHTO_00813 and RHTO_01107';'RHTO_06924';'RHTO_06930 and RHTO_00443';'RHTO_02861';'RHTO_02861';'RHTO_01494 and RHTO_06911';'RHTO_04169';'RHTO_02122';'RHTO_04876';'RHTO_00803';'RHTO_01116';'RHTO_00143 or RHTO_01136 or RHTO_07034';'(RHTO_02963 and RHTO_06542) or (RHTO_01875 and RHTO_06542)';'RHTO_02073';'RHTO_02564';'RHTO_03868';'RHTO_04142';'RHTO_05237';'RHTO_03474';'RHTO_00658';'RHTO_04784 and RHTO_05934 and RHTO_05997';'RHTO_06687 and RHTO_03730 and RHTO_05146 and RHTO_00339 and RHTO_06918 and RHTO_00813 and RHTO_01107';'RHTO_02963';'RHTO_02963';'RHTO_06687 and RHTO_03730 and RHTO_05146 and RHTO_00339 and RHTO_06918 and RHTO_00813 and RHTO_01107';'(RHTO_01432 and RHTO_06834) or (RHTO_01432 and RHTO_06834)';'RHTO_04435';'RHTO_04435';'(RHTO_00723 and RHTO_05714 and RHTO_00534 and RHTO_06068) or (RHTO_00723 and RHTO_00534 and RHTO_05714 and RHTO_06068)';'RHTO_02552';'YFR030W and RHTO_02113';'RHTO_02808 and RHTO_01709 and RHTO_07893 and RHTO_06925';'RHTO_07522';'RHTO_07522';'RHTO_01529 or RHTO_01529';'RHTO_00143';'RHTO_03383';'RHTO_05452';'RHTO_07997 or RHTO_02411';'RHTO_06628';'RHTO_04753';'RHTO_06628';'RHTO_07997 or RHTO_02411';'RHTO_06945';'RHTO_07037';'RHTO_07037';'RHTO_07037';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_07997 or RHTO_02411';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_00398';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_06945';'RHTO_00393 or RHTO_00398 or RHTO_01079 or RHTO_01344 or RHTO_01346 or RHTO_02350 or RHTO_07825 or RHTO_07915';'RHTO_07037';'RHTO_05974';'RHTO_02411';'RHTO_06945';'RHTO_05047';'RHTO_06782';'RHTO_07037';'(RHTO_02579 and RHTO_00085) or (RHTO_02579 and RHTO_00085)';'RHTO_07034';'RHTO_07034';'RHTO_07034';'RHTO_03657';'RHTO_03657';'RHTO_03657';'RHTO_05958';'RHTO_05958';'RHTO_05958';'RHTO_05958';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_04350';'RHTO_07775';'RHTO_07775';'RHTO_07775';'RHTO_07775';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_06195 and RHTO_04105';'RHTO_01116';'RHTO_01116';'RHTO_05332';'RHTO_05332';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_03931';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679';'RHTO_02119 and RHTO_01679'};

model=addRxnsGenesMets(model,modelSce,rxns,grRules,'Identified through homology',2);

cd('..')
newCommit(model);
cd('reconstruction')