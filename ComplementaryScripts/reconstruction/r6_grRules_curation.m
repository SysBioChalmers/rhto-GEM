%% getModelFromHomology left some OLD_sce genes that it could not find
% orthologs for. Additionally, the reactions that were added by MENECO and
% manually added for coenzyme A are still annotated with the Sce gene (not
% prefixed by OLD_sce_!) Try to find the responsible R. toruloides genes.
load('../../scrap/model_r5.mat');

% All Sce genes have a Y in the name, while Rhto genes do not.
rxnIdx=strfind(model.grRules,'Y');
rxnIdx=~cellfun('isempty',rxnIdx); % Get reaction indexes
out=cell(length(find(rxnIdx)),3);
out(:,1)=model.rxns(rxnIdx);
out(:,2)=model.rxnNames(rxnIdx);
out(:,3)=model.grRules(rxnIdx);

% From this list, through manual curation define the following new grRules
model=changeGeneAssoc(model,'r_0438','(COX1 and COX2 and COX3 and RHTO_00755 and RHTO_05208 and RHTO_04577 and RHTO_01605 and RHTO_03666 and RHTO_01854 and RHTO_06415 and RHTO_06298 and RHTO_04910) or (COX1 and COX2 and COX3 and RHTO_00755 and RHTO_05208 and RHTO_04577 and RHTO_01605 and RHTO_03666 and RHTO_06415 and RHTO_06298 and RHTO_04910 and RHTO_01854) or (COX1 and COX2 and COX3 and RHTO_00755 and RHTO_04577 and RHTO_01605 and RHTO_03666 and RHTO_01854 and RHTO_05208 and RHTO_06415 and RHTO_06298 and RHTO_04910) or (COX1 and COX2 and COX3 and RHTO_00755 and RHTO_04577 and RHTO_01605 and RHTO_03666 and RHTO_05208 and RHTO_06415 and RHTO_06298 and RHTO_04910 and RHTO_01854)',true);
model=changeGeneAssoc(model,'r_1021','(RHTO_00723 and RHTO_05714 and RHTO_00534 and RHTO_06068) or (RHTO_00723 and RHTO_00534 and RHTO_05714 and RHTO_06068)',true);
% UTR4 (YEL038W in S.cerevisiae) seems to have no homologue in Rhto
model=changeGeneAssoc(model,'r_0013','RHTO_05673',true);
% VHS3 (YOR054C in S.cerevisiae) seems to have no homologue in Rhto
model=changeGeneAssoc(model,'r_0906','RHTO_07357',true);
model=changeGeneAssoc(model,'r_1027','RHTO_02113 and RHTO_06769',true);
% YDC1 and YPC1 (YBL078W and YBR183W in S.cerevisiae) seem to have no homologue in Rhto
model=changeGeneAssoc(model,{'r_0340','r_0342'},'',true);
% xylulokinase gene unknown
model=changeGeneAssoc(model,'r_1094','',true);
%model=removeGenes(model,{'YPL087W','YBR183W','YGR194C'},false,false,true);
model=deleteUnusedGenes(model);

%% Remove unconnected non-gene associated reactions
subGraphs=getAllSubGraphs(model);

% Find which reactions have no gene associated
rxnToRemove=[];
for i=1:size(subGraphs,2)
    metIdx = subGraphs(:,i);
    rxnIdx = model.S(metIdx,:);
    [~,col,~] = find(rxnIdx);
    col = unique(col);
    grRules=model.grRules(col);
    if isempty(grRules{1})
        rxnToRemove = [rxnToRemove; col];
    end
end

model = removeReactions(model,rxnToRemove,true,true,true);

save('../../scrap/model_r6.mat','model');
cd('..'); newCommit(model); cd('reconstruction')
