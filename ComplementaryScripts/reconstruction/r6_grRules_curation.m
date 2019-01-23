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

%% Curate unlogical grRules
model=changeGeneAssoc(model,'r_0250','RHTO_04703 or RHTO_05549 or RHTO_06321',true);
model=changeGeneAssoc(model,'r_0362','RHTO_02130 or RHTO_02306 or RHTO_07144',true);
model=changeGeneAssoc(model,'r_0658','RHTO_01289 or RHTO_01290 or RHTO_06717',true);
model=changeGeneAssoc(model,{'r_0886','r_0887'},'RHTO_00494',true);
model=changeGeneAssoc(model,'r_0906','RHTO_07357',true);
model=changeGeneAssoc(model,'r_0916','RHTO_02591 or RHTO_04328',true);
model=changeGeneAssoc(model,'r_0961','(RHTO_07250 and RHTO_01852 and RHTO_07893 and RHTO_01754 and RHTO_03543) or RHTO_03059',true);
model=changeGeneAssoc(model,'r_0001','(RHTO_06352 and RHTO_05208) or (RHTO_05208 and RHTO_02645)',true);
model=changeGeneAssoc(model,'r_0002','(RHTO_02645 and RHTO_05208)',true);
model=changeGeneAssoc(model,'r_0004','(RHTO_05208 and RHTO_00251)',true);
model=changeGeneAssoc(model,'r_0437','(RHTO_05208 and RHTO_06193)',true);
model=changeGeneAssoc(model,'r_0550','(RHTO_05117 and RHTO_03771)',true);
model=changeGeneAssoc(model,'r_0552','(RHTO_03771 and RHTO_06286)',true);
model=changeGeneAssoc(model,'r_0883','(RHTO_03771 and RHTO_06542)',true);
model=changeGeneAssoc(model,'r_1021','(RHTO_00723 and RHTO_05714 and RHTO_00534 and RHTO_06068)',true);
model=changeGeneAssoc(model,'r_0510','(RHTO_04065 and RHTO_05749)',true);

%% Remove 'sce' from subsystems
model.subSystems=cellfun(@(x) regexprep(x,'sce[0-9]+ +',''),model.subSystems, 'UniformOutput', 0);

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
