if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% getModelFromHomology left some OLD_sce genes that it could not find
% orthologs for. Additionally, the reactions that were added by MENECO and
% manually added for coenzyme A are still annotated with the Sce gene (not
% prefixed by OLD_sce_!) Try to find the responsible R. toruloides genes.
load([root '/scrap/model_r6.mat']);

% All Sce genes have a Y in the name, while Rhto genes do not.
rxnIdx      = strfind(model.grRules,'Y');
rxnIdx      = ~cellfun('isempty',rxnIdx); % Get reaction indexes
out         = cell(length(find(rxnIdx)),3);
out(:,1)    = model.rxns(rxnIdx);
out(:,2)    = model.rxnNames(rxnIdx);
out(:,3)    = model.grRules(rxnIdx);

% From this list, through manual curation define the following new grRules
fid         = fopen([data '/reconstruction/updateGrRules.txt']);
loadedData  = textscan(fid,'%s %s','delimiter','\t'); fclose(fid);
rxns        = loadedData{1};
grRules     = loadedData{2};

model = changeGrRules(model,rxns,grRules,true);

% Correct faulty grRules where the same complex is representated multiple
% times
for n = 1:length(model.grRules)
    if any(model.grRules{n})
        noAnd = strfind(model.grRules(n),'and');
        noAnd = any(vertcat(noAnd{:})); % Give 0 if no 'and' is present.
        if noAnd == 0
            geneList = transpose(cell(unique(regexp(model.grRules{n},'[)(]*|( and )*|( or )*','split'))));
            geneList = regexprep(geneList,'[(*)*]','');
            if length(geneList) == 1
                newgrRule = geneList;
            else
                newgrRule = geneList{1};
                for k = 2:length(geneList)
                    newgrRule = [newgrRule ' or ' geneList{k}];
                end
            end
            model.grRules(n) = cellstr(newgrRule);
        end
    end
end
exportToExcelFormat(model,'rhto.xlsx');
tmp = haveFlux(model);
modelCb=ravenCobraWrapper(model);
model = removeUnusedGenes(model);

save([root '/scrap/model_r7.mat'],'model');
%cd('..'); newCommit(model); cd('reconstruction')
