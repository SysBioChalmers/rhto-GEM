clear;clc;if ~exist('scripts') | ~endsWith(scripts,'ComplementaryScripts'); run('../../init_rhtoGEM.m'); end
%% Theoretical simulations
% Take rhtoGEM model as it is
cd([scripts '/experimental/']);
model = importModel([root '\ModelFiles\xml\rhto.xml']);
model.rev(getIndexes(model,{'r_1634','r_1718','r_1808'},'rxns')) = 1;

% Allow excretion of lipid chains and backbones. This allows for optimizing
% production of lipid pseudometabolite
chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
tmp = setParam(model,'ub',[chainExIdx,backbExIdx],1000);
tmp = setParam(tmp,'obj','r_4062',1); % Maximize lipid pseudometabolite production
tmp = setParam(tmp,'eq','r_4046',0); % Set NGAM to 0

% Get optimal solutions for each carbon source
tmp = setParam(tmp,'lb',{'r_1634','r_1808','r_1718','r_1714'},[-1/2,0,0,0]);
solAce = solveLP(tmp,1);

tmp = setParam(tmp,'lb',{'r_1634','r_1808','r_1718','r_1714'},[0,-1/3,0,0]);
solGly = solveLP(tmp,1);

tmp = setParam(tmp,'lb',{'r_1634','r_1808','r_1718','r_1714'},[0,0,-1/5,0]);
solXyl = solveLP(tmp,1);

tmp = setParam(tmp,'eq','t_0182',0); % Disable phosphoketolase
solXylnoPK = solveLP(tmp,1);

% Print maximum lipid yields
fprintf(['Lipid pseudometabolite produced per C from each carbon source:\n' ...
    'Acetate: ' num2str(-solAce.f) ' Glycerol: ' num2str(-solGly.f) ...
    ' Xylose: ' num2str(-solXyl.f) ' Xylose (no PK-PTA): ' num2str(-solXylnoPK.f) '\n'])

% Write file with all fluxes
out = [solAce.x, solGly.x, solXyl.x solXylnoPK.x];
out = [model.rxnNames num2cell(out)];
fid = fopen([data '/Lahtvee2019/theor_maxLipid_allFluxes.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["reaction" "Acetate" "Glycerol" "Xylose" ...
    "Xylose (no PK-PTA)"]);
for j=1:length(out)
    fprintf(fid,'%s\t%d\t%d\t%d\t%d\n',out{j,:});
end
fclose(fid);

%% Write file with exchange fluxes
[~,idx] = getExchangeRxns(model,'both');
out = out(idx,:);

% Remove exchange reactions not carrying flux
rmIdx = sum(cell2mat(out(:,2:end)),2) == 0;
out = out(~rmIdx,:);

fid = fopen([data '/Lahtvee2019/theor_maxLipid_exchFluxes.tsv'],'w');
fprintf(fid,['%s' repmat('\t%s',1,4) '\n'],["reaction" "acetate" ...
    "glycerol" "xylose" "xylose (no PK-PTA)"]);
for j=1:length(out)
    fprintf(fid,['%s' repmat('\t%d',1,4) '\n'],out{j,:});
end
fclose(fid);
