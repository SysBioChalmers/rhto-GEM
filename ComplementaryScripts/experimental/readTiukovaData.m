%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readTiukovaData(i)
%
%   Reads data from Tiukova et al. data on lipid classes, lipid chain
%   lengths and exchange fluxes, as provided in ComplementaryData. If no
%   sample index is specified, the first sample is selected.
%
% 2019-01-21    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readTiukovaData(model,i)

if nargin<2
    i=1;
end
%Lipid data:
fid = fopen('../../ComplementaryData/data/lipidData_Tiukova2018.csv');
lipidData = textscan(fid,['%s %s %s %f32 %s' repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.lipidData.metAbbrev = lipidData{1};
data.lipidData.metNames  = lipidData{2};
data.lipidData.formulas  = lipidData{3};
data.lipidData.MW        = lipidData{4};
data.lipidData.comp      = lipidData{5};
data.lipidData.abundance = lipidData{4+2*i};
data.lipidData.std       = lipidData{5+2*i};
for j=1:length(data.lipidData.metNames)
    metIds = model.mets(getIndexes(model,[data.lipidData.metNames{j} ...
        '[' data.lipidData.comp{j} ']'],'metscomps'));
    data.lipidData.metIds(j,1) = metIds;
end
fclose(fid);

%Chain data:
fid = fopen('../../ComplementaryData/data/chainData_Tiukova2018.csv');
chainData = textscan(fid,['%s %s %s %f32 %s' repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.chainData.chain     = chainData{1};
data.chainData.FA        = chainData{2};
data.chainData.formulas  = chainData{3};
data.chainData.MW        = chainData{4};
data.chainData.comp      = chainData{5};
data.chainData.abundance = chainData{4+2*i};
data.chainData.std       = chainData{5+2*i};
for j=1:length(data.chainData.chain)
    metIds = model.mets(getIndexes(model,['C' data.chainData.chain{j} ...
        ' chain[' data.chainData.comp{j} ']'],'metscomps'));
    data.chainData.metIds(j,1) = metIds;
end
fclose(fid);

%Flux data:
fid = fopen('../../ComplementaryData/data/fluxData_Tiukova2018.csv');
fluxData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.fluxData.rxnIDs   = fluxData{2};
data.fluxData.averages = fluxData{i+2};%fluxData{2*i+1};
%data.fluxData.stdevs   = fluxData{2*i+2};
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%