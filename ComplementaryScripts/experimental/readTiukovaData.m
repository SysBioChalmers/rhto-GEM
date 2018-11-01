%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readTiukovaData(i)
%
% Eduard Kerkhoven. Last update: 2018-07-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readTiukovaData(i)

if nargin<1
    i=1;
end
%Lipid data:
fid = fopen('../../ComplementaryData/data/lipidData_Tiukova2018.csv');
lipidData = textscan(fid,['%s %s %s %f32 %s' repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.lipidData.metAbbrev = lipidData{1};
data.lipidData.metNames  = lipidData{2};
data.lipidData.formulas  = lipidData{3};
data.lipidData.MW        = lipidData{4};
data.lipidData.metIDs    = lipidData{5};
data.lipidData.abundance = lipidData{4+2*i};
data.lipidData.std       = lipidData{5+2*i};
fclose(fid);

%Chain data:
fid = fopen('../../ComplementaryData/data/chainData_Tiukova2018.csv');
chainData = textscan(fid,['%s %s %s %f32 %s' repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.chainData.chain     = chainData{1};
data.chainData.FA        = chainData{2};
data.chainData.formulas  = chainData{3};
data.chainData.MW        = chainData{4};
data.chainData.metIDs    = chainData{5};
data.chainData.abundance = chainData{4+2*i};
data.chainData.std       = chainData{5+2*i};
fclose(fid);

% %Other composition data:
% fid = fopen('compData_Lahtvee2017.csv');
% otherData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,9])],'Delimiter',',','HeaderLines',1);
% data.otherData.metIDs    = otherData{2};
% data.otherData.abundance = otherData{2+i};
% fclose(fid);

%Flux data:
fid = fopen('../../ComplementaryData/data/fluxData_Tiukova2018.csv');
fluxData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.fluxData.rxnIDs   = fluxData{2};
data.fluxData.averages = fluxData{2*i+1};
data.fluxData.stdevs   = fluxData{2*i+2};
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%