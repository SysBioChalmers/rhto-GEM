%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readTiukovaData(i)
%
% Eduard Kerkhoven. Last update: 2018-07-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readTiukovaData(i)

%Lipid data:
fid = fopen('../ComplementaryData/lipidData_Tiukova2018.csv');
lipidData = textscan(fid,[repmat('%s ',[1,3]) repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.lipidData.metAbbrev = lipidData{1};
data.lipidData.metNames  = lipidData{2};
data.lipidData.metIDs    = lipidData{3};
data.lipidData.abundance = lipidData{2+2*i};
data.lipidData.std       = lipidData{3+2*i};
fclose(fid);

%Chain data:
fid = fopen('../ComplementaryData/chainData_Tiukova2018.csv');
chainData = textscan(fid,[repmat('%s ',[1,3]) repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.chainData.metNames  = chainData{1};
data.chainData.formulas  = chainData{2};
data.chainData.metIDs    = chainData{3};
data.chainData.abundance = chainData{2+2*i};
data.chainData.std       = chainData{3+2*i};
fclose(fid);

% %Other composition data:
% fid = fopen('compData_Lahtvee2017.csv');
% otherData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,9])],'Delimiter',',','HeaderLines',1);
% data.otherData.metIDs    = otherData{2};
% data.otherData.abundance = otherData{2+i};
% fclose(fid);

%Flux data:
fid = fopen('../ComplementaryData/fluxData_Tiukova2018.csv');
fluxData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,8])],'Delimiter',',','HeaderLines',1);
data.fluxData.rxnIDs   = fluxData{2};
data.fluxData.averages = fluxData{2*i+1};
data.fluxData.stdevs   = fluxData{2*i+2};
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%