function data = readLopesData()
fid = fopen('../../ComplementaryData/data/fluxData_Lopes2019.csv');
types = textscan(fid,repmat('%s',[1,12]),1,'Delimiter',',');
dat = textscan(fid,['%s' repmat(' %f32 ',[1,11])],'Delimiter',',','HeaderLines',1);
fclose(fid);
data.types      = horzcat(types{2:end});
data.sample     = dat{1};
data.rates      = cat(2,dat{[2:10]});
data.rates(isnan(data.rates)) = 0;
data.lipids     = dat{11};
data.protein    = dat{12};
end
