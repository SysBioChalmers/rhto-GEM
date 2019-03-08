%% This script compares experimentally measured with model predicted growth rates

model       = importModel('../../ModelFiles/xml/rhto.xml');
model       = setParam(model,'rev',{'r_1808','r_1718'},1); %glycerol, xylose
model       = setParam(model,'eq',{'r_1808','r_1718','r_1714'},0);
model       = setParam(model,'obj','r_2111',1);

% Load file
fid         = fopen('../../ComplementaryData/data/validationData.csv');
fluxData    = textscan(fid,'%f32 %f32 %s %s','Delimiter',',','HeaderLines',1);
growth      = fluxData{1};
rate        = fluxData{2};
source      = fluxData{3};
fclose(fid);

for i = 1:length(growth)
    if strcmp(source(i),'glucose')
        modelTmp = setParam(model,'lb','r_1714',-rate(i));
    elseif strcmp(source(i),'glycerol')
        modelTmp = setParam(model,'lb','r_1808',-rate(i));
	elseif strcmp(source(i),'xylose')
        modelTmp = setParam(model,'lb','r_1718',-rate(i));
    end
    sol=solveLP(modelTmp);
    out(i) = -sol.f;
end
out=transpose(out);

plot(growth,out,'o')
hold on; plot([0 0.3],[0 0.3])
hold off
