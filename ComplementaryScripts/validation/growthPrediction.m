%% This script compares experimentally measured with model predicted growth rates

model       = importModel('../../ModelFiles/xml/rhto.xml');
model       = setParam(model,'rev',{'r_1808','r_1718'},1); %glycerol, xylose
model       = setParam(model,'eq',{'r_1808','r_1718','r_1714'},0);
model       = setParam(model,'obj','r_2111',1);

% Load file
fid         = fopen('../../ComplementaryData/validation/bioreactor_growth.csv');
fluxData    = textscan(fid,'%f32 %f32 %s %s','Delimiter',',','HeaderLines',1);
growth      = fluxData{1};
rate        = fluxData{2};
source      = fluxData{3};
reference   = fluxData{4};
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

[~,~,ic] = unique(reference);

cols = [228,26,28; 55,126,184; 77,175,74; 152,78,163; 255,127,0; 255,255,51; 166,86,40; 247,129,191; 153,153,153];
cols = cols/255;
for i=1:9
    plot(growth(ic == i), out(ic == i), 'o', 'LineWidth', 2.5,...
        'Color', cols(i,:));
    hold on;
end
plot([0 0.3],[0 0.3],':k')
hold off
title('Model-predicted growth rate')
xlabel('Measured growth, h^-^1')
ylabel('Predicted growth, h^-^1')
legend(regexprep(unique(reference),'(.*) et al.* ([0-9]{4}).*','$1 $2'),...
    'Location','eastoutside')%,'NumColumns',2)
set(gca,'FontSize',12) % Creates an axes and sets its FontSize to 18
print('validation.pdf','-dpdf')

