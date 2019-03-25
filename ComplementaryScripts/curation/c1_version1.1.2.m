%% Various curations to move from version 1.1.1 to version 1.1.2
model = importModel('../../ModelFiles/xml/rhto.xml');

% Set xylose exchange to reversible, to allow xylose uptake
model = setParam(model,'rev','r_1718',1); 

save('../../scrap/model_c1.mat','model');
cd('..'); newCommit(model); cd('curation')